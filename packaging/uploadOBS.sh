
    version="3.4.6"
     obsDir="../coriolis-obs/home:jpc-lip6/tas-yagle"

 printHelp () {
   echo ""
   echo "  Usage: uploadOBS.sh [--sources] [--commit] [--run]"
   echo ""
   echo "  Options:"
   echo "    [--sources] : Build an archive from the HEAD of the current branch."
   echo "    [--commit]  : Push the files (commit) on the remote builder repository."
   echo "                  This will effectively triggers the rebuild of the packages."
   echo "                  OBS local repository is hardwired to:"
   echo "                      \"${obsDir}\""
   echo "    [--run]     : Perform all actions at once."
   echo ""

 }

 if [ $# -eq 0 ]; then printHelp; fi

    githash=`git log -1 --pretty=format:%h`
  doSources="false"
   doCommit="false"
 badAgument=""
 while [ $# -gt 0 ]; do
   case $1 in
     --sources) doSources="true";;
     --commit)  doCommit="true";;
     --run)     doSources="true"
                doCommit="true";;
     *)         badArgument="$1";;
   esac
   shift
 done
 if [ ! -z "${badArgument}" ]; then
   echo "[ERROR] patchenv.sh: Unknown argument \"${badArgument}\"."
   exit 1
 fi

 echo "Running uploadOBS.sh"
 echo "* Using HEAD githash as release: ${githash}."
 if [ "${doSources}" = "true" ]; then
   echo "* Making source file archive from Git HEAD ..."
   ./packaging/git-archive-all.sh -v --prefix tas-yagle-${version}/ --format tar.gz tas-yagle-${version}.tar.gz
   #git archive --prefix=tas-yagle-${version}/ --format=tar.gz -o tas-yagle-${version}.tar.gz HEAD
 fi

 echo "* Update files in OBS project directory."
 echo "  OBS package directory: \"${obsDir}\"."
 for distribFile in packaging/tas-yagle.spec               \
                    packaging/tas-yagle-rpmlintrc          \
                    packaging/flex-2.5.4-exit.patch        \
                    tas-yagle-${version}.tar.gz            \
		    packaging/tas-yagle.dsc                \
		    packaging/debian.copyright             \
		    packaging/debian.changelog             \
		    packaging/debian.control               \
		    packaging/debian.rules                 \
		    packaging/debian.tas-yagle.install     \
		    packaging/debian.tas-yagle-doc.install \
		    ; do
   if [ ! -f "${distribFile}" ]; then continue; fi
   if [[ "${distribFile}" == packaging* ]]; then
     echo "  - copy ${distribFile}."
     cp ${distribFile} ${obsDir}
   else
     echo "  - move ${distribFile}."
     mv ${distribFile} ${obsDir}
   fi
 done
 
 sed -i "s,^Release: *1,Release:        <CI_CNT>.<B_CNT>.${githash}," ${obsDir}/tas-yagle.spec
 sed -i "s,^%define docGithash .*,%define docGithash ${docGithash},"  ${obsDir}/tas-yagle.spec
 if [ "${doCommit}" = "true" ]; then
   pushd ${obsDir}
   osc commit
   popd
 fi

