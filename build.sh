#!/usr/bin/env bash
set -ex


# Build patched flex, if needed
mkdir -p localinstall
pushd localinstall
localInstall="`pwd`"
if [[ ! -d "flex-2.5.4" ]]; then
	tar -xf ../distrib_extras/flex-2.5.4_patch.tar.gz
fi
pushd flex-2.5.4
./configure --prefix=${localInstall}
make install
popd
popd

# Create build subdirectories
pushd distrib
buildDir="`pwd`"
buildDirs="api_include api_lib bin lib include man/man3 share/doc"
for dir in ${buildDirs}; do
	mkdir -p ${dir}
done
ln -fs sources obj
ln -fs share/etc etc

pushd obj
make WITH_FLEXLM=NOFLEX            \
  ALLIANCE_TOP=${buildDir}         \
  AVERTEC_TOP=${buildDir}          \
  AVERTEC_OS=Linux                 \
  AVERTEC_LICENSE=AVERTEC_DUMMY    \
  AVT_LICENSE_SERVER=house         \
  AVT_LICENSE_FILE=27009@house     \
  AVT_COMPILATION_TYPE=distrib     \
  AVT_DISTRIB_DIR=${buildDir}      \
  PACKAGING_TOP=${localInstall}    \
  LEX=${localInstall}/bin/flex     \
  JAVA_HOME=/usr/lib/jvm/default   \
  CFLAGS="-g -O3" CXXFLAGS=-"-g -O3" STRIP=true \
  SAXON="java -jar ${buildDir}/../distrib_extras/saxon9.jar"
popd
popd


# 'Install'
mkdir -p install
installDir="`pwd`/install"

mkdir -p ${installDir}/bin
mkdir -p ${installDir}/share/tasyag/etc

for conf in avt.slib avttools.dtb Xtas Xyagle trmodel.cfg; do
  cp distrib/share/etc/$conf ${installDir}/share/tasyag/etc
done

for tool in avt_shell avtman xtas xyagle ttvdiff ttvren; do
  cp distrib/bin/${tool} ${installDir}/bin
done

cp -r distrib/share/tcl ${installDir}/share/tasyag

echo "AVERTEC_TOP=`pwd`/install/share/tasyag" > "${installDir}/avt_env.sh"
echo 'PATH=${AVERTEC_TOP}/tcl:${PATH}' >> "${installDir}/avt_env.sh"
echo 'export AVERTEC_TOP' >> "${installDir}/avt_env.sh"
