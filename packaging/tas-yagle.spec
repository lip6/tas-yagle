
%define _unpackaged_files_terminate_build 1
%define debug_package                     %{nil}

Name:           tas-yagle
Version:        3.4p6
Release:        1
Summary:        Tas/Yagle - Static Timing Analyser
Group:          Applications/Engineering
License:        GPL-2.0-only         
URL:            https://coriolis.lip6.fr/
Source0:        %{name}-%{version}.tar.gz
Patch0:         flex-2.5.4-exit.patch 
#Packager:       Jean-Paul Chaput <Jean-Paul.Chaput@lip6.fr>
#BuildRequires:  texlive-scheme-tetex
BuildRequires:  libstdc++-devel
BuildRequires:  gcc-c++
BuildRequires:  tcsh
BuildRequires:  tcl-devel >= 8.5.3
BuildRequires:  motif-devel
BuildRequires:  libXt-devel
BuildRequires:  libXp-devel
BuildRequires:  libXpm-devel
BuildRequires:  ncurses-devel
BuildRequires:  byacc
BuildRequires:  bison
BuildRequires:  swig >= 1.3.27
BuildRequires:  lua
#BuildRequires:  fop >= 0.95 
#BuildRequires:  saxon9 >= 9.3.0.4
BuildRequires:  libedit-devel >= 2.11
Requires:       libedit >= 2.11

# All RHEL clones.
%if 0%{?rhel} || 0%{?fedora}
%global  build_variant  Linux.el9

BuildRequires:  lsb_release
%endif

%if 0%{?rhel} >= 9 || 0%{?fedora} >= 39
BuildRequires:  java-1.8.0-openjdk-headless
%endif

# All openSUSE.
%if 0%{?is_opensuse}
BuildRequires:  lsb-release
BuildRequires:  java-1_8_0-openjdk-headless
%endif

#opensuse_leap-15.6
%if 0%{?sle_version} == 150600 && 0%{?is_opensuse}
%global  build_variant  Linux.openSUSE

BuildRequires:  termcap
%endif

#opensuse_Tumbleweed
%if 0%{?suse_version} > 1600 && 0%{?is_opensuse}
%global  build_variant  Linux.openSUSE

BuildRequires:  libtermcap2
%endif


%description
STATIC TIMING ANALYSIS
  The  advent   of  semiconductor  fabrication  technologies   now  allows  high
performance in complex integrated circuits.
  With the increasing complexity of these circuits, static timing analysis (STA)
has  revealed  itself  as  the  only  feasible  method  ensuring  that  expected
performances are actually obtained.
  In addition, signal integrity (SI) issues due to crosstalk play a crucial role
in performance and reliability of these  systems, and must be taken into account
during the timing analysis.
  However, performance  achievement not  only lies in  fabrication technologies,
but also  in the way circuits  are designed.  Very high  performance designs are
obtained with semi or full-custom designs techniques.
  The HITAS platform provides advanced STA and SI solutions at transistor level.
It has been  built-up in order to allow engineers to  ensure complete timing and
SI coverage on their digital custom  designs, as well as IP-reuse through timing
abstraction.
  Furthermore,  hierarchy  handling  through  transparent  timing  views  allows
full-chip verification, with virtually no limit of capacity in design size.


%prep
%setup


%build
# Compile the local version of Flex.
 tar zxf ./distrib_extras/flex-2.5.4_patch.tar.gz
 localInstall=`pwd`
 pushd flex-2.5.4
 patch -p1 < %{_sourcedir}/flex-2.5.4-exit.patch
 ./configure --prefix=${localInstall}
 make install
 popd
 PATH=${localInstall}/bin:$PATH
 export PATH

# Setting up the build layout.
# tas/yag uses an old fashioned Makefile.
 pushd distrib
 buildDir=`pwd`
 AVERTEC_OS="Linux"

 buildDirs="api_include api_lib bin lib include man/man3 share/doc"
 for dir in ${buildDirs}; do
   mkdir -p ${dir}
 done
 ln -s sources obj
 ln -s share/etc etc

 makefileEnv=""

 pushd obj
 make WITH_FLEXLM=NOFLEX               \
      BUILD_VARIANT=%{build_variant}   \
      LIB_SUFFIX=64                    \
      LIB_SUFFIX_=_64                  \
      ALLIANCE_TOP=${buildDir}         \
      AVERTEC_TOP=${buildDir}          \
      AVERTEC_OS=$AVERTEC_OS           \
      AVERTEC_LICENSE=AVERTEC_DUMMY    \
      AVT_LICENSE_SERVER=house         \
      AVT_LICENSE_FILE=27009@house     \
      AVT_COMPILATION_TYPE=distrib     \
      AVT_DISTRIB_DIR=${buildDir}      \
      PACKAGING_TOP=${localInstall}

 popd

# Build documentation.
 PATH=`pwd`/bin:$PATH
 pushd docxml2
%ifarch x86_64 aarch64
# wrap_nolicense makes core dump under 64 bits.
# Uses the precompiled version of the doc (from RHEL6 32bits).
 tar jxf docxml2-compiled.tar.bz2
%else
 make AVERTEC_TOP=${buildDir}/share
%endif
 popd

 popd


%install
 rm -rf %{buildroot}
 mkdir -p %{buildroot}%{_sysconfdir}/profile.d
 mkdir -p %{buildroot}%{_bindir}
 mkdir -p %{buildroot}%{_datadir}/tasyag/etc
 mkdir -p %{buildroot}%{_mandir}/man3

 for conf in avt.slib avttools.dtb Xtas Xyagle trmodel.cfg; do
   cp distrib/share/etc/$conf %{buildroot}%{_datadir}/tasyag/etc
 done

 for tool in avt_shell avtman xtas xyagle ttvdiff ttvren; do
   cp distrib/bin/$tool %{buildroot}%{_bindir}
 done

 cp distrib/share/etc/avt_env.sh %{buildroot}%{_sysconfdir}/profile.d
 cp -r distrib/share/tcl %{buildroot}%{_datadir}/tasyag
 cp distrib/man/man3/* %{buildroot}%{_mandir}/man3
 mv distrib/docxml2/compiled/{docavertec.html,docpdf,dochtml} .
 mv distrib/share/tutorials .
 find . -name '*.gif' | xargs chmod a-x


%files
%defattr(-,root,root,-)
%doc docavertec.html dochtml docpdf tutorials LICENSE.rst
%dir %{_sysconfdir}/profile.d
%dir %{_bindir}
%dir %{_mandir}/man3
%dir %{_datadir}/tasyag
%dir %{_datadir}/tasyag/etc
%dir %{_datadir}/tasyag/tcl
%config %{_sysconfdir}/profile.d/avt_env.sh
%{_bindir}/*
%config %{_datadir}/tasyag/etc/*
%{_datadir}/tasyag/tcl/*
%{_mandir}/man3/*



%changelog
* Fri Dec 27 2024 Jean-Paul Chaput <Jean-Paul.Chaput@lip6.fr> - 3.5p5.1
- Ported to openSUSE Build System (OBS).

* Tue Apr 17 2012 Jean-Paul Chaput <Jean-Paul.Chaput@lip6.fr> - 3.5p5.1
- Initial packaging.
