Format:            1.0
Source:            tas-yagle
Binary:            tas-yagle, tas-yagle-dev, tas-yagle-doc
Architecture:      any
Version:           3.4.6
Maintainer:        Jean-Paul.Chaput <Jean-Paul.Chaput@lip6.fr>
Homepage:          https://coriolis.lip6.fr/
Standards-Version: 3.4.6
Build-Depends:     debhelper-compat (= 13),
                   tcsh,
		   build-essential,
		   pkg-config,
                   bison,
                   gawk,
                   libedit-dev,
		   libstdc++6,
		   libboost-all-dev,
		   libbz2-dev,
		   libxml2-dev,
		   libxpm-dev,
		   libmotif-common,
		   libmotif-dev,
		   libxm4,
		   tcl-dev,
		   swig,
		   texlive-latex-recommended,
		   graphviz,
		   xpdf,
                   openjdk-17-jre-headless,
Package-List:
 tas-yagle     deb Science/Electronics optional arch=any
 tas-yagle-doc deb Science/Electronics optional arch=any
DEBTRANSFORM-RELEASE:   1
DEBTRANSFORM-TAR:       tas-yagle-3.4.6.tar.gz
