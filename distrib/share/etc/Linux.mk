# -*- Mode: Makefile -*-
#
####---------------------------------------------------------###
# description  : Alliance include file for Makefiles
# architecture : Linux avertec.lip6.fr 2.2.5-15 #2 ven oct 1 10:39:45 CEST 1999 i686 unknown
# date         : Mon Oct  4 17:38:18 CEST 1999
# file         : Linux.mk
#

# The variables $ALLIANCE_* are set by
# alc_env.[c]sh script or libraries.mk

$(info Processing Linux.mk)

ifeq ($(BUILD_VARIANT),)
  $(info Performing OS recognition based on uname.)
  UNAME_S = $(shell uname -s)
  UNAME_R = $(shell uname -r)
  UNAME_M = $(shell uname -m)
  
  LIB_SUFFIX  = ""
  LIB_SUFFIX_ = ""
  ifeq ($(UNAME_M),x86_64)
    LIB_SUFFIX  = 64
    LIB_SUFFIX_ = _64
  endif
  
  BUILD_VARIANT = Linux
  ifeq ($(UNAME_S),Linux)
    $(info Configuring for a Linux system "$(UNAME_R)".)
    ifneq ($(findstring .el6.,$(UNAME_R)),)
      BUILD_VARIANT    = Linux.slsoc6x
      $(info Configuring for RHEL 6 (and clones).)
    endif
    ifneq ($(findstring .slsoc6.,$(UNAME_R)),)
      BUILD_VARIANT    = Linux.slsoc6x
      $(info Configuring for LIP6 / RHEL 6 (and clones).)
    endif
    ifneq ($(findstring .el7.,$(UNAME_R)),)
      BUILD_VARIANT    = Linux.el7
      $(info Configuring for RHEL 7 (and clones).)
    endif
    ifneq ($(findstring .el8,$(UNAME_R)),)
      BUILD_VARIANT    = Linux.el8
      $(info Configuring for RHEL 8 (and clones).)
    endif
    ifneq ($(findstring .el9,$(UNAME_R)),)
      BUILD_VARIANT    = Linux.el9
      $(info Configuring for RHEL 9 (and clones).)
    endif
    ifneq ($(findstring ubuntu.,$(UNAME_R)),)
      BUILD_VARIANT    = Linux.ubuntu
      $(info Configuring for Ubuntu.)
    endif
    ifneq ($(findstring openSUSE.,$(shell lsb_release -i)),)
      BUILD_VARIANT    = Linux.openSUSE
      $(info Configuring for openSUSE.)
    endif
  endif
else
  $(info Using user supplied values:)
  $(info * BUILD_VARIANT=${BUILD_VARIANT})
  $(info * LIB_SUFFIX   =${LIB_SUFFIX})
  $(info * LIB_SUFFIX_  =${LIB_SUFFIX_})
endif


PROGRAM_SUFFIX   =

GNU_LIB          = /usr/lib
GNU_INCLUDE      = /usr/include

X11_LIB          = /usr/lib 
X11_INCLUDE      = /usr/include

MOTIF_LIB        = /usr/lib64 -L/usr/lib
MOTIF_INCLUDE    = /usr/include

XPM_LIB          = /usr/lib
XPM_INCLUDE      = /usr/include

SHELL            = /bin/sh
CSH              = /bin/csh
CP               = /bin/cp
CAT              = /bin/cat
MV               = /bin/mv
RM               = /bin/rm
MKDIR            = /bin/mkdir
FIND             = /usr/bin/find
SED              = /bin/sed
ifeq ($(BUILD_VARIANT),Linux.ubuntu)
AWK              = /usr/bin/gawk
else			 
AWK              = /bin/awk
endif			 
TR               = /usr/bin/tr
TOUCH            = /bin/touch
ECHO             = /bin/echo
STRIP            = /usr/bin/strip
RANLIB           = /usr/bin/ranlib

MAKE             = /usr/bin/make
MAKEFLAGS        = 

CC               = /usr/bin/gcc -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -DHAVE_UNISTD_H
SCC              = /usr/bin/gcc -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -DHAVE_UNISTD_H
CPLUSPLUS        = /usr/bin/g++ -DHAVE_UNISTD_H
CFLAGS           =
CPPFLAGS         =

CC              += -I/usr/include/tcl
SCC             += -I/usr/include/tcl
CPLUSPLUS       += -I/usr/include/tcl

ifeq ($(PACKAGING_TOP),)
  CC            += -I${HOME}/softs/$(BUILD_VARIANT)$(LIB_SUFFIX_)/install/include
  SCC           += -I${HOME}/softs/$(BUILD_VARIANT)$(LIB_SUFFIX_)/install/include
  CPLUSPLUS     += -I${HOME}/softs/$(BUILD_VARIANT)$(LIB_SUFFIX_)/install/include
else
  CC            += -I${PACKAGING_TOP}/include
  SCC           += -I${PACKAGING_TOP}/include
  CPLUSPLUS     += -I${PACKAGING_TOP}/include
endif

ifeq ($(shell uname -m),x86_64)
  AVT_COMPILATION_64BIT = yes
endif

OPTIM            = -O3

ENABLE_STATIC    = -Xlinker -Bstatic
DISABLE_STATIC   = -Xlinker -Bdynamic

PURIFY           = purify

YACC             = /usr/bin/bison
YACCFLAGS        = -y 

#LEX             = flex
ifeq ($(PACKAGING_TOP),)
LEX              = ${HOME}/softs/$(BUILD_VARIANT)$(LIB_SUFFIX_)/install/bin/flex
else
LEX              = ${PACKAGING_TOP}/bin/flex
endif
LEXFLAGS         =

AR               = /usr/bin/ar
ARFLAGS          = rv

SWIG             = /usr/bin/swig

WHOLE            = -Xlinker --whole-archive
NOWHOLE          = -Xlinker --no-whole-archive

ifeq ($(BUILD_VARIANT),Linux.ubuntu)
TCL_L            = $(shell pkg-config --libs tcl8.6)
TCL_INCLUDES     = -I/usr/include/tcl8.6
$(info Debian, looking for tcl8.6)
$(info -> $(TCL_L))
else
TCL_L            = $(shell pkg-config --libs tcl)
$(info Others, looking for tcl)
$(info -> $(TCL_L))
TCL_INCLUDES     = 
endif

# EOF
