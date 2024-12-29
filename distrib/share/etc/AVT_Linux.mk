RPC_L = -lcrypt
RPCGEN = /usr/bin/rpcgen
RPCGENFLAGS = -C

LDL = -ldl
BIN_EXT=

OSCFLAGS = -fno-inline

DYNAMIC = -rdynamic
STATIC  = -rstatic

GENDOCS     = /users10/chaos1/avertec/gendocs
LIB_TERMCAP = -ltermcap

ifeq ($(BUILD_VARIANT)$(LIB_SUFFIX_),Linux.slsoc6x_64)
  JAVA_HOME = /usr/lib/jvm/java-1.6.0-openjdk.x86_64
  JAVA      = $(JAVA_HOME)/bin/java
  SAXON     = $(JAVA) -jar /usr/share/java/saxon9.jar
else
  ifeq ($(BUILD_VARIANT)$(LIB_SUFFIX_),Linux.slsoc6x)
    JAVA_HOME = /usr/lib/jvm/java-1.6.0-openjdk
    JAVA      = $(JAVA_HOME)/bin/java
    SAXON     = $(JAVA) -jar /usr/share/java/saxon9.jar
  else
    ifeq ($(BUILD_VARIANT),Linux.el7)
      JAVA_HOME = /usr/lib/jvm/java-1.7.0-openjdk
      JAVA      = $(JAVA_HOME)/bin/java
      SAXON     = $(JAVA) -jar /usr/share/java/saxon.jar
    else
      ifeq ($(BUILD_VARIANT),Linux.el8)
        JAVA_HOME = /usr/lib/jvm/java-1.8.0-openjdk
        JAVA      = $(JAVA_HOME)/bin/java
        SAXON     = $(JAVA) -jar /usr/share/java/saxon.jar
      else
        ifeq ($(BUILD_VARIANT),Linux.el9)
          JAVA_HOME = /usr
          JAVA      = $(JAVA_HOME)/bin/java
          SAXON     = $(JAVA) -jar ${AVERTEC_TOP}/../distrib_extras/saxon9.jar
        else
          ifeq ($(BUILD_VARIANT),Linux.ubuntu)
            JAVA_HOME = /usr/lib/jvm/java-1.6.0-openjdk
            JAVA      = $(JAVA_HOME)/bin/java
            SAXON     = saxonb-xslt -ext:on
           #SAXON     = CLASSPATH=/usr/share/java/saxonb.jar $(JAVA) net.sf.saxon.Transform
          else
            ifeq ($(BUILD_VARIANT),Linux.openSUSE)
              JAVA_HOME   = /usr
              JAVA        = $(JAVA_HOME)/bin/java
              SAXON       = $(JAVA) -jar ${AVERTEC_TOP}/../distrib_extras/saxon9.jar
              LIB_TERMCAP = /usr/lib64/libtermcap.so.2
            else
              JAVA_HOME = /usr/lib/jvm/java-1.6.0-openjdk
              JAVA      = $(JAVA_HOME)/bin/java
              SAXON     = $(JAVA) -jar /usr/share/java/saxon9.jar
            endif
          endif
        endif
      endif
    endif
  endif
endif

FOP = export JAVA_HOME=$(JAVA_HOME);/usr/bin/fop


ifndef FLEX_HEARTBEAT
FLEXOSLIBS = -lpthread
else
FLEXOSLIBS = 
endif

LICENSE_API =license_api

