include $(AVERTEC_TOP)/etc/$(AVERTEC_OS).mk

all : lib
	@cd xyagle ; $(MAKE) -f Makefile

bin : lib
	@cd fcl ; $(MAKE) -f Makefile
	@cd genius ; $(MAKE) -f Makefile
	@cd yagle ; $(MAKE) -f Makefile
	@cd slib ; $(MAKE) -f Makefile
	@cd cgv ; $(MAKE) -f Makefile
	@cd xyagle ; $(MAKE) -f Makefile
	@$(MAKE) -C gnshr

lib :
	@cd fcl ; $(MAKE) -f Makefile
	@cd genius ; $(MAKE) -f Makefile
	@cd yagle ; $(MAKE) -f Makefile lib
	@cd slib ; $(MAKE) -f Makefile
	@cd cgv ; $(MAKE) -f Makefile
	@cd xyagle ; $(MAKE) -f Makefile lib

headers :
	@cd fcl ; $(MAKE) -f Makefile headers
	@cd genius ; $(MAKE) -f Makefile headers
	@cd yagle ; $(MAKE) -f Makefile headers
	@cd slib ; $(MAKE) -f Makefile headers
	@cd cgv ; $(MAKE) -f Makefile headers
	@cd xyagle ; $(MAKE) -f Makefile headers

clean :
	@cd fcl ; $(MAKE) -f Makefile clean
	@cd genius ; $(MAKE) -f Makefile clean
	@$(MAKE) -C gnshr clean
	@cd yagle ; $(MAKE) -f Makefile clean
	@cd slib ; $(MAKE) -f Makefile clean
	@cd cgv ; $(MAKE) -f Makefile clean
	@cd xyagle ; $(MAKE) -f Makefile clean

realclean :
	@cd fcl ; $(MAKE) -f Makefile realclean
	@cd genius ; $(MAKE) -f Makefile realclean
	@$(MAKE) -C gnshr realclean
	@cd yagle ; $(MAKE) -f Makefile realclean
	@cd slib ; $(MAKE) -f Makefile realclean
	@cd cgv ; $(MAKE) -f Makefile realclean
	@cd xyagle ; $(MAKE) -f Makefile realclean

dep :
	@cd fcl ; $(MAKE) -f Makefile dep
	@cd genius ; $(MAKE) -f Makefile dep
	@$(MAKE) -C gnshr dep
	@cd yagle ; $(MAKE) -f Makefile dep
	@cd slib ; $(MAKE) -f Makefile dep
	@cd cgv ; $(MAKE) -f Makefile dep
	@cd xyagle ; $(MAKE) -f Makefile dep

build_gnstest:
	@$(MAKE) -C genius build_gnstest
