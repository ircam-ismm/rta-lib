
doc:
	doxygen distlib.doxyfile; $(MAKE) tex
tex:
	make -C doc/latex/

mif:	mifclean miftest
	./miftest

miftest:	mif.c mif.h ../rta_util.c
	gcc -g -Wall -I. -I.. -DMIF_BUILD_TEST -o $@ $^

mifclean:
	rm miftest

mifdoc:
	doxygen mif.doxyfile; $(MAKE) miftex
miftex:
	make -C mifdoc/latex/

.PHONY: doc mifdoc
