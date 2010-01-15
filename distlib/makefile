
doc:
	doxygen distlib.doxyfile; $(MAKE) tex
tex:
	make -C doc/latex/

MIFVERSION = 0.5

mif:	mifclean mifquery
	/usr/bin/time ./mifquery -1 -1 -1 -1 DescFile.4001 DescFile.4001 100 5

mifprofile:
	/usr/bin/time ./mifquery -1 -1 -1 -1 DescFile.4001 DescFile.4001 100 5 result.txt \
	| & tee log.txt | awk '/\#bytes.*in index/ { print $7 } /\#bytes.*in optimal index/ { print $8 }'

mif-test:	libmif.o mifclean miftest
	./miftest

mif.o:	mif.c mif.h rta.h
	gcc -c -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -o $@ $<

../rta_util.o:	../rta_util.c ../rta_util.h ../rta.h
	gcc -c -g -Wall -I. -I.. -o $@ $<

libmif.o:	mif.o rta_util.o
	libtool -static -o $@ $^

miftest:	mif.c mif.h rta_util.c
	gcc -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -DMIF_BUILD_TEST -o $@ $^

mifquery:	mifquery.c mif.c mif.h rta_util.c
	gcc -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -lm -o $@ $^

mifclean:
	-rm miftest mifquery

mifdoc:
	doxygen mif.doxyfile; $(MAKE) miftex
miftex:
	make -C mifdoc/latex/

mifdist:
	tar cvfz mif-$(MIFVERSION)-src.tar.gz mif.c mif.h mifquery.c makefile rta_configuration.h -C .. rta.h rta_stdio.h rta_stdlib.h rta_util.{c,h} rta_types.h rta_math.h
#	zip -r mif-0.3-src.zip mif.c mif.h mifquery.c makefile ../rta.h ../rta_stdio.h ../rta_stdlib.h ../rta_util.[ch] rta_configuration.h

PUBLISH_ROOT      = DISCO/mifdoc
PUBLISH_FILES	  = mifdoc/html/\*
PUBLISH_HOST      = recherche.ircam.fr
PUBLISH_USER      = wwwas

mifpub:	mifdoc
	echo mput $(PUBLISH_FILES) | ncftp -u$(PUBLISH_USER) -p$(PW) ftp://$(PUBLISH_HOST)/$(PUBLISH_ROOT)


.PHONY: doc mifdoc
