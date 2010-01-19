
doc:
	doxygen distlib.doxyfile; $(MAKE) tex
tex:
	make -C doc/latex/

MIFVERSION = 0.6

mif:	mifclean discotrunc mifquery
	/usr/bin/time ./mifquery -1 -1 -1 -1 DescFile.4001 DescFile.4001 100 20

mif10:	mifquery
	/usr/bin/time ./mifquery -1 -1 -1 -1 DESCRIPTOR_SubSet_10000.4001 QueryData.4001 100 20 result10.log 2>&1| tee log10.txt
mif100:	mifquery
	/usr/bin/time ./mifquery -1 -1 -1 -1 DESCRIPTOR_SubSet_100000.4001 QueryData.4001 100 20 result100.log 2>&1| tee log100.txt

mifprofile:
	/usr/bin/time ./mifquery -1 -1 -1 -1 DescFile.4001 DescFile.4001 100 5 result.txt \
	| & tee log.txt | awk '/\#bytes.*in index/ { print $7 } /\#bytes.*in optimal index/ { print $8 }'

mif-test:	libmif.o mifclean miftest
	./miftest

#### compilation

mif.o:	mif.c mif.h rta.h
	gcc -c -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -o $@ $<

discofile.o:	discofile.c discofile.h
	gcc -c -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -o $@ $<

../rta_util.o:	../rta_util.c ../rta_util.h ../rta.h
	gcc -c -g -Wall -I. -I.. -o $@ $<

libmif.o:	mif.o discofile.o rta_util.o 
	libtool -static -o $@ $^


#### executables

miftest:	mif.c mif.h rta_util.c
	gcc -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -DMIF_BUILD_TEST -o $@ $^

mifquery:	mifquery.c mif.c mif.h discodist.c discodist.h discofile.c discofile.h rta_util.c
	gcc -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -lm -o $@ $^

discotrunc:	discotrunc.c discofile.c discofile.h
	gcc -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -o $@ $^

#### profiling

mifquery.prof:	mifquery.c mif.c mif.h discofile.c discofile.h rta_util.c
#mac:	gcc -finstrument-functions -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -lm -o $@ $^ /usr/lib/libSaturn.dylib
	gcc -pg -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -lm -o $@ $^
	/usr/bin/time ./$@ -1 -1 -1 -1 DESCRIPTOR_SubSet_10000.4001 QueryData.4001 10 20 result10.log 2>&1| tee log10k.txt
	gprof mifquery.prof gmon.out > gprof-$(MIFVERSION)-10k.txt

#### misc

mifclean:
	-rm miftest mifquery mifquery.prof

mifdoc:
	doxygen mif.doxyfile; $(MAKE) miftex
miftex:
	make -C mifdoc/latex/

mifdist:
	svn copy -m 'tagging metric inverted file index distribution $(MIFVERSION)' . http://trac.ircam.fr/repos/imtr-libs/tags/mif-$(MIFVERSION)
	tar cvfz mif-$(MIFVERSION)-src.tar.gz mif.c mif.h mifquery.c makefile rta_configuration.h -C .. rta.h rta_stdio.h rta_stdlib.h rta_util.{c,h} rta_types.h rta_math.h
#	zip -r mif-0.3-src.zip mif.c mif.h mifquery.c makefile ../rta.h ../rta_stdio.h ../rta_stdlib.h ../rta_util.[ch] rta_configuration.h

PUBLISH_ROOT      = DISCO/mifdoc
PUBLISH_FILES	  = mifdoc/html/\*
PUBLISH_HOST      = recherche.ircam.fr
PUBLISH_USER      = wwwas

mifpub:	mifdoc
	echo mput $(PUBLISH_FILES) | ncftp -u$(PUBLISH_USER) -p$(PW) ftp://$(PUBLISH_HOST)/$(PUBLISH_ROOT)


.PHONY: doc mifdoc
