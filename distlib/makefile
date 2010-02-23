
doc:
	doxygen distlib.doxyfile; $(MAKE) tex
tex:
	make -C doc/latex/

MIFVERSION = 0.8
z	   = z
CFLAGS	   = -g -O3 -Wall -DMIFVERSION=$(MIFVERSION) -DUSE_ZLIB=1 -I. -I..
data       = data

mif:	mifclean discotrunc mifindex mifquery mifindex-test mifquery-test

mifindex-test:	mifindex
	/usr/bin/time ./mifindex -1 -1 testindex.db$z DescFile.4001 DescFile.4001 DescFile.4001

miftest-mini: mifindex mifquery
	/usr/bin/time ./mifindex -1 -1 testmini.db$z DescFile-20.4001
	/usr/bin/time ./mifquery -1 -1 testmini.db$z DescFile-20.4001 10 5

mifindex-all:	mifindex
	/usr/bin/time ./mifindex -1 -1  ${data}/testindex10k.db$z   DESCRIPTOR_SubSet_10000.4001   2>&1| tee testindex10k$z.log   &
	/usr/bin/time ./mifindex -1 -1  ${data}/testindex100k.db$z  DESCRIPTOR_SubSet_100000.4001  2>&1| tee testindex100k$z.log  &
	/usr/bin/time ./mifindex -1 100 ${data}/testindex500k.db$z  DESCRIPTOR_SubSet_500000.4001  2>&1| tee testindex500k$z.log  &
	/usr/bin/time ./mifindex -1 150 ${data}/testindex1000k.db$z DESCRIPTOR_SubSet_1000000.4001 2>&1| tee testindex1000k$z.log &

mifindex-10:	mifindex
	/usr/bin/time ./mifindex -1 -1  ${data}/testindex10k.db$z   DESCRIPTOR_SubSet_10000.4001   2>&1| tee testindex10k$z.log
mifindex-100:	mifindex
	/usr/bin/time ./mifindex -1 -1  ${data}/testindex100k.db$z  DESCRIPTOR_SubSet_100000.4001  2>&1| tee testindex100k$z.log
mifindex-500:	mifindex
	/usr/bin/time ./mifindex -1 100 ${data}/testindex500k.db$z  DESCRIPTOR_SubSet_500000.4001  2>&1| tee testindex500k$z.log
mifindex-1000:	mifindex
	/usr/bin/time ./mifindex -1 150 ${data}/testindex1000k.db$z DESCRIPTOR_SubSet_1000000.4001 2>&1| tee testindex1000k$z.log

mifquery-test:	mifquery
	/usr/bin/time ./mifquery -1 -1 testindex.db$z DescFile.4001 10 5

mifquery-all:	mifquery
	/usr/bin/time ./mifquery -1 -1 ${data}/testindex10k.db$z   QueryData.4001 100 20 result10$z.txt   2>&1| tee testquery10k$z.log   &
	/usr/bin/time ./mifquery -1 -1 ${data}/testindex100k.db$z  QueryData.4001 100 20 result100$z.txt  2>&1| tee testquery100k$z.log  &
	/usr/bin/time ./mifquery 100 5 ${data}/testindex500k.db$z  QueryData.4001 100 20 result500$z.txt  2>&1| tee testquery500k$z.log  &
	/usr/bin/time ./mifquery 150 5 ${data}/testindex1000k.db$z QueryData.4001 100 20 result1000$z.txt 2>&1| tee testquery1000k$z.log &


mif10:	mifquery								       
	/usr/bin/time ./mifquery -1 -1  ${data}/testindex10k.db$z   QueryData.4001 100 20 result10$z.txt   2>&1| tee testquery10$z.log
mif100:	mifquery}			
	/usr/bin/time ./mifquery -1 -1  ${data}/testindex100k.db$z  QueryData.4001 100 20 result100$z.txt  2>&1| tee testquery100$z.log
mif500:	mifquery			
	/usr/bin/time ./mifquery  100 5 ${data}/testindex500k.db$z  QueryData.4001 100 20 result500$z.txt  2>&1| tee testquery500$z.log
mif1000: mifquery			
	/usr/bin/time ./mifquery  150 5 ${data}/testindex1000k.db$z QueryData.4001 100 20 result1000$z.txt 2>&1| tee testquery1000$z.log

mifprof mifprofile:
	/usr/bin/time ./mifquery -1 -1 DescFile.4001 DescFile.4001 100 5 result.txt \
	| & tee log.txt | awk '/\#bytes.*in index/ { print $7 } /\#bytes.*in optimal index/ { print $8 }'

mif-test:	libmif.o mifclean miftest
	./miftest

#### compilation

mif.o:	mif.c mif.h rta.h
	gcc -c -O3 -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -o $@ $<

discofile.o:	discofile.c discofile.h
	gcc -c -O3 -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -o $@ $<

../rta_util.o:	../rta_util.c ../rta_util.h ../rta.h
	gcc -c -O3 -g -Wall -I. -I.. -o $@ $<

libmif.o:	mif.o discofile.o rta_util.o 
	libtool -static -o $@ $^


#### executables

miftest:	mif.c mif.h rta_util.c
	gcc $(CFLAGS) -DMIF_BUILD_TEST -o $@ $^

mifindex:	mifindex.c mif.c mif.h mifhash.h mifdb.h mifdbsqlite3.c discodist.c discodist.h discofile.c discofile.h rta_util.c
	gcc $(CFLAGS) -lm -lsqlite3 -lz -o $@ $^

mifquery:	mifquery.c mif.c mif.h mifhash.h mifdb.h mifdbsqlite3.c discodist.c discodist.h discofile.c discofile.h rta_util.c

	gcc $(CFLAGS) -lm -lsqlite3 -lz -o $@ $^

discotrunc:	discotrunc.c discofile.c discofile.h
	gcc $(CFLAGS) -o $@ $^

#### profiling

mifquery.prof:	mifquery.c mif.c mif.h discodist.c discodist.h discofile.c discofile.h rta_util.c
#mac:	gcc -finstrument-functions -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -lm -o $@ $^ /usr/lib/libSaturn.dylib
	gcc -pg -O3 -g -Wall -DMIFVERSION=$(MIFVERSION) -I. -I.. -lm -o $@ $^
	/usr/bin/time ./$@ -1 -1 -1 -1 DESCRIPTOR_SubSet_10000.4001 QueryData.4001 10 20 result-$(MIFVERSION)-10k$z.log 2>&1| tee log-$(MIFVERSION)-10k.txt
	gprof mifquery.prof gmon.out > gprof-$(MIFVERSION)-10k.txt

#### misc

tags etags:
	etags *.[ch]
mifclean:
	-rm miftest mifquery mifquery.prof

mifdoc:
	doxygen mif.doxyfile; $(MAKE) miftex
miftex:
	TEXINPUTS='' make -C mifdoc/latex/

MIFREPOSITORY = http://trac.ircam.fr/repos/imtr-libs/tags/mif-$(MIFVERSION)
miftag:
	svn remove -m 'overwrite tag $(MIFVERSION)' $(MIFREPOSITORY)
	svn copy -m 'tagging metric inverted file index distribution $(MIFVERSION)' . $(MIFREPOSITORY)

mifdist: miftag
	tar cvfz mif-$(MIFVERSION)-src.tar.gz mif.c mif.h mifhash.h mifquery.c makefile rta_configuration.h -C .. rta.h rta_stdio.h rta_stdlib.h rta_util.{c,h} rta_types.h rta_math.h

PUBLISH_ROOT      = DISCO/mifdoc
PUBLISH_FILES	  = mifdoc/html/\*
PUBLISH_HOST      = recherche.ircam.fr
PUBLISH_USER      = wwwas
PW = $(shell perl -ne 'm/^username=(.*)/ && do {$$u=$$1}; m/^password=(.*)/ && ($$u eq "$(PUBLISH_USER)") && do { print $$1; exit; }' $(HOME)/.gftp/bookmarks)

mifpub:	mifdoc
	echo mput $(PUBLISH_FILES) | ncftp -u$(PUBLISH_USER) -p$(PW) ftp://$(PUBLISH_HOST)/$(PUBLISH_ROOT)


.PHONY: doc mif mifdoc mifindex
