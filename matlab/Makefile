CC		= mex
RTA_ROOT	= ..
RTA_STD_HEADERS	= $(RTA_ROOT)/rta.h \
	$(RTA_ROOT)/rta_stdlib.h $(RTA_ROOT)/rta_types.h

.PHONY: clean

all: rta_yin

rta_yin: rta_yin.o rta_correlation.o rta_mex_yin.o
	$(CC) -cxx $^ -o rta_yin

rta_yin.o: $(RTA_STD_HEADERS) $(RTA_ROOT)/rta_yin.h $(RTA_ROOT)/rta_yin.c 
	$(CC) -c -I. -I$(RTA_ROOT) $(RTA_ROOT)/rta_yin.c

rta_mex_yin.o: rta_mex_yin.c
	$(CC) -c -I. -I$(RTA_ROOT) $<

rta_correlation.o: $(RTA_STD_HEADERS) \
	$(RTA_ROOT)/rta_correlation.h $(RTA_ROOT)/rta_correlation.c
	$(CC) -c -I. -I$(RTA_ROOT) $(RTA_ROOT)/rta_correlation.c

clean:
	rm -rf *.o *.mexmac *.mexa64
