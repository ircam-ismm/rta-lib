CC		= mex
CFLAGS		=
RTA_ROOT	= ..
RTA_STD_HEADERS	= $(RTA_ROOT)/rta.h \
	$(RTA_ROOT)/rta_stdlib.h $(RTA_ROOT)/rta_types.h

ifdef DEBUG
CFLAGS += -g
endif

.PHONY: clean 

all: rta_yin rta_yin_setup_new rta_yin_setup_delete

rta_yin: rta_yin.o rta_correlation.o rta_mex_yin.o
	$(CC) $(CFLAGS) -cxx $^ -o rta_yin

rta_yin_setup_new: rta_yin.o rta_mex_yin_setup_new.o
	$(CC) $(CFLAGS) -cxx $^ -o rta_yin_setup_new

rta_yin_setup_delete: rta_yin.o rta_mex_yin_setup_delete.o
	$(CC) $(CFLAGS) -cxx $^ -o rta_yin_setup_delete

rta_yin.o: $(RTA_STD_HEADERS) $(RTA_ROOT)/rta_yin.h $(RTA_ROOT)/rta_yin.c 
	$(CC) $(CFLAGS) -c -I. -I$(RTA_ROOT) $(RTA_ROOT)/rta_yin.c

rta_mex_yin.o: rta_mex_yin.c
	$(CC) $(CFLAGS) -c -I. -I$(RTA_ROOT) $<

rta_mex_yin_setup_new.o: rta_mex_yin_setup_new.c
	$(CC) $(CFLAGS) -c -I. -I$(RTA_ROOT) $<

rta_mex_yin_setup_delete.o: rta_mex_yin_setup_delete.c
	$(CC) $(CFLAGS) -c -I. -I$(RTA_ROOT) $<

rta_correlation.o: $(RTA_STD_HEADERS) \
	$(RTA_ROOT)/rta_correlation.h $(RTA_ROOT)/rta_correlation.c
	$(CC) $(CFLAGS) -c -I. -I$(RTA_ROOT) $(RTA_ROOT)/rta_correlation.c

clean:
	rm -rf *.o *.mexmac *.mexa64
