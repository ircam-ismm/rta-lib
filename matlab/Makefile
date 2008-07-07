CC		= mex
CFLAGS		=
RTA_ROOT	= ..

ifdef DEBUG
CFLAGS += -g
endif

.PHONY: clean 

all: rta_downsample_int_mean rta_preemphasis \
	rta_yin rta_yin_setup_new rta_yin_setup_delete


rta_%_mex.o: rta_%_mex.c
	$(CC) $(CFLAGS) -c -I. -I$(RTA_ROOT) $<

rta_%.o: $(RTA_ROOT)/rta_%.h $(RTA_ROOT)/rta_%.c
	$(CC) $(CFLAGS) -c -I. -I$(RTA_ROOT) $^

rta_%: rta_%.o rta_%_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $<


rta_downsample_int_mean: rta_resample.o rta_downsample_int_mean_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o rta_downsample_int_mean

rta_yin: rta_yin.o rta_correlation.o rta_yin_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o rta_yin

rta_yin_setup_new: rta_yin.o rta_yin_setup_new_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o rta_yin_setup_new

rta_yin_setup_delete: rta_yin.o rta_yin_setup_delete_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o rta_yin_setup_delete

clean:
	rm -rf *.o *.mexmac *.mexa64
