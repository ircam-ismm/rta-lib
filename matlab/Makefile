CC		= mex
CFLAGS		=
RTA_ROOT	= ..
MEXEXT		= $(shell mexext)
ifdef DEBUG
CFLAGS += -g
endif

.PHONY: clean 

all: rta_downsample_int_mean.$(MEXEXT) \
	rta_preemphasis.$(MEXEXT) \
	rta_fft.$(MEXEXT) \
	rta_fft_setup_new.$(MEXEXT) rta_fft_setup_delete.$(MEXEXT) \
	rta_window_weights.$(MEXEXT) rta_window_apply.$(MEXEXT) \
	rta_yin.$(MEXEXT) \
	rta_yin_setup_new.$(MEXEXT) rta_yin_setup_delete.$(MEXEXT)

rta_%_mex.o: rta_%_mex.c
	$(CC) $(CFLAGS) -c -I. -I$(RTA_ROOT) $<

rta_%.o: $(RTA_ROOT)/rta_%.c $(RTA_ROOT)/rta_%.h 
	$(CC) $(CFLAGS) -c -I. -I$(RTA_ROOT) $<


rta_downsample_int_mean.$(MEXEXT): rta_resample.o rta_downsample_int_mean_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_fft.$(MEXEXT): rta_fft.o rta_fft_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_fft_setup_delete.$(MEXEXT): rta_fft.o rta_fft_setup_delete_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_fft_setup_new.$(MEXEXT): rta_fft.o rta_complex.o rta_int.o rta_fft_setup_new_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_preemphasis.$(MEXEXT): rta_preemphasis.o rta_preemphasis_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_window_apply.$(MEXEXT): rta_window.o rta_window_apply_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_window_weights.$(MEXEXT): rta_window.o rta_window_weights_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_yin.$(MEXEXT): rta_yin.o rta_correlation.o rta_yin_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_yin_setup_delete.$(MEXEXT): rta_yin.o rta_yin_setup_delete_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_yin_setup_new.$(MEXEXT): rta_yin.o rta_yin_setup_new_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@


clean:
	rm -rf *.o *.$(MEXEXT)
