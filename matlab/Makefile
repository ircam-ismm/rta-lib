CC		= mex
CFLAGS		= 
RTA_ROOT	= ..
#MEXEXT		= mexmaci
#MEXEXT		= mexglx
#MEXEXT		= mexa64
#MEXEXT		= $(shell mexext)

ifdef DEBUG
CFLAGS += -g
else
CFLAGS += -O
endif

.PHONY: clean 

all: rta_bands_weights \
	rta_biquad rta_biquad_coefs \
	rta_dct_apply rta_dct_weights \
	rta_delta_apply rta_delta_weights \
	rta_downsample_int_mean \
	rta_fft \
	rta_fft_setup_delete rta_fft_setup_new \
	rta_ifft \
	rta_ifft_setup_delete rta_ifft_setup_new \
	rta_lifter_apply rta_lifter_weights \
	rta_lpc \
	rta_moments \
	rta_onepole \
	rta_preemphasis \
	rta_selection \
	rta_spectrum_to_bands \
	rta_svd \
	rta_var \
	rta_window_apply rta_window_weights \
	rta_yin \
	rta_yin_setup_delete rta_yin_setup_new


rta_%_mex.o: rta_%_mex.c
	$(CC) $(CFLAGS) -c -I. -I$(RTA_ROOT) $<

rta_%.o: $(RTA_ROOT)/rta_%.c $(RTA_ROOT)/rta_%.h 
	$(CC) $(CFLAGS) -c -I. -I$(RTA_ROOT) $<


rta_bands_weights: rta_bands.o rta_mel.o rta_bands_weights_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_biquad: rta_biquad.o rta_biquad_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_biquad_coefs: rta_biquad.o rta_biquad_coefs_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_dct_apply: rta_dct.o rta_dct_apply_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_dct_weights: rta_dct.o rta_dct_weights_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_delta_apply: rta_delta.o rta_delta_apply_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_delta_weights: rta_delta.o rta_delta_weights_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_downsample_int_mean: rta_resample.o rta_downsample_int_mean_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_fft: rta_fft.o rta_int.o rta_fft_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_fft_setup_delete: rta_fft.o rta_int.o rta_fft_setup_delete_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_fft_setup_new: rta_fft.o rta_int.o rta_fft_setup_new_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_ifft: rta_fft.o rta_int.o rta_ifft_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_ifft_setup_delete: rta_fft.o rta_int.o rta_ifft_setup_delete_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_ifft_setup_new: rta_fft.o rta_int.o rta_ifft_setup_new_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_lifter_apply: rta_lifter.o rta_lifter_apply_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_lifter_weights: rta_lifter.o rta_lifter_weights_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_lpc: rta_lpc.o rta_correlation.o rta_lpc_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_moments: rta_moments.o rta_moments_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_onepole: rta_onepole.o rta_onepole_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_preemphasis: rta_preemphasis.o rta_preemphasis_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_selection: rta_selection.o rta_selection_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_spectrum_to_bands: rta_bands.o rta_mel.o \
	rta_spectrum_to_bands_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_svd: rta_svd.o rta_svd_mex.o rta_int.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_var: rta_mean_variance.o rta_var_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_window_apply: rta_window.o rta_window_apply_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_window_weights: rta_window.o rta_window_weights_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_yin: rta_yin.o rta_correlation.o rta_yin_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_yin_setup_delete: rta_yin.o rta_correlation.o rta_yin_setup_delete_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@

rta_yin_setup_new: rta_yin.o rta_correlation.o rta_yin_setup_new_mex.o
	$(CC) $(CFLAGS) -cxx $^ -o $@


clean:
	rm -rf *.o