function rta_compile()
%function rta_compile()
%     alternative to Makefile to create mex files on any platform
%     supported by Matlab
%

% RTA_ROOT = '..' ;

mex -O -I. -I.. ../rta_bands.c ../rta_mel.c rta_bands_weights_mex.c -o rta_bands_weights
mex -O -I. -I.. ../rta_biquad.c rta_biquad_mex.c -o rta_biquad
mex -O -I. -I.. ../rta_biquad.c rta_biquad_coefs_mex.c -o rta_biquad_coefs
mex -O -I. -I.. ../rta_dct.c rta_dct_apply_mex.c -o rta_dct_apply
mex -O -I. -I.. ../rta_dct.c rta_dct_weights_mex.c -o rta_dct_weights
mex -O -I. -I.. ../rta_delta.c rta_delta_apply_mex.c -o rta_delta_apply
mex -O -I. -I.. ../rta_delta.c rta_delta_weights_mex.c -o rta_delta_weights
mex -O -I. -I.. ../rta_resample.c rta_downsample_int_mean_mex.c -o rta_downsample_int_mean
mex -O -I. -I.. ../rta_fft.c rta_fft_mex.c -o rta_fft
mex -O -I. -I.. ../rta_fft.c rta_fft_setup_delete_mex.c -o rta_fft_setup_delete
mex -O -I. -I.. ../rta_fft.c ../rta_int.c rta_fft_setup_new_mex.c -o rta_fft_setup_new
mex -O -I. -I.. ../rta_fft.c rta_ifft_mex.c -o rta_ifft
mex -O -I. -I.. ../rta_fft.c rta_ifft_setup_delete_mex.c -o rta_ifft_setup_delete
mex -O -I. -I.. ../rta_fft.c ../rta_int.c rta_ifft_setup_new_mex.c -o rta_ifft_setup_new
mex -O -I. -I.. ../rta_lifter.c rta_lifter_apply_mex.c -o rta_lifter_apply
mex -O -I. -I.. ../rta_lifter.c rta_lifter_weights_mex.c -o rta_lifter_weights
mex -O -I. -I.. ../rta_lpc.c ../rta_correlation.c rta_lpc_mex.c -o rta_lpc
mex -O -I. -I.. ../rta_moments.c rta_moments_mex.c -o rta_moments
mex -O -I. -I.. ../rta_onepole.c rta_onepole_mex.c -o rta_onepole
mex -O -I. -I.. ../rta_preemphasis.c rta_preemphasis_mex.c -o rta_preemphasis
mex -O -I. -I.. ../rta_selection.c rta_selection_mex.c -o rta_selection
mex -O -I. -I.. ../rta_bands.c ../rta_mel.c rta_spectrum_to_bands_mex.c -o rta_spectrum_to_bands
mex -O -I. -I.. ../rta_svd.c rta_svd_mex.c ../rta_int.c -o rta_svd
mex -O -I. -I.. ../rta_mean_variance.c rta_var_mex.c -o rta_var
mex -O -I. -I.. ../rta_window.c rta_window_apply_mex.c -o rta_window_apply
mex -O -I. -I.. ../rta_window.c rta_window_weights_mex.c -o rta_window_weights
mex -O -I. -I.. ../rta_yin.c ../rta_correlation.c rta_yin_mex.c -o rta_yin
mex -O -I. -I.. ../rta_yin.c rta_yin_setup_delete_mex.c -o rta_yin_setup_delete
mex -O -I. -I.. ../rta_yin.c rta_yin_setup_new_mex.c -o rta_yin_setup_new
