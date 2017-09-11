/**
 * @file   rta_fft_setup_mex.h
 * @author Jean-Philippe.Lambert@m1267.ircam.fr
 * @date   Thu Jul 10 16:42:20 2008
 * 
 * @brief  
 * 
 * 
 */

struct rta_fft_setup_mex
{
  rta_fft_setup_t * fft_setup;
  rta_complex_t * fft;
  unsigned int fft_size;
  rta_real_t nyquist;
  rta_real_t scale;
};

typedef struct rta_fft_setup_mex rta_fft_setup_mex_t;

