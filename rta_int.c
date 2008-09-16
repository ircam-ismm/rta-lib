/**
 * @file   rta_int.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu Sep 12 18:10:41 2007
 * 
 * @brief  Integer mathematical functions
 * 
 *
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */
#include "rta_int.h"


int 
rta_ilog2(unsigned int n)
{
  unsigned int log2 = 0;
  int i;
  
  for(i=n>>1; i; i>>=1)
  {
    log2++;
  }

  return log2;
}

inline int rta_imax(int m, int n)
{
  return (m > n) ? m : n;
}

inline int rta_imin(int m, int n)
{
  return (m < n) ? m : n;
}

unsigned int
rta_inextpow2(unsigned int n)
{
  unsigned int pow2 = 2;
  int i;

  for(i=((n-1)>>1); i>0; i>>=1)
  {
    pow2 <<= 1;
  }

  return pow2;
}


