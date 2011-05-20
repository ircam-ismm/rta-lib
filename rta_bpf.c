
#include "rta_bpf.h"

static int rta_bf_get_index (rta_bpf_t *bpf, double time)
{
  int size  = rta_bpf_get_size(bpf);
  int index = bpf->index;
  
  if (index >= size - 1)
    index = size - 2;
  
  /* search index */
  if (time >= rta_bpf_get_time(bpf, index + 1))
  {
    index++;
	
    while (time >= rta_bpf_get_time(bpf, index + 1))
      index++;
  }
  else if (time < rta_bpf_get_time(bpf, index))
  {
    index--;
	
    while (time < rta_bpf_get_time(bpf, index))
      index--;
  }
  else if (rta_bpf_get_slope(bpf, index) == DBL_MAX)
  {
    index++;
	
    while (rta_bpf_get_slope(bpf, index) == DBL_MAX)
      index++;
  }
  
  bpf->index = index;
  
  return index;
}


double rta_bpf_get_interpolated (rta_bpf_t *self, double time)
{
  double duration = rta_bpf_get_duration(self);
  
  if (time <= rta_bpf_get_time(self, 0))
    return rta_bpf_get_value(self, 0);
  else if (time >= duration)
  {
    int size = rta_bpf_get_size(self);
    return rta_bpf_get_value(self, size - 1);
  }
  else
  {
    int index = rta_bf_get_index(self, time); 
    return rta_bpf_get_value(self, index) + (time - rta_bpf_get_time(self, index)) * rta_bpf_get_slope(self, index);
  }
}
