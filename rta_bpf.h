
#ifndef _RTA_BPF_H_
#define _RTA_BPF_H_

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif

/** 
 * Break-point function (or time-tagged values).
 *
 * For the moment read only, since construction is done by binary-compatible bpf_t from FTM
 */

/* a single bpf point */
typedef struct rta_bpf_point
{
  double time; /**< absolute break point time */
  double value; /**< break point value */
  double slope; /**< slope to next value */
} rta_bpf_point_t;

/** the break-point function itself 
 *
 * ATTENTION: must be binary-compatible with FTM struct bpfunc in ftmlib/classes/bpf.c!!!
 */
typedef struct _rta_bpf
{
  rta_bpf_point_t *points;	/**< break points ... */
  int alloc;			/**< alloc ... */
  int size;			/**< size ... */
  int index;			/**< index cache for get_interpolated method */
} rta_bpf_t;


#define rta_bpf_get_size(b) ((b)->size)
#define rta_bpf_get_time(b, i) ((b)->points[i].time)
#define rta_bpf_get_value(b, i) ((b)->points[i].value)
#define rta_bpf_get_slope(b, i) ((b)->points[i].slope)
#define rta_bpf_get_duration(b) ((b)->points[(b)->size - 1].time)

double rta_bpf_get_interpolated (rta_bpf_t *bpf, double time);

#ifdef __cplusplus
}
#endif

#endif /* _RTA_BPF_H_ */
