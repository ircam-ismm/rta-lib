/*
 * FTM 
 * Copyright (C) 1994, 1995, 1998, 1999, 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * 
 * See file COPYING.LIB for further informations on licensing terms.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * 
 */

#ifdef WIN32
#include <malloc.h>
#else
#include <alloca.h>
#endif
#include <math.h>

#include "fts.h"


#define MAX_FLOAT 0x7FFFFFFF


const char *kdtree_dmodestr[] = { "orthogonal", "hyperplane", "pca" };
const char *kdtree_mmodestr[] = { "mean", "middle", "median" };

fts_symbol_t kdtree_symbol = 0;
fts_class_t *kdtree_class  = 0;

fts_symbol_t s_mean;
fts_symbol_t s_hyper;
fts_symbol_t s_vec;
fts_symbol_t s_profile;
fts_symbol_t s_getknni;


#if PROFILE
void profile_clear (kdtree_t *t)
{
    t->profile.v2v   	  = 0; 
    t->profile.v2n   	  = 0;
    t->profile.mean  	  = 0;
    t->profile.hyper 	  = 0;
    t->profile.searches   = 0;
    t->profile.neighbours = 0;
    t->profile.maxstack   = 0;
}
#endif

static double 
vec_dist (fmat_t *a, fmat_t *b)
{
    int    m 	= fmat_get_m(a);
    int    n 	= fmat_get_n(a);
    float *l 	= fmat_get_ptr(a);
    float *r 	= fmat_get_ptr(b);
    int    size = m * n;
    double dist = 0;
    int    i;
    
    for (i = 0; i < size; i++)
		dist += (l[i] - r[i]) * (l[i] - r[i]);

    return dist;    
}


static float 
distfmatV2H(fmat_t* vect, fmat_t* hplane) 
{
    //standard algebra computing
    int	  i;
    float tmp1 = 0;
    float tmp2 = 0;
    int   m  = fmat_get_m(vect);
    int   n  = fmat_get_n(vect);
    float *v = fmat_get_ptr(vect);
    float *h = fmat_get_ptr(hplane);
    int   size = m * n;
	
    for(i = 0; i < size; i++) 
    {
	tmp1 += v[i]*h[i];
	tmp2 += h[i]*h[i];
    }
	
    return (tmp1 + h[size])/sqrt(tmp2);
}


void vec_post (float *v, int stride, int n, const char *suffix)
{
    int i, ii;

    for (i = 0, ii = 0; i < n; i++, ii += stride) 
    { 	
	fts_post("%s%.1f", (ii == 0  ?  "["  :  ", "), v[ii]);
    }
    fts_post("]%s", suffix);
}

void kdtree_info_display (kdtree_t* t) 
{
#   define MB(b)  ((float) (b) / (float) (1024 * 1024))
#   define FLT(b) ((b) * sizeof(float))

    float mbdata  = MB(FLT(t->dataalloc));
    float mbindex = MB(t->dataalloc   * sizeof(int));
    float mbstack = MB(t->stack.alloc * sizeof(kdtree_stack_elem_t));
    float mbnodes = MB(t->nnodes      * sizeof(node_t));
    /* inner nodes' mean vectors and splitplanes 
       (these only in hyperplane mode) */
    float mbinner = MB(t->ninner * FLT(t->ndim) * 
		       (t->dmode == dmode_orthogonal ? 1 : 2));

    fts_post("\nTree Info:\n");
    fts_post("ndim        = %d\n", t->ndim);
    fts_post("ndata       = %d  (%.3f MB extern alloc size)\n", t->ndata, mbdata);
    fts_post("nalloc      = %d  (%.3f MB index)\n",  t->dataalloc, mbindex);
    fts_post("maxheight   = %d\n", t->maxheight);
    fts_post("givenheight = %d\n", t->givenheight);
    fts_post("height      = %d\n", t->height);
    fts_post("nnodes      = %d  (%.3f MB node struct)\n", t->nnodes, mbnodes);
    fts_post("inner nodes = %d  (%.3f MB node vectors)\n", t->ninner, mbinner);
    fts_post("stack       = %d  (%.3f MB)\n", t->stack.alloc, mbstack);
    fts_post("total size  = %.3f MB\n", 
	     MB(sizeof(kdtree_t)) + mbnodes + mbinner + mbindex + mbstack);
    fts_post("sort mode     = %d\n", t->sort);
    fts_post("decomposition = %s\n", kdtree_dmodestr[t->dmode]);
    fts_post("mean vector   = %s\n", kdtree_mmodestr[t->mmode]);
}

void kdtree_raw_display (kdtree_t* t) 
{
    int i;

    if (t->height == 0 || t->ndata == 0) fts_post("Empty Tree\n");

    for (i = 0; i < t->ndata; i++) 
    {
	fts_post("raw data vec %-3i = ", i);
	vec_post(kdtree_get_row_ptr(t, i), 1, t->ndim, "\n");
    }
}

void kdtree_data_display(kdtree_t* t, int print_data) 
{
    float plane[t->ndim];
    int l, n, i;

    fts_post("\nTree Data:\n");
    if (t->height == 0 || t->ndata == 0) fts_post("Empty Tree\n");
	
    for (l = 0; l < t->height; l++) 
    {
	fts_post("Level #%d  nodes %d..%d  splitdim %d\n", 
		 l, pow2(l) - 1, pow2(l+1) - 2, t->nodes[pow2(l) - 1].splitdim);
	for (n = pow2(l) - 1; n < pow2(l+1) - 1; n++) 
	{
	  node_t *node = &t->nodes[n];

	  if (n < t->ninner)
	  {
	    fts_post("  inner node %d size %d <%d..%d> splitdim %d splitplane ",
		     n, node->size, node->startind, node->endind, node->splitdim);
		
	    if (t->dmode == dmode_orthogonal)
	    {	/* splitplane implicit orthogonal to splitdim */
		bzero(plane, t->ndim * sizeof(float));
		plane[node->splitdim] = 1;
		vec_post(plane, 1, t->ndim, "");
	    }
	    else
		vec_post(fmat_get_ptr(node->split), 1, t->ndim, "");
	  }
	  else
	    fts_post("  leaf node %d size %d <%d..%d> ",
		     n, node->size, node->startind, node->endind);

	  if (print_data)
	  {
	    fts_post(" = (");
	    for (i = node->startind; i <= node->endind; i++) 
	    {
		fts_post("%svec %d = ", (print_data >= 2  ?  "\n    " : ""),
			 t->dataindex[i]);
		vec_post(kdtree_get_vector(t, i), 1, t->ndim,
			 i < node->endind ? ", " : "");
	    }
	    fts_post(")");
	  }
	  fts_post("\n");
	} /* end for nodes n */
    } /* end for level l */
}


static void
kdtree_getelem_function(fts_object_t *o, int ac, const fts_atom_t *at, fts_atom_t *ret)
{  
  if(ac > 0 && fts_is_number(at))
  {
  }
}


/* free nodes */
void kdtree_clear_nodes (kdtree_t *self)
{
  int i;

  for (i = 0; i < self->nnodes; i++)
  {
      if (self->nodes[i].mean)	  fts_object_release(self->nodes[i].mean);
      if (self->nodes[i].split)   fts_object_release(self->nodes[i].split);
#if DEBUG
      self->nodes[i].mean  = NULL;
      self->nodes[i].split = NULL;
#endif
  }
}


void kdtree_set (kdtree_t *self, fmat_t *data, fmat_t *sigma, int use_sigma)
{
  int maxheight   = floor(log2(fmat_get_m(data)));
  int givenheight = self->givenheight;
  int height      = givenheight > 0  ?  givenheight :  maxheight + givenheight;

  /* clip height */
  if (height > maxheight)
      height = maxheight;
  if (height < 1)	
      height = 1;	/* minimum: just one node, linear search */
  self->maxheight = maxheight;

  kdtree_clear_nodes(self);

  if (self->data)
      fts_object_release((fts_object_t *) self->data);

  self->data  = data;
  fts_object_refer((fts_object_t *) self->data);

  kdtree_set_sigma(self, sigma);
  kdtree_init_data(self, height, fmat_get_m(data), fmat_get_n(data));
  kdtree_build(self, use_sigma); 
}

/* update non-zero sigma index list */
int kdtree_update_sigmanz (kdtree_t *self)
{
    float *sigmaptr = fmat_get_ptr(self->sigma);
    int j, nnz = 0;

    for (j = 0; j < self->ndim; j++)
	if (sigmaptr[j] != 0)
	    self->sigma_indnz[nnz++] = j;
    self->sigma_nnz = nnz;
    return nnz;
}

void kdtree_set_sigma (kdtree_t *self, fmat_t *sigma)
{
    if (self->sigma)
      fts_object_release((fts_object_t *) self->sigma);
 
    if (sigma)
	self->sigma = sigma;
    else
    {
	self->sigma = fmat_create(1, fmat_get_n(self->data));
	fmat_set_const(self->sigma, 1);
    }

    fts_object_refer((fts_object_t *) self->sigma);

    self->sigma_indnz = fts_realloc(self->sigma_indnz, 
                            fmat_get_n(self->sigma) * sizeof(*self->sigma));
    kdtree_update_sigmanz(self);
}


/******************************************************************************
 *
 * user methods
 *
 */

static fts_method_status_t _kdtree_set (fts_object_t *o, fts_symbol_t s, 
				      int ac, const fts_atom_t *at, 
				      fts_atom_t *ret)
{
  kdtree_t *self = (kdtree_t *) o;
  fmat_t   *data = (fmat_t *)   fts_get_object(at);

  kdtree_set(self, data, self->sigma, 1); 

  fts_set_object(ret, o);
  return fts_ok;
}


static fts_method_status_t _kdtree_add (fts_object_t *o, fts_symbol_t s, int ac, const fts_atom_t *at, fts_atom_t *ret)
{
  kdtree_t *self = (kdtree_t *) o;

  /* rebuild */
  /* kdtree_set(self, ac, at); 
     kdtree_build(self);
  */

  fts_set_object(ret, self);
  fts_object_refer(self);

  return fts_ok;
}



static fts_method_status_t _kdtree_getknn (fts_object_t *o, fts_symbol_t s, 
					 int ac, const fts_atom_t *at, 
					 fts_atom_t *ret)
{
  kdtree_t *self = (kdtree_t *) o;
  fmat_t *x = NULL;	/* search vector */
  int     k = 1, n;
  float  *result;
  float  *d;

  if (ac > 0  &&  fts_is_a(at, fmat_class))
      x = (fmat_t *) fts_get_object(at);

  if (ac > 1  &&  fts_is_number(at+1))
      k = fts_get_number_int(at+1);

  result = alloca(k * sizeof(float));
  d      = alloca(k * sizeof(float));

  n = kdtree_search_knn(self, fmat_get_ptr(x), 1, k, 0, result, d);

  if (s == fts_s_get) /* geti: nearest index? */
  {   /* return one fvec */
      fvec_t *vec = fvec_create_row(self->data);
      fvec_set_index(vec, result[0]);
      fts_set_object(ret, vec);
  }
  else
  {   /* return n-tuple of fvec or index */
      fts_tuple_t *tup = (fts_tuple_t *) fts_object_create(fts_tuple_class, 0, NULL);
      fts_atom_t  *at;
      int i;

      fts_tuple_set_size(tup, n);
      at = fts_tuple_get_atoms(tup); 
      
      for (i = 0; i < n; i++)
      {
	  if (s == s_getknni)
	      fts_set_int(at + i, result[i]);
	  else /* getknn */
	  {
	      fvec_t *vec = fvec_create_row(self->data);
	      fvec_set_index(vec, result[i]);
	      
	      fts_object_refer((fts_object_t *) vec);
	      fts_set_object(at + i, vec);
	  }
      }

      fts_set_object(ret, tup);
  }

  return fts_ok;
}


static fts_method_status_t _kdtree_get_node (fts_object_t *o, fts_symbol_t s, 
					   int ac, const fts_atom_t *at, 
					   fts_atom_t *ret)
{
  kdtree_t *t = (kdtree_t *) o;
  int     n = 0;	/* node to inspect */

  if (ac > 1  &&  fts_is_number(at + 1))
  {
      n = fts_get_number_int(at + 1);
      if (n >= t->nnodes)
	  n = t->nnodes - 1;
  }

  if (t->nnodes > 0  &&  ac > 0  &&  fts_is_symbol(at))
  {
      fts_symbol_t  func = fts_get_symbol(at);
      node_t       *node = &t->nodes[n];

      if (func == fts_s_num)
	  fts_set_int(ret, t->nnodes);
      else if (func == s_mean)
	  fts_set_object(ret, node->mean);
      else if (func == s_hyper)
	  fts_set_object(ret, node->split);
      else if (func == fts_s_size)
	  fts_set_int(ret, node->size);
      else if (ac > 2  &&  fts_is_number(at + 2))
      {
	  int i = fts_get_number_int(at + 2); /* node-data vector to inspect */
      
	  if (func == s_vec)
	  {
	      fvec_t *vec = fvec_create_row(t->data);

	      fvec_set_index(vec, t->dataindex[node->startind + i]);
	      fts_set_object(ret, vec);
	  }
	  else if (func == fts_s_index)
	      fts_set_int(ret, t->dataindex[node->startind + i]);
      }
      else if (func == s_profile)
      {
	  fts_post("build profile:\n" 
		   "vector to vector distances:\t%d\n"
		   "vector to node distances:  \t%d\n"
		   "mean vector calculations:  \t%d\n",     
		   "split plane calculations:  \t%d\n",     
		   "searches performed:        \t%d\n",     
		   "neighbours found:          \t%d\n",     
		   t->profile.v2v, t->profile.v2n, 
		   t->profile.mean, t->profile.hyper,
		   t->profile.searches, t->profile.neighbours);
	  fts_set_object(ret, t);
	  profile_clear(t);
      }
  }

  return fts_ok;
}


static fts_method_status_t _kdtree_print (fts_object_t *o, fts_symbol_t s, 
					int ac, const fts_atom_t *at, 
					fts_atom_t *ret)
{
  kdtree_t *self = (kdtree_t *) o;

  kdtree_info_display(self);
  kdtree_data_display(self, 1);

  return fts_ok;
}



/**********************************************************
*
*  class
*
*/

static fts_method_status_t
_kdtree_init(fts_object_t *o, fts_symbol_t s, int ac, const fts_atom_t *at, fts_atom_t *ret)
{
  kdtree_t *self = (kdtree_t *) o;
      
  /* Init */
  self->dmode       = dmode_orthogonal;
  self->mmode       = mmode_mean;
  self->sort        = 1;
  self->ndata  	    = 0;
  self->dataalloc   = 0;
  self->height 	    = 0;
  self->maxheight   = 0;
  self->givenheight = -1;	/* -1 gives less comparisons than -2 */
  self->nnodes 	    = 0;
  self->ndim   	    = 0;
  self->dataindex   = NULL;
  self->nodes  	    = NULL;
  self->data   	    = fmat_create(0, 0);	/* always an fmat here */
  self->sigma  	    = fmat_create(0, 0);	/* always an fmat here */
  self->sigma_nnz   = 0;
  self->sigma_indnz = NULL;

  fts_object_refer((fts_object_t *) self->data);
  fts_object_refer((fts_object_t *) self->sigma);

  kdtree_stack_init(&self->stack, 0);

#if PROFILE_BUILD
  profile_clear(self);
#endif

  if (ac > 0  &&  fts_is_a(at, fmat_class))
  {   /* set data arg and build tree */
      _kdtree_set(o, s, ac, at, ret);
  }

  return fts_ok;
}


static fts_method_status_t
_kdtree_delete(fts_object_t *o, fts_symbol_t s, int ac, const fts_atom_t *at, fts_atom_t *ret)
{
  kdtree_t *self = (kdtree_t *) o;
  int i;
  
  kdtree_clear_nodes(self);

  /* free data */
  fts_object_release(self->data);	/* release stored vectors */

  /* free structure */
  if (self->dataindex)	fts_free(self->dataindex);
  if (self->nodes)	fts_free(self->nodes);
  if (self->sigma_indnz != NULL) fts_free(self->sigma_indnz);

  kdtree_stack_free(&self->stack);

#if DEBUG
  self->data      = NULL;
  self->dataindex = NULL;
  self->nodes     = NULL;
#endif

  return fts_ok;
}


static void
_kdtree_instantiate(fts_class_t *cl)
{
  fts_class_init(cl, sizeof(kdtree_t), _kdtree_init, _kdtree_delete, "[<fmat: data>>] - vector space index tree");
  
/*  
  fts_class_set_copy_function(cl, kdtree_copy_function);
  fts_class_set_array_function(cl, kdtree_array_function);
  fts_class_set_dump_function(cl, kdtree_dump_function);
*/
  fts_class_set_getelem_function(cl, kdtree_getelem_function);

  fts_class_message_varargs(cl, fts_s_print, _kdtree_print, "print list of entries");

/*
  fts_class_message_void(cl, fts_s_clear, _kdtree_clear, "- erase all entries");
  fts_class_message_varargs(cl, fts_new_symbol("rebuild"), _kdtree_rebuild);  
  fts_class_message_varargs(cl, fts_s_remove, _kdtree_remove, "<any: key> ... - remove given entries");
*/

  fts_class_message(cl, fts_s_set,	       fmat_class, _kdtree_set, "<fmat: data> - set data matrix(ndata, ndim)");
/* no:  fts_class_message(cl, fts_new_symbol("add"), fmat_class, _kdtree_add, "<fmat: vector> - add one vector"); */

  fts_class_message_varargs(cl, fts_s_get,                _kdtree_getknn, 
    "<fmat: x> - get frow vector of data matrix nearest to vector x");
  fts_class_message_varargs(cl, fts_new_symbol("getknn"), _kdtree_getknn, 
    "<fmat: x> [<int: k>] - get k nearest neighbour vectors to vector x (default k=1)");
  fts_class_message_varargs(cl, s_getknni,                _kdtree_getknn, 
    "<fmat: x> [<int: k>] - get row indices in data matrix of k nearest neighbours to vector x (default k=1) ");

  /* inspection */ 
 fts_class_message_varargs(cl, fts_new_symbol("node"),  _kdtree_get_node,  
   "<num|mean|hyper|size|start|end|profile> <int: index> - get information about a node");

}


FTS_MODULE_INIT(kdtree)
{
  kdtree_symbol = fts_new_symbol("kdtree");
  s_mean  = fts_new_symbol("mean");
  s_hyper = fts_new_symbol("hyper");
  s_vec   = fts_new_symbol("vec");
  s_profile = fts_new_symbol("profile");
  s_getknni = fts_new_symbol("getknni");

  kdtree_class = fts_class_install(kdtree_symbol, _kdtree_instantiate);
}
