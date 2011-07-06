#include "rta_unispring.h"
#define DIM RTA_UNISPRING_NDIM

using namespace UniSpringSpace ;

/**
 To launch when setting up Qhull for the first time.
 */
void UniSpring::setupQhull(){
	
	rows = new coordT*[mNpoints];
	ismalloc= False;    // True if qhull should free points in qh_freeqhull() or reallocation
	flags = new char[250];          // option flags for qhull, see qh_opt.htm
	outfile = NULL;    // output from qh_produce_output(). use NULL to skip qh_produce_output(). Otherwise use stdout
	errfile= stderr;    // error messages from qhull code
	sprintf (flags, "qhull d Qt Qbb Qc %s", ""); // Options used in Distmesh	
	
#if qh_QHpointer  // see user.h (Qhull distribution)
	if (qh_qh){
		printf ("QH6233: Qhull link error.  The global variable qh_qh was not initialized\n\
				to NULL by global.c.  Please compile user_eg.c with -Dqh_QHpointer_dllimport\n\
				as well as -Dqh_QHpointer, or use libqhullstatic, or use a different tool chain.\n\n");
	}
#endif
	
}


/**
 To launch when performing the first delaunay triangulation.
 */
void UniSpring::triangulate(){
	
	int i;
	for (i=mNpoints; i--; )
		rows[i]= mPoints+DIM*i;
	//qh_printmatrix (outfile, "input", rows, mNpoints, DIM); // Unused, print input points
	exitcode= qh_new_qhull (DIM, mNpoints, mPoints, ismalloc,
							flags, outfile, errfile);
	
}


/**
 To launch when performing subsequent delaunay triangulations.
 */
void UniSpring::retriangulate(){
	
	int i;
	for (i=mNpoints; i--; )
		rows[i]= mPoints+DIM*i;
	//oldqh = qh_save_qhull(); // unused, provided for example: to save intermediate complex hull
	exitcode= qh_new_qhull2 (DIM, mNpoints, mPoints, ismalloc,
							 flags, outfile, errfile);
	
	//printf("Retriangulated.\n");
	
}


/**
 Print a summary of Qhull's output.
 */
void UniSpring::print_summary (void) {
	facetT *facet;
	int k;
	
	printf ("\n%d vertices and %d facets with normals:\n",
			qh num_vertices, qh num_facets);
	FORALLfacets {
		for (k=0; k < qh hull_dim; k++)
			printf ("%6.2g ", facet->normal[k]);
		printf ("\n");
	}
}


/**
 Get delaunay triangulation edges by visiting each facet ridge once.
 Only edges belonging to a facet whose centroid is inside the target region boundaries are kept.
 */
void UniSpring::getEdgeVector(){
	
	if (!exitcode) { // if no error
		//'qh facet_list' contains the convex hull
				
		if (!qh_qh){
			printf ("qh_qh non initialisé");
		}
		
		qh visit_id++;
		
		FORALLfacets {
			
			if (!facet->upperdelaunay) {
				
				facet->visitid= qh visit_id;
				qh_makeridges(facet);
				std::vector< std::vector<int> > Edges_temp; // Current facet edges
				std::vector<double> centroid(2); // Centroid coordinates
				
				// Check if centroid facet is inside
				FOREACHridge_(facet->ridges) {
					
					std::vector<int> edge_temp;
					
					FOREACHvertex_(ridge->vertices) {
						
						//printf("id : %d ",qh_pointid (vertex->point) ); // obention de l'ID du point
						edge_temp.push_back(qh_pointid (vertex->point));
						
					}
					
					Edges_temp.push_back(edge_temp);
					
				}
				
				// Compute centroid
				centroid[0] = (mPoints[Edges_temp[0][0]*DIM] + mPoints[Edges_temp[0][1]*DIM] + mPoints[Edges_temp[1][0]*DIM] + mPoints[Edges_temp[1][1]*DIM] + mPoints[Edges_temp[2][0]*DIM] + mPoints[Edges_temp[2][1]*DIM])/6; //Check
				centroid[1] = (mPoints[Edges_temp[0][0]*DIM+1] + mPoints[Edges_temp[0][1]*DIM+1] + mPoints[Edges_temp[1][0]*DIM+1] + mPoints[Edges_temp[1][1]*DIM+1] + mPoints[Edges_temp[2][0]*DIM+1] + mPoints[Edges_temp[2][1]*DIM+1])/6;
				
				// if centroid is inside, record vertex indices of current facet				
				if (mShape->fd_compute2(centroid[0],centroid[1])<-GEPS) {
					
					FOREACHridge_(facet->ridges) {
						
						std::vector<int> edge_temp;
						neighbor= otherfacet_(ridge, facet);
						
						if (neighbor->visitid != qh visit_id) { // Check if ridge has already been visited (in neighbor facet)
							
							FOREACHvertex_(ridge->vertices) {
								
								edge_temp.push_back(qh_pointid (vertex->point));  
								
							}
							
							mEdges.push_back(edge_temp);
							
						}
						
					}
					
				}
				
			}
			
		}
		
	}
		
}


/**
 Used in UniSpring::retriangulate. Use previously allocate qh_qh pointer.
 */
 int UniSpring::qh_new_qhull2(int dim, int numpoints, coordT *points, boolT ismalloc,
						   char *qhull_cmd, FILE *outfile, FILE *errfile) {
	int exitcode, hulldim;
	boolT new_ismalloc;
	static boolT firstcall = True;
	coordT *new_points;
	
	if (firstcall) {
		qh_meminit(errfile);
		firstcall= False;
	}
	if (strncmp(qhull_cmd,"qhull ", (size_t)6)) {
		qh_fprintf(errfile, 6186, "qhull error (qh_new_qhull): start qhull_cmd argument with \"qhull \"\n");
		qh_exit(qh_ERRinput);
	}
	qh_initstatistics();
	qh_initqhull_start2(NULL, outfile, errfile); // Does not reallocate qh_qh
	trace1((qh ferr, 1044, "qh_new_qhull: build new Qhull for %d %d-d points with %s\n", numpoints, dim, qhull_cmd));
	exitcode = setjmp(qh errexit);
	if (!exitcode)
	{
		qh NOerrexit = False;
		qh_initflags(qhull_cmd);
		if (qh DELAUNAY)
			qh PROJECTdelaunay= True;
		if (qh HALFspace) { // Points is an array of halfspaces, the last coordinate of each halfspace is its offset
			hulldim= dim-1;
			qh_setfeasible(hulldim);
			new_points= qh_sethalfspace_all(dim, numpoints, points, qh feasible_point);
			new_ismalloc= True;
			if (ismalloc)
				qh_free(points);
		}else {
			hulldim= dim;
			new_points= points;
			new_ismalloc= ismalloc;
		}
		qh_init_B(new_points, numpoints, hulldim, new_ismalloc);
		qh_qhull();
		qh_check_output();
		if (outfile)
			qh_produce_output();
		if (qh VERIFYoutput && !qh STOPpoint && !qh STOPcone)
			qh_check_points();
	}
	qh NOerrexit = True;
	return exitcode;
}


/**
 Free memory allocate by UniSpring::qh_new_qhull2.
 */
void UniSpring::freeQhullMemory(){
	
	qh_freeqhull2(!qh_ALL);                 // Free long memory
	qh_memfreeshort (&curlong, &totlong);  // Free short memory and memory allocator
	if (curlong || totlong)
		fprintf (errfile, "qhull internal warning (user_eg, #2): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
	
}
