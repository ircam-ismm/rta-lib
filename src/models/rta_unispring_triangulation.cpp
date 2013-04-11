#include "rta_unispring.h"
#define DIM RTA_UNISPRING_NDIM

using namespace UniSpringSpace ;


/**
 To launch when performing the first delaunay triangulation.
 */
/*
void UniSpring::triangulate(){
	
	int i;
	for (i=mNpoints; i--; )
		rows[i]= mPoints+DIM*i;
	//qh_printmatrix (outfile, "input", rows, mNpoints, DIM); // Unused, print input points
	exitcode= qh_new_qhull (DIM, mNpoints, mPoints, ismalloc,
							flags, outfile, errfile);
	
}*/


/**
 To launch when performing subsequent delaunay triangulations.
 */
/*
void UniSpring::retriangulate(){
	
	int i;
	for (i=mNpoints; i--; )
		rows[i]= mPoints+DIM*i;
	//oldqh = qh_save_qhull(); // unused, provided for example: to save intermediate complex hull
	exitcode= qh_new_qhull2 (DIM, mNpoints, mPoints, ismalloc,
							 flags, outfile, errfile);
	
	//printf("Retriangulated.\n");
	
}*/



/**
 Get delaunay triangulation edges by visiting each facet ridge once.
 Only edges belonging to a facet whose centroid is inside the target region boundaries are kept.
 */
void UniSpring::getEdgeVector(){
	
	const list<hed::Edge*>& leadingEdges = triang.getLeadingEdges();
	list<hed::Edge*>::const_iterator it2;
	
	mEdges.clear(); // Reset		

	// iterate over all triangles
	for (it2 = leadingEdges.begin(); it2 != leadingEdges.end(); ++it2) {
		hed::Edge* edge = *it2;
		
		std::vector<int> tmp_ids(3,0);
		
		// get all nodes in triangle
		for (int i = 0; i < 3; ++i) {
			
			hed::Node* node = edge->getSourceNode();
			tmp_ids[i] = node->id();
					
			edge = edge->getNextEdgeInFace();
			
		}
		
		std::vector<int> edge_temp1(2);
		std::vector<int> edge_temp2(2);
		std::vector<int> edge_temp3(2);
		
		edge_temp1[0] = tmp_ids[0];
		edge_temp1[1] = tmp_ids[1];
		edge_temp2[0] = tmp_ids[0];
		edge_temp2[1] = tmp_ids[2];
		edge_temp3[0] = tmp_ids[1];
		edge_temp3[1] = tmp_ids[2];
		
		// put smaller indices first (needed by removeDuplicateEdges)
		std::sort(edge_temp1.begin(), edge_temp1.end());
		std::sort(edge_temp2.begin(), edge_temp2.end());
		std::sort(edge_temp3.begin(), edge_temp3.end());
		
		if (euclDistance(edge_temp1[0], edge_temp1[1]) < MAX_EDGE_LENGTH) {
			mEdges.push_back(edge_temp1);
		}
		if (euclDistance(edge_temp2[0], edge_temp2[1]) < MAX_EDGE_LENGTH) {
			mEdges.push_back(edge_temp2);
		}
		if (euclDistance(edge_temp3[0], edge_temp3[1]) < MAX_EDGE_LENGTH) {
			mEdges.push_back(edge_temp3);
		}
		
	}
	
	removeDuplicateEdges();
		
}

void UniSpring::getEdgeVector_3D(){ // TODO: check if it works in 3D
	
	const list<hed::Edge*>& leadingEdges = triang.getLeadingEdges();
	list<hed::Edge*>::const_iterator it2;
	
	mEdges.clear(); // Reset		

	// iterate over all triangles
	for (it2 = leadingEdges.begin(); it2 != leadingEdges.end(); ++it2) {
		hed::Edge* edge = *it2;
		
		std::vector<int> tmp_ids(3,0);
		
		// get all nodes in triangle
		for (int i = 0; i < 3; ++i) {
			
			hed::Node* node = edge->getSourceNode();
			tmp_ids[i] = node->id();
			edge = edge->getNextEdgeInFace();
			
		}
		
		// Compute centroid
		std::vector<double> centroid(2,0);
		
		centroid[0] = (mPoints[tmp_ids[0]].x() + mPoints[tmp_ids[1]].x() + mPoints[tmp_ids[2]].x())/3;
		centroid[1] = (mPoints[tmp_ids[0]].y() + mPoints[tmp_ids[1]].y() + mPoints[tmp_ids[2]].y())/3;

		if (mShape->fd_compute(centroid[0],centroid[1])<-GEPS) {
			
			std::vector<int> edge_temp1(2);
			std::vector<int> edge_temp2(2);
			std::vector<int> edge_temp3(2);
			
			edge_temp1[0] = tmp_ids[0];
			edge_temp1[1] = tmp_ids[1];
			edge_temp2[0] = tmp_ids[0];
			edge_temp2[1] = tmp_ids[2];
			edge_temp3[0] = tmp_ids[1];
			edge_temp3[1] = tmp_ids[2];
			
			// put smaller indices first
			std::sort(edge_temp1.begin(), edge_temp1.end());
			std::sort(edge_temp2.begin(), edge_temp2.end());
			std::sort(edge_temp3.begin(), edge_temp3.end());

			if (euclDistance(edge_temp1[0], edge_temp1[1]) < MAX_EDGE_LENGTH) {
				mEdges.push_back(edge_temp1);
			}
			if (euclDistance(edge_temp2[0], edge_temp2[1]) < MAX_EDGE_LENGTH) {
				mEdges.push_back(edge_temp2);
			}
			if (euclDistance(edge_temp3[0], edge_temp3[1]) < MAX_EDGE_LENGTH) {
				mEdges.push_back(edge_temp3);
			}
			
		}
		
	}
	
	//removeDuplicateEdges();
	
	
}

void UniSpring::removeDuplicateEdges() {
	
	std::sort(mEdges.begin(), mEdges.end()); //sort uses operator<, overloaded below
	mEdges.erase(std::unique(mEdges.begin(), mEdges.end()), mEdges.end()); // std::unique uses std::adjacent_find which in turn uses operator==, overloaded below	
}


bool operator<(std::vector<int> const& v1, std::vector<int> const& v2) {
	
	if (v1[0] < v2[0]) return true;
	else return false;
	
}

bool operator==(std::vector<int> const& v1, std::vector<int> const& v2) { // assume that vectors have 2 elements (edge indices)
	
	if (v1[0] == v2[0] && v1[1] == v2[1]) return true;
	else return false;
	
}
