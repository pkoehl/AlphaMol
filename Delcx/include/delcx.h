/* ====================================================================
   delcx
 
	This class computes the regular triangulation of a set
	of N weighted points in 3D, using the incremental flipping
	algorithm of Edelsbrunner.

	This implementation is based on the algorithm published in
	H. Edelsbrunner and N.R. Shah, Algorithmica (1996) 15: 223-241

	1) Algorithm:
	*************

	Briefly, the algorithm works as follows:

	- first, a large tetrahedron initialises the program. All
	four vertices of this tetrahedron are set at "infinite"

	- All N points are added one by one.

	- For each point:

		- localize the tetrahedron in the current regular
		triangulation that contains this point

		- test if the point is redundant; if yes, remove

		- If the point is not redundant, insert in the
		tetrahedron : this is a "1-4" flip

		- collect all "link facets" (i.e. all triangles
		in tetrahedron containing the new point, that face
		this new point) that are not regular.

		- for each non-regular link facet, check if it
		is "flippable". If yes, perform a "2-3", "3-2"
		or "1-4" flip. Add new link facets in the list,
		if needed.

		- when link facet list if empty, move to next 
		point

	- Remove "infinite" tetrahedron, i.e. tetrahedron with
	one vertice at "infinite"

	- collect all remaining tetrahedron, define convex hull,
	and exit.

	2) Data structure:
	******************

	I maintain a minimal data structure that includes only 
	the tetrahedrons of the triangulation (triangles,
	edges and vertices are implicit).

	For each tetrahedron, I store:

	- the index of its four vertices
	- pointers to its neighbours (4 maximum).
		neighbor(i) is the tetrahedron that shares
		all vertices of the tetrahedron considered, except i
		(0 if the corresponding face is on the convex hull)
	- its status: 1 "active" (i.e. part of the triangulation), 0 inactive
	- its orientation


	2) number representation:
	**************************

	I use double precision floating points. However, if one of
	the geometricaly test becomes "imprecise", I switch to
	arbitrary precision arithmetics (using the gmp package).

	All predicates have therefore a floating point filter

 ==================================================================== */

#ifndef DELCX_H
#define DELCX_H

#include <vector>
#include <queue> 
#include <algorithm>
#include "Vertex.h"
#include "Vector.h"
#include "Tetrahedron.h"
#include "sos_gmp.h"
#include "sort_tools.h"

  SOS sos;
  SortingTools mysort;

/*********************************************************************************
  Delcx class
 *********************************************************************************/

  class DELCX {

  public:
	// Setup calculations
	void setup(int npoints, double *coord, double *radii, double *coefS, double *coefV,
	double *coefM, double *coefG, std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra);

	// Compute 3D weighted Delaunay triangulation
	void regular3D(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra);

	// generate list of edges
	void delaunayEdges(std::vector<Tetrahedron>& tetra, std::vector<std::pair<int, int> >& edges);

  private:

	// locate tetrahedron in which point is inserted
	void locate_jw(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra, 
			int ipoint, int *tetra_loc, int *iredundant);

	// go over full link_facet
	void flip(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra);

	// Check if a point is inside a tetrahedron
	void inside_tetra(std::vector<Vertex>& vertices, int p, int a, int b, int c, int d, 
		int iorient, bool *is_in, bool *redundant, int *ifail);
 
	// Check if a facet connects two tetrehedra that are convex
	void regular_convex(std::vector<Vertex>& vertices, int a, int b, int c, 
		int p, int o, int itest_abcp, bool *regular, bool *convex, bool *test_abpo, 
		bool *test_bcpo, bool *test_capo);

	// sign associated with the missing inf point
	void missinf_sign(int i, int j, int k, int *l, int *sign);

	// flip 1-4: from 1 to 4 tetrahedra
	void flip_1_4(std::vector<Tetrahedron>& tetra, int ipoint, int itetra, 
		int *tetra_last);

	// flip 4-1: from 4 to 1 tetrahedron
	void flip_4_1(std::vector<Vertex>& vertice, std::vector<Tetrahedron>& tetra, 
	int itetra, int jtetra, int ktetra, int ltetra, int *vertices,
	int idp, int jdp, int kdp, int ldp, bool test_acpo, int *ierr, int *tetra_last);

	// flip 2-3: from 2 to 3 tetrahedra
	void flip_2_3(std::vector<Tetrahedron>& tetra, int itetra, int jtetra, 
	int *vertices, int *facei, int *facej, bool test_abpo, bool test_bcpo, 
	bool test_capo, int *ierr, int *tetra_last);

	// flip 3-2: from 3 to 2 tetrahedra
	void flip_3_2(std::vector<Tetrahedron>& tetra, int itetra, int jtetra, 
		int ktetra, int *vertices, int *edgei, int *edgej, int *edgek,
		bool test_bcpo, bool test_acpo, int *ierr, int *tetra_last);

	// Define facet between two tetrahedron
	void define_facet(std::vector<Tetrahedron>& tetra, int itetra, int jtetra, 
		int idx_o, int *facei, int *facej);

	// info about tetrahedron
	void find_tetra(std::vector<Tetrahedron>& tetra, int itetra, int idx_c, 
	int a, int b, int o, int *ifind, int *tetra_loc, int *idx_a, int *idx_b);

	// reorder vertices of tetrahedra in increasing order
	void reorder_tetra(std::vector<Tetrahedron>& tetra);

	// remove "infinite" tetrahedron
	void remove_inf(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra);

	// mark tetrahedron
	void mark_zero(std::vector<Tetrahedron>& tetra, int itetra, int ivertex);

	// remove flat tetrahedra
	void peel(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra);

	// Computes the volume of a tetrahedron
	double tetra_vol(double *a, double *b, double *c, double *d);

	// find edge in tetra defined by 2 vertices
	int findEdge(Tetrahedron t, int i1, int i2);

	// Add bogus points as needed so that we have at least 4 points
	void addBogus(int npoints, double *coord, double *radii, double *bcoord, double *brad); 

  protected:

	std::queue<std::pair<int, int> > link_facet;
	std::queue<std::pair<int, int> > link_index;
	std::stack<int> free;
	std::vector<int> kill;

	double eps = 1.e-5;

	int inf4_1[4] = {1, 1, 0, 0};
	int sign4_1[4] = {-1, 1, 1, -1};
	int inf4_2[4][4] = {
		{ -1, 1, 2, 2},
		{ 1, -1, 2, 2},
		{ 2, 2, -1, 0},
		{ 2, 2, 0, -1}
	};
	int sign4_2[4][4] = {
		{ 0, 1, -1, 1},
		{ -1, 0, 1, -1},
		{ 1, -1, 0, 1},
		{ -1, 1, -1, 0}
	};
	int sign4_3[4] = {-1, 1, -1, 1};
	int inf5_2[4][4] = {
		{ -1, 1, 0, 0},
		{ 1, -1, 0, 0},
		{ 0, 0, -1, 0},
		{ 0, 0, 0, -1}
	};
	int sign5_2[4][4] = {
		{ 0, -1, -1, 1},
		{ 1, 0, -1, 1},
		{ 1, 1, 0, 1},
		{ -1, -1, -1, 0}
	};
	int inf5_3[4] = {0, 0, 2, 2};
	int sign5_3[4] = {1, 1, -1, 1};
	int order1[4][3] = {
		{ 2, 1, 3},
		{ 0, 2, 3},
		{ 1, 0, 3},
		{ 0, 1, 2}
	};

	int ord_rc[3][3] = {
		{0, 1, 2},
		{2, 0, 1},
		{1, 2, 0},
	};

	int order2[6][2] = {
		{ 2, 3},
		{ 3, 1},
		{ 1, 2},
		{ 0, 3},
		{ 2, 0},
		{ 0, 1}
	};
	int order3[6][2] = {
		{ 0, 1},
		{ 0, 2},
		{ 0, 3},
		{ 1, 2},
		{ 1, 3},
		{ 2, 3}
	};
	int idxList[4][3] = {
		{ 0, 0, 0},
		{ 0, 1, 1},
		{ 1, 1, 2},
		{ 2, 2, 2}
	};
	int table32[3][3] = {
		{ 0, 1, 2},
		{ 0, 2, 1},
		{ 2, 0, 1}
	};
	int table32_2[3][2] = {
		{ 0, 1},
		{ 0, 2},
		{ 1, 2}
	};
	int table41[3][3] = {
		{ 1, 0, 2},
		{ 0, 1, 2},
		{ 0, 2, 1}
	};
	int table41_2[3][2] = {
		{ 0, 0},
		{ 1, 0},
		{ 1, 1}
	};
	int order[3][2] = {
		{ 1, 2},
		{ 2, 0},
		{ 0, 1}
	};
	int other[4][3] = {
		{ 1, 2, 3},
		{ 0, 2, 3},
		{ 0, 1, 3},
		{ 0, 1, 2}
	};
	int other2[4][4][2] = {
		{
			{ -1, -1},
			{ 2, 3},
			{ 1, 3},
			{ 1, 2}
		},
		{
			{ 2, 3},
			{ -1, -1},
			{ 0, 3},
			{ 0, 2}
		},
		{
			{ 1, 3},
			{ 0, 3},
			{ -1, -1},
			{ 0, 1}
		},
		{
			{ 1, 2},
			{ 0, 2},
			{ 0, 1},
			{ -1, -1}
		},
	};

  };

/* ====================================================================
	This program gets the coordinates of the N points considered, and stores
	these into the structures (common blocks) used in all suite of programs
	Regular3D.f
 ==================================================================== */

  void DELCX::setup(int npoints, double *coord, double *radii, double *coefS, double *coefV,
	double *coefM, double *coefG, std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra)
  {

	sos.init_sos_gmp();

/* ====================================================================
	Initialisation:
		- clear up all data structures
 ==================================================================== */

	vertices.clear();
	tetra.clear();

	while (!link_facet.empty())
	{
		link_facet.pop();
	}

	while (!link_index.empty())
	{
		link_facet.pop();
	}

	while (!free.empty())
	{
		free.pop();
	}

	kill.clear();

/* ====================================================================
	Initialisation:
		- set four "infinite" points
		- copy atoms into vertex list
 ==================================================================== */

	double zero=0.;
	for(int i = 0; i < 4; i++) {
		Vertex vert(zero, zero, zero, zero, zero, zero, zero, zero);
		vert.info[0] = 1;
		vert.status = 0;
		vertices.push_back(vert);
	}

/* ====================================================================
	Add remaining points
 ==================================================================== */

	double x, y, z, r, c, d, e, f;
	for(int i = 0; i < npoints; i++) {
		x = coord[3*i]; y = coord[3*i+1]; z = coord[3*i+2];
		r = radii[i]; c = coefS[i];
		d = coefV[i]; e = coefM[i]; f = coefG[i];
		Vertex vert(x, y, z, r, c, d, e, f);
		vert.info[0] = 1;
		vert.status = 1;
		vertices.push_back(vert);
	}

/* ====================================================================
	If number of points is smaller than 4, add "bogus" points
 ==================================================================== */

	if(npoints < 4) {
		int new_points = 4-npoints;
		double *bcoord = new double[3*new_points];
		double *brad   = new double[new_points];
		addBogus(npoints, coord, radii, bcoord, brad); 
		for(int i = 0; i < new_points; i++) {
			x = bcoord[3*i]; y = bcoord[3*i+1]; z = bcoord[3*i+2];
			r = brad[i]; c = 1.; d = 1; e = 1; f = 1;
			Vertex vert(x, y, z, r, c, d, e, f);
			vert.info[0] = 1;
			vert.status = 0;
			vertices.push_back(vert);
		}
	}

/* ====================================================================
	Initialisation:
		- create an "infinite" tetrahedron
 ==================================================================== */

	Tetrahedron t;
	t.init();

	t.Vertices[0] = 0;
	t.Vertices[1] = 1;
	t.Vertices[2] = 2;
	t.Vertices[3] = 3;

	t.Neighbours[0] = -1;
	t.Neighbours[1] = -1;
	t.Neighbours[2] = -1;
	t.Neighbours[3] = -1;

	t.nindex[0] = -1;
	t.nindex[1] = -1;
	t.nindex[2] = -1;
	t.nindex[3] = -1;

	t.info[0] = 0;
	t.info[1] = 1;

	tetra.push_back(t);

  }

/* ====================================================================
	Implements a 3D weighted delaunay triangulation

	Input:
	*******

	vertices:	vertices (with vertices 0-3 representing
			an infinite tetrahedron
	Output:
	********

	tetra:		all tetrahedra of the regular triangulation
 ==================================================================== */

  void DELCX::regular3D(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra)
  {

	int iredundant, tetra_loc;
	int ipoint;

/* ====================================================================
	Build regular triangulation
 ==================================================================== */

	int tetra_last = -1;
	int npoint = vertices.size()-4;

	for(int i = 0;  i < npoint; i++)
	{
		ipoint = i+4;
		if(vertices[ipoint].info[1] == 0) continue;

/* ====================================================================
		first locate the point in the list of known tetrahedra
 ==================================================================== */

		tetra_loc = tetra_last;
		locate_jw(vertices, tetra, ipoint, &tetra_loc, &iredundant);

/* ====================================================================
		If the point is redundant, move to next point
 ==================================================================== */

		if(iredundant==1) {
			vertices[ipoint].info[0] = 0;
			continue;
		}

/* ====================================================================
		Otherwise, add point to tetrahedron : 1-4 flip
 ==================================================================== */

		int tetra_last;
		flip_1_4(tetra, ipoint, tetra_loc, &tetra_last);

/* ====================================================================
		Now scan link_facet list, and flip till list is empty
 ==================================================================== */

		flip(vertices, tetra);

/* ====================================================================
		At this stage, I should have a regular triangulation
		of the i+4 points (i-th real points+4 "infinite" points)
		Now add another point
 ==================================================================== */
		
	}

/* ====================================================================
	Reorder the tetrahedra, such that vertices are in increasing order
 ==================================================================== */

	reorder_tetra(tetra);

/* ====================================================================
	I have the regular triangulation: I need to remove the
	simplices including infinite points, and define the
	convex hull
 ==================================================================== */

	remove_inf(vertices, tetra);

/* ====================================================================
	Now I peel off flat tetrahedra at the boundary of the DT
 ==================================================================== */

//	peel(vertices, tetra);

	sos.clear_sos_gmp();
  }

/* =====================================================================================================

 This procedure tests if a point p is inside a tetrahedron defined by four points (a,b,c,d) 
 (with orientation "iorient"). If p is found inside the tetrahedron, it also checks if it
 is redundant

 Computations are performed  in floating point, but they are switched to multiple precision if the 
 results are imprecise

 ===================================================================================================== */

  void DELCX::inside_tetra(std::vector<Vertex>& vertices, int p, int a, int b, int c, int d, 
	int iorient, bool *is_in, bool *redundant, int *ifail)
  {

	int i, j, k, l;
	int ia = 0, ib = 0, ic = 0, id = 0, ie = 0, idx;
	int ic1, ic5, ic1_k, ic1_l, sign, sign5, sign_k, sign_l;

	int nswap = 0, iswap = 0, ninf;
	int val = 0;

	int list[4];
	int infPoint[4];

	double xa, ya, xb, yb, xc, yc;
	double ra, rb, rc, rd, re;
	double coord_a[3], coord_b[3], coord_c[3], coord_d[3], coord_e[3];

	double Sij_1, Sij_2, Sij_3, Skl_1, Skl_2, Skl_3;
	double det_pijk, det_pjil, det_pkjl, det_pikl, det_pijkl;

	double detij[3];
	double i_p[4];
	double j_p[4];
	double k_p[4];
	double l_p[4];
	bool test_pijk, test_pjil, test_pkjl, test_pikl;

	/*
	If (i,j,k,l) is the tetrahedron in positive orientation, we need to test:
	(p,i,j,k)
	(p,j,i,l)
	(p,k,j,l)
	(p,i,k,l)
	If all four are positive, than p is inside the tetrahedron.
	All four tests relies on the sign of the corresponding 4x4
	determinant. Interestingly, these four determinants share
	some common lines, which can be used to speed up the computation.

	Let us consider or example:

	det(p,i,j,k) =  | p(1) p(2) p(3) 1|
			| i(1) i(2) i(3) 1|
			| j(1) j(2) j(3) 1|
			| k(1) k(2) k(3) 1|

	p appears in each determinant. The corresponding line can therefore
	be substraced from all 3 other lines . Using the example above,
	we find:

	det(i,j,k,l) = -|ip(1) ip(2) ip(3)|
			|jp(1) jp(2) jp(3)|
			|kp(1) kp(2) kp(3)|

	where :xp(m) = x(m) - p(m) for x = i,j,k and m = 1,2,3

	Now we notice that the first two lines of det(p,i,j,k) and det(p,i,j,l) are the same.

	Let us define: Sij_3 = |ip(1) ip(2)| Sij_2 = |ip(1) ip(3)| and Sij_1 = |ip(2) ip(3)|
			   	|jp(1) jp(2)|         |jp(1) jp(3)|             |jp(2) jp(3)|

	We find:
	det(p,i,j,k) = - kp(1)*Sij_1 + kp(2)*Sij_2 - kp(3)*Sij_3
	and:
	det(p,j,i,l) =   lp(1)*Sij_1 - lp(2)*Sij_2 + lp(3)*Sij_3

	Similarly, if we define:

	Skl_3 = |kp(1) kp(2)|	Skl_2 = |kp(1) kp(3)|	Skl_1 = |kp(2) kp(3)|
		|lp(1) lp(2)|		|lp(1) lp(3)|		|lp(2) lp(3)|

	We find:
	det(p,k,j,l) = jp(1)*Skl_1 - jp(2)*Skl_2 + jp(3)*Skl_3
	and:
	det(p,i,k,l) = -ip(1)*Skl_1 + ip(2)*Skl_2 - ip(3)*Skl_3

	Furthermore:
	det(p,i,j,k,l) = -ip(4)*det(p,k,j,l)-jp(4)*det(p,i,k,l)
			-kp(4)*det(p,j,i,l)-lp(4)*det(p,i,j,k)

	The equations above hold for the general case; special care is
	required to take in account infinite points (see below)
	 */

	*is_in = false;
	*redundant = false;

	list[0] = a;
	list[1] = b;
	list[2] = c;
	list[3] = d;

	memset(infPoint, 0, 4*sizeof(int));
	if(a<4) infPoint[0] = 1;
	if(b<4) infPoint[1] = 1;
	if(c<4) infPoint[2] = 1;
	if(d<4) infPoint[3] = 1;

	ninf = infPoint[0] + infPoint[1] + infPoint[2] + infPoint[3];

	if (ninf == 0) //no infinite points
	{
		for (i = 0; i < 3; i++) {
			i_p[i] = vertices[a].Coordinates[i] - vertices[p].Coordinates[i];
			j_p[i] = vertices[b].Coordinates[i] - vertices[p].Coordinates[i];
			k_p[i] = vertices[c].Coordinates[i] - vertices[p].Coordinates[i];
			l_p[i] = vertices[d].Coordinates[i] - vertices[p].Coordinates[i];
		}

		//Now compute 2x2 determinants Sij and Skl
		Sij_1 = i_p[1] * j_p[2] - i_p[2] * j_p[1];
		Sij_2 = i_p[0] * j_p[2] - i_p[2] * j_p[0];
		Sij_3 = i_p[0] * j_p[1] - i_p[1] * j_p[0];

		Skl_1 = k_p[1] * l_p[2] - k_p[2] * l_p[1];
		Skl_2 = k_p[0] * l_p[2] - k_p[2] * l_p[0];
		Skl_3 = k_p[0] * l_p[1] - k_p[1] * l_p[0];

		//Now perform tests

		//We check all other four determinants
		det_pijk = -k_p[0] * Sij_1 + k_p[1] * Sij_2 - k_p[2] * Sij_3;
		det_pijk = det_pijk * iorient;
		test_pijk = (std::abs(det_pijk) > eps);
		if (test_pijk && det_pijk > 0) {
			*ifail = 3;
			return;
		}

		det_pjil = l_p[0] * Sij_1 - l_p[1] * Sij_2 + l_p[2] * Sij_3;
		det_pjil = det_pjil * iorient;
		test_pjil = (std::abs(det_pjil) > eps);
		if (test_pjil && det_pjil > 0) {
			*ifail = 2;
			return;
		}

		det_pkjl = j_p[0] * Skl_1 - j_p[1] * Skl_2 + j_p[2] * Skl_3;
		det_pkjl = det_pkjl * iorient;
		test_pkjl = (std::abs(det_pkjl) > eps);
		if (test_pkjl && det_pkjl > 0) {
			*ifail = 0;
			return;
		}

		det_pikl = -i_p[0] * Skl_1 + i_p[1] * Skl_2 - i_p[2] * Skl_3;
		det_pikl = det_pikl * iorient;
		test_pikl = (std::abs(det_pikl) > eps);
		if (test_pikl && det_pikl > 0) {
			*ifail = 1;
			return;
		}
		/*
		At this stage, either all four determinants are positive,
		or one of the determinant is not precise enough, and
		we need to switch the MP
		In this case, since we may need SoS, we have to rank
		the indices
		 */
		if (!test_pijk) {
			mysort.valsort4(p, a, b, c, &ia, &ib, &ic, &id, &nswap);
			for(int i = 0; i < 3; i++) {
				coord_a[i] = vertices[ia].Coordinates[i];
				coord_b[i] = vertices[ib].Coordinates[i];
				coord_c[i] = vertices[ic].Coordinates[i];
				coord_d[i] = vertices[id].Coordinates[i];
			}
			sos.sos_minor4_gmp(coord_a, coord_b, coord_c, coord_d, &val);
			val = val * nswap * iorient;
			if (val == 1) {
				*ifail = 3;
				return;
			}
		}

		if (!test_pjil) {
			mysort.valsort4(p, b, a, d, &ia, &ib, &ic, &id, &nswap);
			for(int i = 0; i < 3; i++) {
				coord_a[i] = vertices[ia].Coordinates[i];
				coord_b[i] = vertices[ib].Coordinates[i];
				coord_c[i] = vertices[ic].Coordinates[i];
				coord_d[i] = vertices[id].Coordinates[i];
			}
			sos.sos_minor4_gmp(coord_a, coord_b, coord_c, coord_d, &val);
			val = val * nswap * iorient;
			if (val == 1) {
				*ifail = 2;
				return;
			}
		}

		if (!test_pkjl) {
			mysort.valsort4(p, c, b, d, &ia, &ib, &ic, &id, &nswap);
			for(int i = 0; i < 3; i++) {
				coord_a[i] = vertices[ia].Coordinates[i];
				coord_b[i] = vertices[ib].Coordinates[i];
				coord_c[i] = vertices[ic].Coordinates[i];
				coord_d[i] = vertices[id].Coordinates[i];
			}
			sos.sos_minor4_gmp(coord_a, coord_b, coord_c, coord_d, &val);
			val = val * nswap * iorient;
			if (val == 1) {
				*ifail = 0;
				return;
			}
		}

		if (!test_pikl) {
			mysort.valsort4(p, a, c, d, &ia, &ib, &ic, &id, &nswap);
			for(int i = 0; i < 3; i++) {
				coord_a[i] = vertices[ia].Coordinates[i];
				coord_b[i] = vertices[ib].Coordinates[i];
				coord_c[i] = vertices[ic].Coordinates[i];
				coord_d[i] = vertices[id].Coordinates[i];
			}
			sos.sos_minor4_gmp(coord_a, coord_b, coord_c, coord_d, &val);
			val = val * nswap * iorient;
			if (val == 1) {
				*ifail = 1;
				return;
			}
		}

		//If we have gone that far, p is inside the tetrahedron
		*is_in = 1; //true

		//Now we check if p is redundant
		i_p[3] = vertices[a].Weight - vertices[p].Weight;
		j_p[3] = vertices[b].Weight - vertices[p].Weight;
		k_p[3] = vertices[c].Weight - vertices[p].Weight;
		l_p[3] = vertices[d].Weight - vertices[p].Weight;

		det_pijkl = -i_p[3] * det_pkjl - j_p[3] * det_pikl - k_p[3] * det_pjil - l_p[3] * det_pijk;

		//No need to multiply by iorient, since all minors contais iorient
		if (std::abs(det_pijkl) < eps) {
			mysort.valsort5(p, a, b, c, d, &ia, &ib, &ic, &id, &ie, &nswap);
			for(int i = 0; i < 3; i++) {
				coord_a[i] = vertices[ia].Coordinates[i];
				coord_b[i] = vertices[ib].Coordinates[i];
				coord_c[i] = vertices[ic].Coordinates[i];
				coord_d[i] = vertices[id].Coordinates[i];
				coord_e[i] = vertices[ie].Coordinates[i];
			}
			ra = vertices[ia].Radius;
			rb = vertices[ib].Radius;
			rc = vertices[ic].Radius;
			rd = vertices[id].Radius;
			re = vertices[ie].Radius;
			sos.sos_minor5_gmp(coord_a, ra, coord_b, rb, coord_c, rc,
			coord_d, rd, coord_e, re, &val);
			det_pijkl = val * nswap * iorient;
		}
		*redundant = (det_pijkl < 0) ? 1 : 0;

	} else if (ninf == 1) { /*
				We know that one of the 4 vertices a,b,c or d is infinite
				To find which one it is, we use a map between
				(inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
				is 1 if i is infinite, 0 otherwise, and X = 0,1,2,3
				if a,b,c or d are infinite, respectively.
				A good mapping function is:
				X = 2 - inf(a) - inf(a) -inf(b) + inf(d)
				*/

		idx = (2 - infPoint[0] - infPoint[0] - infPoint[1] + infPoint[3]);
		l = list[idx];

		i = list[order1[idx][0]];
		j = list[order1[idx][1]];
		k = list[order1[idx][2]];

		ic1 = inf4_1[l];
		sign = sign4_1[l];
		/*
		let us look at the four determinant we need to compute:
		det_pijk	: unchanged
		det_pjil	: 1 infinite point (l), becomes det3_pji where
		det3_pij =	    |p(ic1) p(ic2) 1|
						  |i(ic1) i(ic2) 1|
						  |j(ic1) j(ic2) 1|
		and ic1 and ic2 depends on which infinite ( ic2 is always 3)
		point is considered
		det_pkjl	: 1 infinite point (l), becomes det3_pkj
		det_pikl	: 1 infinite point (l), becomes det3_pik
		 */

		//get Coordinates

		for (int _i = 0; _i < 3; _i++) {
			i_p[_i] = vertices[i].Coordinates[_i] - vertices[p].Coordinates[_i];
			j_p[_i] = vertices[j].Coordinates[_i] - vertices[p].Coordinates[_i];
			k_p[_i] = vertices[k].Coordinates[_i] - vertices[p].Coordinates[_i];
		}

		detij[0] = i_p[0] * j_p[2] - i_p[2] * j_p[0];
		detij[1] = i_p[1] * j_p[2] - i_p[2] * j_p[1];
		detij[2] = i_p[0] * j_p[1] - i_p[1] * j_p[0];

		*is_in = 0; //false

		det_pijk = -k_p[0] * detij[1] + k_p[1] * detij[0] - k_p[2] * detij[2];
		det_pijk = det_pijk * iorient;
		test_pijk = (std::abs(det_pijk) > eps);
		if (test_pijk && det_pijk > 0) {
			*ifail = idx;
			return;
		}

		det_pjil = -detij[ic1] * sign * iorient;
		test_pjil = (std::abs(det_pjil) > eps);
		if (test_pjil && det_pjil > 0) {
			*ifail = order1[idx][2];
			return;
		}

		det_pkjl = k_p[ic1] * j_p[2] - k_p[2] * j_p[ic1];
		det_pkjl = sign * det_pkjl * iorient;
		test_pkjl = (std::abs(det_pkjl) > eps);
		if (test_pkjl && det_pkjl > 0) {
			*ifail = order1[idx][0];
			return;
		}

		det_pikl = i_p[ic1] * k_p[2] - i_p[2] * k_p[ic1];
		det_pikl = sign * det_pikl * iorient;
		test_pikl = (std::abs(det_pikl) > eps);
		if (test_pikl && det_pikl > 0) {
			*ifail = order1[idx][1];
			return;
		}
		/*
		At this stage, either all four determinants are positive,
		or one of the determinant is not precise enough, and
		we need to switch the MP
		 */

		if (!test_pijk) {
			mysort.valsort4(p, i, j, k, &ia, &ib, &ic, &id, &nswap);
			for(int i = 0; i < 3; i++) {
				coord_a[i] = vertices[ia].Coordinates[i];
				coord_b[i] = vertices[ib].Coordinates[i];
				coord_c[i] = vertices[ic].Coordinates[i];
				coord_d[i] = vertices[id].Coordinates[i];
			}
			sos.sos_minor4_gmp(coord_a, coord_b, coord_c, coord_d, &val);
			val = val * nswap * iorient;
			if (val == 1) {
				*ifail = idx;
				return;
			}
		}

		if (!test_pjil) {
			mysort.valsort3(p, j, i, &ia, &ib, &ic, &nswap);
			int temp = 2;
			xa = vertices[ia].Coordinates[ic1];
			ya = vertices[ia].Coordinates[temp];
			xb = vertices[ib].Coordinates[ic1];
			yb = vertices[ib].Coordinates[temp];
			xc = vertices[ic].Coordinates[ic1];
			yc = vertices[ic].Coordinates[temp];
			sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
			val = val * sign * nswap * iorient;
			if (val == 1) {
				*ifail = order1[idx][2];
				return;
			}
		}

		if (!test_pkjl) {
			mysort.valsort3(p, k, j, &ia, &ib, &ic, &nswap);
			int temp = 2;
			xa = vertices[ia].Coordinates[ic1];
			ya = vertices[ia].Coordinates[temp];
			xb = vertices[ib].Coordinates[ic1];
			yb = vertices[ib].Coordinates[temp];
			xc = vertices[ic].Coordinates[ic1];
			yc = vertices[ic].Coordinates[temp];
			sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
			val = val * sign * nswap * iorient;
			if (val == 1) {
				*ifail = order1[idx][0];
				return;
			}
		}

		if (!test_pikl) {
			mysort.valsort3(p, i, k, &ia, &ib, &ic, &nswap);
			int temp = 2;
			xa = vertices[ia].Coordinates[ic1];
			ya = vertices[ia].Coordinates[temp];
			xb = vertices[ib].Coordinates[ic1];
			yb = vertices[ib].Coordinates[temp];
			xc = vertices[ic].Coordinates[ic1];
			yc = vertices[ic].Coordinates[temp];
			sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
			val = val * sign * nswap * iorient;
			if (val == 1) {
				*ifail = order1[idx][1];
				return;
			}
		}

		//If we have gone so far, p is inside the tetrahedron
		*is_in = 1; //true
		/*
		Now we check if p is redundant
		since det_pijkl = det_pijk > 1
		p cannot be redundant !
		 */
		*redundant = 0; //false

	} else if (ninf == 2) { /*

				We know that two of the 4 vertices a,b,c or d are infinite
				To find which one it is, we use a map between
				(inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
				is 1 if i is infinite, 0 otherwise, and X = 1,2,3,4,5,6
				if (a,b), (a,c), (a,d), (b,c), (b,d), or (c,d) are
				infinite, respectively
				A good mapping function is:
				X = 2 - inf(a) - inf(a) +inf(c) + inf(d) + inf(d)
				*/

		idx = (2 - infPoint[0] - infPoint[0] + infPoint[2] + 2 * infPoint[3]);

		//The two infinite points :
		k = list[order3[idx][0]];
		l = list[order3[idx][1]];

		//The two finite points
		i = list[order2[idx][0]];
		j = list[order2[idx][1]];

		ic1_k = inf4_1[k];
		ic1_l = inf4_1[l];
		sign_k = sign4_1[k];
		sign_l = sign4_1[l];
		ic1 = inf4_2[l][k];
		sign = sign4_2[l][k];

		//Get coordinates
		for (int _i = 0; _i < 3; _i++) {
			i_p[_i] = vertices[i].Coordinates[_i] - vertices[p].Coordinates[_i];
			j_p[_i] = vertices[j].Coordinates[_i] - vertices[p].Coordinates[_i];
		}

		//Perform test; first set is_in .false.
		*is_in = 0; //false

		//det_pijk is now det3_pij with k as infinite point
		det_pijk = i_p[ic1_k] * j_p[2] - i_p[2] * j_p[ic1_k];
		det_pijk = det_pijk * sign_k * iorient;
		test_pijk = (std::abs(det_pijk) > eps);
		if (test_pijk && det_pijk > 0) {
			*ifail = order3[idx][1];
			return;
		}

		//det_pjil is now det3_pji with l as infinite point
		det_pjil = i_p[2] * j_p[ic1_l] - i_p[ic1_l] * j_p[2];
		det_pjil = det_pjil * sign_l * iorient;
		test_pjil = (std::abs(det_pjil) > eps);
		if (test_pjil && det_pjil > 0) {
			*ifail = order3[idx][0];
			return;
		}

		//det_pkjl is now -det2_pj (k,l infinite)
		det_pkjl = j_p[ic1] * sign * iorient;
		test_pkjl = (std::abs(det_pkjl) > eps);
		if (test_pkjl && det_pkjl > 0) {
			*ifail = order2[idx][0];
			return;
		}

		//det_pikl is now det2_pi (k,l infinite)
		det_pikl = -i_p[ic1] * sign * iorient;
		test_pikl = (std::abs(det_pikl) > eps);
		if (test_pikl && det_pikl > 0) {
			*ifail = order2[idx][1];
			return;
		}
		/*
		At this stage, either all four determinants are positive,
		or one of the determinant is not precise enough, and
		we need to switch the MP
		 */
		if (!test_pijk) {
			mysort.valsort3(p, i, j, &ia, &ib, &ic, &nswap);
			int temp = 2;
			xa = vertices[ia].Coordinates[ic1_k];
			ya = vertices[ia].Coordinates[temp];
			xb = vertices[ib].Coordinates[ic1_k];
			yb = vertices[ib].Coordinates[temp];
			xc = vertices[ic].Coordinates[ic1_k];
			yc = vertices[ic].Coordinates[temp];
			sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
			val = val * sign_k * nswap * iorient;
			if (val == 1) {
				*ifail = order3[idx][1];
				return;
			}
		}

		if (!test_pjil) {
			mysort.valsort3(p, j, i, &ia, &ib, &ic, &nswap);
			int temp = 2;
			xa = vertices[ia].Coordinates[ic1_l];
			ya = vertices[ia].Coordinates[temp];
			xb = vertices[ib].Coordinates[ic1_l];
			yb = vertices[ib].Coordinates[temp];
			xc = vertices[ic].Coordinates[ic1_l];
			yc = vertices[ic].Coordinates[temp];
			sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
			val = val * sign_l * nswap * iorient;
			if (val == 1) {
				*ifail = order3[idx][0];
				return;
			}
		}

		if (!test_pkjl) {
			mysort.valsort2(p, j, &ia, &ib, &nswap);
			xa = vertices[ia].Coordinates[ic1];
			xb = vertices[ib].Coordinates[ic1];
			sos.sos_minor2_gmp(xa, xb, &val);
			val = -val * sign * nswap * iorient;
			if (val == 1) {
				*ifail = order2[idx][0];
				return;
			}
		}

		if (!test_pikl) {
			mysort.valsort2(p, i, &ia, &ib, &nswap);
			xa = vertices[ia].Coordinates[ic1];
			xb = vertices[ib].Coordinates[ic1];
			sos.sos_minor2_gmp(xa, xb, &val);
			val = val * sign * nswap * iorient;
			if (val == 1) {
				*ifail = order2[idx][1];
				return;
			}
		}

		//Again, if we have gone so far, p is inside the tetrahedron
		*is_in = 1; //true

		//Now we check if p is redundant
		//det_pijkl becomes det3_pij
		ic5 = inf5_2[l][k];
		sign5 = sign5_2[l][k];
		det_pijkl = i_p[ic5] * j_p[2] - i_p[2] * j_p[ic5];
		if (fabs(det_pijkl) < eps) {
			mysort.valsort3(p, i, j, &ia, &ib, &ic, &nswap);
			int temp = 2;
			xa = vertices[ia].Coordinates[ic5];
			ya = vertices[ia].Coordinates[temp];
			xb = vertices[ib].Coordinates[ic5];
			yb = vertices[ib].Coordinates[temp];
			xc = vertices[ic].Coordinates[ic5];
			yc = vertices[ic].Coordinates[temp];
			sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
			det_pijkl = val*nswap;
		}
		det_pijkl = det_pijkl * sign5 * iorient;

		*redundant = (det_pijkl < 0) ? 1 : 0;

	} else if (ninf == 3) { /*
				We know that three of the 4 vertices a,b,c or d are infinite
				To find which one is finite, we use a map between
				(inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
				is 1 if i is infinite, 0 otherwise, and X = 0, 1, 2, 3
				if a,b,c or d are finite, respectively.
				A good mapping function is:
				X = inf(a) + inf(a) +inf(b) - inf(d)
				*/

		idx = 2 * infPoint[0] + infPoint[1] - infPoint[3];

		i = list[idx];
		j = list[order1[idx][0]];
		k = list[order1[idx][1]];
		l = list[order1[idx][2]];

		//Index of the "missing" infinite point (i.e. the fourth infinite point)
		missinf_sign(j, k, l, &ie, &iswap);

		//Get coordinates
		for (int _i = 0; _i < 3; _i++) {
			i_p[_i] = vertices[i].Coordinates[_i] - vertices[p].Coordinates[_i];
		}

		//Perform test; first set is_in to .false.
		*is_in = 0; //false

		//det_pijk is now - det2_pi (missing j,k)
		det_pijk = i_p[inf4_2[k][j]] * iorient * sign4_2[k][j];
		test_pijk = (std::abs(det_pijk) > eps);
		if (test_pijk && det_pijk > 0) {
			*ifail = order1[idx][2];
			return;
		}

		//det_pjil is now det2_pi (missing j,l)
		det_pjil = -i_p[inf4_2[l][j]] * iorient * sign4_2[l][j];
		test_pjil = (std::abs(det_pjil) > eps);
		if (test_pjil && det_pjil > 0) {
			*ifail = order1[idx][1];
			return;
		}

		//det_pkjl is now det1_p
		det_pkjl = iorient * iswap * sign4_3[ie];
		if (det_pkjl > 0) {
			*ifail = idx;
			return;
		}

		//det_ikl is now - det2_pi (missing k,l)
		det_pikl = i_p[inf4_2[l][k]] * iorient * sign4_2[l][k];
		test_pikl = (std::abs(det_pikl) > eps);
		if (test_pikl && det_pikl > 0) {
			*ifail = order1[idx][0];
			return;
		}
		/*
		At this stage, either all four determinants are positive,
		or one of the determinant is not precise enough, and
		we need to switch the MP
		 */
		if (!test_pijk) {
			mysort.valsort2(p, i, &ia, &ib, &nswap);
			xa = vertices[ia].Coordinates[inf4_2[k][j]];
			xb = vertices[ib].Coordinates[inf4_2[k][j]];
			sos.sos_minor2_gmp(xa, xb, &val);
			val = -val * sign4_2[k][j] * iorient * nswap;
			if (val == 1) {
				*ifail = order1[idx][2];
				return;
			}
		}

		if (!test_pjil) {
			mysort.valsort2(p, i, &ia, &ib, &nswap);
			xa = vertices[ia].Coordinates[inf4_2[l][j]];
			xb = vertices[ib].Coordinates[inf4_2[l][j]];
			sos.sos_minor2_gmp(xa, xb, &val);
			val = val * sign4_2[l][j] * iorient * nswap;
			if (val == 1) {
				*ifail = order1[idx][1];
				return;
			}
		}

		if (!test_pikl) {
			mysort.valsort2(p, i, &ia, &ib, &nswap);
			xa = vertices[ia].Coordinates[inf4_2[l][k]];
			xb = vertices[ib].Coordinates[inf4_2[l][k]];
			sos.sos_minor2_gmp(xa, xb, &val);
			val = -val * sign4_2[l][k] * iorient * nswap;
			if (val == 1) {
				*ifail = order1[idx][0];
				return;
			}
		}

		*is_in = 1; //true

		//Now check for redundancy
		//det_pijkl becomes -det2_pi

		ic1 = inf5_3[ie];
		sign5 = sign5_3[ie];
		det_pijkl = -i_p[ic1];
		if (fabs(det_pijkl) < eps) {
			mysort.valsort2(p, i, &ia, &ib, &nswap);
			xa = vertices[ia].Coordinates[ic1];
			xb = vertices[ib].Coordinates[ic1];
			sos.sos_minor2_gmp(xa, xb, &val);
			det_pijkl = val * nswap;
		}
		det_pijkl = -iorient * det_pijkl * sign5 * iswap;
		*redundant = (det_pijkl < 0) ? 1 : 0;

	} else {
		//In the case all four points ia,ib,ic, and id are infinite,
		//then is_in = .true. and redundant = .false.
		*is_in = 1; //true
		*redundant = 0; //false
	}
  }

/* ====================================================================
	Check for local regularity
 ==================================================================== */

  void DELCX::regular_convex(std::vector<Vertex>& vertices, int a, int b, int c, 
	int p, int o, int itest_abcp, bool *regular, bool *convex, bool *test_abpo, 
	bool *test_bcpo, bool *test_capo) 
  {
	int i, j, k, l, m;
	int ia = 0, ib = 0, ic = 0, id = 0, ie = 0;
	int ninf, infp, info, iswap = 0, iswap2 = 0, idx, val = 0;
	int icol1, sign1, icol2, sign2, icol4 = 0, sign4, icol5, sign5;
	int list[3], infPoint[3];

	double det_abpo, det_bcpo, det_capo, det_abcpo, det_abpc;
	double xa, ya, xb, yb, xc, yc;
	double ra, rb, rc, rd, re;
	double a_p[4], b_p[4], c_p[4], o_p[4];
	double i_p[3], j_p[3];
	double Mbo[3], Mca[3], Mjo[3], Mio[3];
	double coord1[3], coord2[3], coord3[3], coord4[3], coord5[3];

	bool testc[3];

	*regular = true;
	*convex  = true;
	*test_abpo = false;
	*test_bcpo = false;
	*test_capo = false;

	/*
	To test if the union of the two tetrahedron is convex, check the position of o w.r.t
	three faces (a,b,p), (b,c,p) and (c,a,p) of (a,b,c,p).evaluate the three determinants:
	det(a,b,p,o)
	det(b,c,p,o)
	det(c,a,p,o)
	If the three determinants are positive, & det(a,b,c,p) is negative,the union is convex
	If the three determinants are negative, & det(a,b,c,p) is positive,the union is convex
	In all other cases, the union is non convex
	The regularity is tested by computing det(a,b,c,p,o)

	first count how many infinite points (except o)
	only a and/or b and/or c can be infinite:
	 */
	infPoint[0] = 0;
	infPoint[1] = 0;
	infPoint[2] = 0;
	if(a < 4) infPoint[0] = 1;
	if(b < 4) infPoint[1] = 1;
	if(c < 4) infPoint[2] = 1;

	ninf = infPoint[0] + infPoint[1] + infPoint[2];

	list[0] = a;
	list[1] = b;
	list[2] = c;

	//general case:no inf points
	if (ninf == 0) {
		//First, simple case:
		//if o is infinite, then det(a,b,c,p,o) = -det(a,b,c,p)
		//and consequently (a,b,c,p,o) is regular:nothing to do!
		if (o < 4) {
			*regular = true;
			return;
		}
		/*
		The three determinants det(a,b,p,o), det(b,c,p,o), and det(c,a,p,o) are "real" 4x4 determinants.
		Subtract the row corresponding to p from the other row, and develop with respect to p.
		The determinants become:
		det(a,b,p,o)= - | ap(1) ap(2) ap(3) |
				| bp(1) bp(2) bp(3) |
				| op(1) op(2) op(3) |

		det(b,c,p,o)=  -| bp(1) bp(2) bp(3) |
				| cp(1) cp(2) cp(3) |
				| op(1) op(2) op(3) |

		det(c,a,p,o)=  -| cp(1) cp(2) cp(3) |
				| ap(1) ap(2) ap(3) |
				| op(1) op(2) op(3) |

		Where ip(j) = i(j) - p(j) for all i in {a,b,c,o} and j in {1,2,3}
				Compute two types of minors:

		Mbo_ij = bp(i)op(j) - bp(j)op(i) and Mca_ij = cp(i)ap(j) - cp(j)op(i)

		Store Mbo_12 in Mbo(3), Mbo_13 in Mbo(2),...
		 */
		//get the coordinates
		for (m = 0; m < 3; m++) {
			a_p[m] = vertices[a].Coordinates[m] - vertices[p].Coordinates[m];
			b_p[m] = vertices[b].Coordinates[m] - vertices[p].Coordinates[m];
			c_p[m] = vertices[c].Coordinates[m] - vertices[p].Coordinates[m];
			o_p[m] = vertices[o].Coordinates[m] - vertices[p].Coordinates[m];
		}
		a_p[3] = vertices[a].Weight - vertices[p].Weight;
		b_p[3] = vertices[b].Weight - vertices[p].Weight;
		c_p[3] = vertices[c].Weight - vertices[p].Weight;
		o_p[3] = vertices[o].Weight - vertices[p].Weight;

		//compute 2x2 determinants Mbo and Mca
		Mbo[0] = b_p[1] * o_p[2] - b_p[2] * o_p[1];
		Mbo[1] = b_p[0] * o_p[2] - b_p[2] * o_p[0];
		Mbo[2] = b_p[0] * o_p[1] - b_p[1] * o_p[0];

		Mca[0] = c_p[1] * a_p[2] - c_p[2] * a_p[1];
		Mca[1] = c_p[0] * a_p[2] - c_p[2] * a_p[0];
		Mca[2] = c_p[0] * a_p[1] - c_p[1] * a_p[0];

		det_abpo = -a_p[0] * Mbo[0] + a_p[1] * Mbo[1] - a_p[2] * Mbo[2];
		det_bcpo = c_p[0] * Mbo[0] - c_p[1] * Mbo[1] + c_p[2] * Mbo[2];
		det_capo = -o_p[0] * Mca[0] + o_p[1] * Mca[1] - o_p[2] * Mca[2];
		det_abpc = -b_p[0] * Mca[0] + b_p[1] * Mca[1] - b_p[2] * Mca[2];
		/*
		To compute:
		det(a,b,c,p,o) =	| a(1) a(2) a(3) a(4) 1 |
					| b(1) b(2) b(3) b(4) 1 |
					| c(1) c(2) c(3) c(4) 1 |
					| p(1) p(2) p(3) p(4) 1 |
					| o(1) o(2) o(3) o(4) 1 |

		First substract row p :
		det(a,b,c,p,o) =      - | ap(1) ap(2) ap(3) ap(4) |
					| bp(1) bp(2) bp(3) bp(4) |
					| cp(1) cp(2) cp(3) cp(4) |
					| op(1) op(2) op(3) op(4) |

		Expand w.r.t last column
		 */
		det_abcpo = -a_p[3] * det_bcpo - b_p[3] * det_capo - c_p[3] * det_abpo + o_p[3] * det_abpc;

		//if (a,b,c,p,o) regular, no need to flip

		if (std::abs(det_abcpo) < eps) {
			mysort.valsort5(a, b, c, p, o, &ia, &ib, &ic, &id, &ie, &iswap);
			for(int i = 0; i < 3; i++) {
				coord1[i] = vertices[ia].Coordinates[i];
				coord2[i] = vertices[ib].Coordinates[i];
				coord3[i] = vertices[ic].Coordinates[i];
				coord4[i] = vertices[id].Coordinates[i];
				coord5[i] = vertices[ie].Coordinates[i];
			}
			ra = vertices[ia].Radius; rb = vertices[ib].Radius;
			rc = vertices[ic].Radius; rd = vertices[id].Radius;
			re = vertices[ie].Radius;
			sos.sos_minor5_gmp(coord1, ra, coord2, rb, coord3, rc, 
			coord4, rd, coord5, re, &val);
			det_abcpo = val * iswap;
		}

		if ((det_abcpo * itest_abcp) < 0) {
			*regular = true;
			return;
		}
		*regular = false;

		//If not regular, we test for convexity
		if (std::abs(det_abpo) < eps) {
			mysort.valsort4(a, b, p, o, &ia, &ib, &ic, &id, &iswap);
			for(int i = 0; i < 3; i++) {
				coord1[i] = vertices[ia].Coordinates[i];
				coord2[i] = vertices[ib].Coordinates[i];
				coord3[i] = vertices[ic].Coordinates[i];
				coord4[i] = vertices[id].Coordinates[i];
			}
			sos.sos_minor4_gmp(coord1, coord2, coord3, coord4, &val);
			det_abpo = val * iswap;
		}
		if (std::abs(det_bcpo) < eps) {
			mysort.valsort4(b, c, p, o, &ia, &ib, &ic, &id, &iswap);
			for(int i = 0; i < 3; i++) {
				coord1[i] = vertices[ia].Coordinates[i];
				coord2[i] = vertices[ib].Coordinates[i];
				coord3[i] = vertices[ic].Coordinates[i];
				coord4[i] = vertices[id].Coordinates[i];
			}
			sos.sos_minor4_gmp(coord1, coord2, coord3, coord4, &val);
			det_bcpo = val * iswap;
		}
		if (std::abs(det_capo) < eps) {
			mysort.valsort4(c, a, p, o, &ia, &ib, &ic, &id, &iswap);
			for(int i = 0; i < 3; i++) {
				coord1[i] = vertices[ia].Coordinates[i];
				coord2[i] = vertices[ib].Coordinates[i];
				coord3[i] = vertices[ic].Coordinates[i];
				coord4[i] = vertices[id].Coordinates[i];
			}
			sos.sos_minor4_gmp(coord1, coord2, coord3, coord4, &val);
			det_capo = val * iswap;
		}

		*test_abpo = (det_abpo > 0);
		*test_bcpo = (det_bcpo > 0);
		*test_capo = (det_capo > 0);

		*convex = false;
		if ((itest_abcp * det_abpo) > 0) return;
		if ((itest_abcp * det_bcpo) > 0) return;
		if ((itest_abcp * det_capo) > 0) return;
		*convex = true;

	} else if (ninf == 1) { 
		/*
		Define X as infinite point, and (i,j) the pair of finite points.
		If X = a, (i,j) = (b,c)
		If X = b, (i,j) = (c,a)
		If X = c, (i,j) = (a,b)
		If we define inf(a) = 1 if a infinite, 0 otherwise,
		then idx_X  = 2 - inf(a) + inf(c)
		*/
		idx = 1 - infPoint[0] + infPoint[2];
		infp = list[idx];
		i = list[order[idx][0]];
		j = list[order[idx][1]];

		//Get the coordinates
		for (m = 0; m < 3; m++) {
			i_p[m] = vertices[i].Coordinates[m] - vertices[p].Coordinates[m];
			j_p[m] = vertices[j].Coordinates[m] - vertices[p].Coordinates[m];
		}

		if (o > 3) {
			//First case:	o is finite
			for (m = 0; m < 3; m++) {
				o_p[m] = vertices[o].Coordinates[m] - vertices[p].Coordinates[m];
			}

			icol1 = inf4_1[infp];
			sign1 = sign4_1[infp];
			/*
			The three 4x4 determinants become:
			-det(i,p,o) [X missing], det(j,p,o) [X missing],det(i,j,p,o)
			And the 5x5 determinant becomes:
			- det(i,j,p,o)
			 */
			Mjo[0] = j_p[0] * o_p[2] - j_p[2] * o_p[0];
			Mjo[1] = j_p[1] * o_p[2] - j_p[2] * o_p[1];
			Mjo[2] = j_p[0] * o_p[1] - j_p[1] * o_p[0];
			/*
			The correspondence between a,b,c and i,j is not essential
			We use corresponce for a infinite; in thetwo other cases
			(b infinite or c infinite),computed determinants are the
			same , but not in the same order
			 */
			det_abpo = i_p[icol1] * o_p[2] - i_p[2] * o_p[icol1];
			if (std::abs(det_abpo) < eps) {
				int temp = 2;
				mysort.valsort3(i, p, o, &ia, &ib, &ic, &iswap);
				xa = vertices[ia].Coordinates[icol1];
				ya = vertices[ia].Coordinates[temp];
				xb = vertices[ib].Coordinates[icol1];
				yb = vertices[ib].Coordinates[temp];
				xc = vertices[ic].Coordinates[icol1];
				yc = vertices[ic].Coordinates[temp];
				sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
				det_abpo = -val * iswap;
			}
			det_abpo = det_abpo * sign1;
			det_capo = -Mjo[icol1];
			if (std::abs(det_capo) < eps) {
				int temp = 2;
				mysort.valsort3(j, p, o, &ia, &ib, &ic, &iswap);
				xa = vertices[ia].Coordinates[icol1];
				ya = vertices[ia].Coordinates[temp];
				xb = vertices[ib].Coordinates[icol1];
				yb = vertices[ib].Coordinates[temp];
				xc = vertices[ic].Coordinates[icol1];
				yc = vertices[ic].Coordinates[temp];
				sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
				det_capo = val * iswap;
			}
			det_capo = det_capo * sign1;
			det_bcpo = -i_p[0] * Mjo[1] + i_p[1] * Mjo[0] - i_p[2] * Mjo[2];
			if (std::abs(det_bcpo) < eps) {
				mysort.valsort4(i, j, p, o, &ia, &ib, &ic, &id, &iswap);
				for(int i = 0; i < 3; i++) {
					coord1[i] = vertices[ia].Coordinates[i];
					coord2[i] = vertices[ib].Coordinates[i];
					coord3[i] = vertices[ic].Coordinates[i];
					coord4[i] = vertices[id].Coordinates[i];
				}
				sos.sos_minor4_gmp(coord1, coord2, coord3, coord4, &val);
				det_bcpo = val * iswap;
			}
			det_abcpo = -det_bcpo;
		} else {
			//Second case: o is infinite
			info = o;
			/*
			The three 4x4 determinants become:
			-det(i,p) [o,X missing]
			det(j,p) [o,X missing]
			det(i,j,p) [o missing]
			And the 5x5 determinant becomes:
			det(i,j,p) [o,X missing]
			 */
			icol1 = inf4_2[infp][info];
			sign1 = sign4_2[infp][info];

			icol2 = inf4_1[info];
			sign2 = sign4_1[info];

			icol5 = inf5_2[infp][info];
			sign5 = sign5_2[infp][info];


			det_abpo = -i_p[icol1] * sign1;
			if (std::abs(det_abpo) < eps) {
				mysort.valsort2(i, p, &ia, &ib, &iswap);
				xa = vertices[ia].Coordinates[icol1];
				xb = vertices[ib].Coordinates[icol1];
				sos.sos_minor2_gmp(xa, xb, &val);
				det_abpo = -val * iswap * sign1;
			}
			det_capo = j_p[icol1] * sign1;
			if (std::abs(det_capo) < eps) {
				mysort.valsort2(j, p, &ia, &ib, &iswap);
				xa = vertices[ia].Coordinates[icol1];
				xb = vertices[ib].Coordinates[icol1];
				sos.sos_minor2_gmp(xa, xb, &val);
				det_capo = val * iswap * sign1;
			}
			det_bcpo = i_p[icol2] * j_p[2] - i_p[2] * j_p[icol2];
			if (std::abs(det_bcpo) < eps) {
				int temp = 2;
				mysort.valsort3(i, j, p, &ia, &ib, &ic, &iswap);
				xa = vertices[ia].Coordinates[icol2];
				ya = vertices[ia].Coordinates[temp];
				xb = vertices[ib].Coordinates[icol2];
				yb = vertices[ib].Coordinates[temp];
				xc = vertices[ic].Coordinates[icol2];
				yc = vertices[ic].Coordinates[temp];
				sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
				det_bcpo = val * iswap;
			}
			det_bcpo = det_bcpo * sign2;
			det_abcpo = i_p[icol5] * j_p[2] - i_p[2] * j_p[icol5];
			if (std::abs(det_abcpo) < eps) {
				int temp = 2;
				mysort.valsort3(i, j, p, &ia, &ib, &ic, &iswap);
				xa = vertices[ia].Coordinates[icol5];
				ya = vertices[ia].Coordinates[temp];
				xb = vertices[ib].Coordinates[icol5];
				yb = vertices[ib].Coordinates[temp];
				xc = vertices[ic].Coordinates[icol5];
				yc = vertices[ic].Coordinates[temp];
				sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
				det_abcpo = val * iswap;
			}
			det_abcpo = det_abcpo * sign5;
		}

		//Test if (a,b,c,p,o) regular, in which case there is no need to flip
		if ((det_abcpo * itest_abcp) < 0) {
			*regular = true;
			return;
		}
		*regular = false;

		//If not regular, we test for convexity

		testc[0] = (det_abpo > 0);
		testc[1] = (det_bcpo > 0);
		testc[2] = (det_capo > 0);
		*test_abpo = testc[ord_rc[idx][0]];
		*test_bcpo = testc[ord_rc[idx][1]];
		*test_capo = testc[ord_rc[idx][2]];

		*convex = false;
		if ((itest_abcp * det_abpo) > 0) return;
		if ((itest_abcp * det_bcpo) > 0) return;
		if ((itest_abcp * det_capo) > 0) return;
		*convex = true;

	} else if (ninf == 2) { /*
				define(k,l) as the two infinite points, and i be finite
				If i = a, (k,l) = (b,c)
				If i = b, (k,l) = (c,a)
				If i = c, (k,l) = (a,b)

				Again: idx = 2 + inf(a) - inf(c)
				*/
		idx = 1 + infPoint[0] - infPoint[2];
		i = list[idx];
		k = list[order[idx][0]];
		l = list[order[idx][1]];

		//Get the coordinates

		for (m = 0; m < 3; m++) {
			i_p[m] = vertices[i].Coordinates[m] - vertices[p].Coordinates[m];
		}

		if (o > 3) {
			//First case: o is finite
			/*
			The three 4x4 determinants become:
			det(i,p,o) [k missing]
			-det(i,p,o) [l missing]
			S*det(p,o) [k,l missing, with S =1 if k<l, -1 otherwise]
			The 5x5 determinants become:
			S*det(i,p,o) [k,l missing, with S=1 if k<l, -1 otherwise]
			 */
			for (m = 0; m < 3; m++) {
				o_p[m] = vertices[o].Coordinates[m] - vertices[p].Coordinates[m];
			}
			icol1 = inf4_1[k];
			sign1 = sign4_1[k];
			icol2 = inf4_1[l];
			sign2 = sign4_1[l];
			icol4 = inf4_2[l][k];
			sign4 = sign4_2[l][k];
			icol5 = inf5_2[l][k];
			sign5 = sign5_2[l][k];

			Mio[0] = i_p[0] * o_p[2] - i_p[2] * o_p[0];
			Mio[1] = i_p[1] * o_p[2] - i_p[2] * o_p[1];
			Mio[2] = i_p[0] * o_p[1] - i_p[1] * o_p[0];
			/*
			The correspondence between a,b,c and i,j,k is not essential
			use the correspondence for a finite; in the two other cases
			(b finite or c finite),have computed the same determinants,
			but not in the same order
			 */
			det_abpo = -Mio[icol1] * sign1;
			if (std::abs(det_abpo) < eps) {
				int temp = 2;
				mysort.valsort3(i, p, o, &ia, &ib, &ic, &iswap);
				xa = vertices[ia].Coordinates[icol1];
				ya = vertices[ia].Coordinates[temp];
				xb = vertices[ib].Coordinates[icol1];
				yb = vertices[ib].Coordinates[temp];
				xc = vertices[ic].Coordinates[icol1];
				yc = vertices[ic].Coordinates[temp];
				sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
				det_abpo = val * iswap * sign1;
			}
			det_capo = Mio[icol2] * sign2;
			if (std::abs(det_capo) < eps) {
				int temp = 2;
				mysort.valsort3(i, p, o, &ia, &ib, &ic, &iswap);
				xa = vertices[ia].Coordinates[icol2];
				ya = vertices[ia].Coordinates[temp];
				xb = vertices[ib].Coordinates[icol2];
				yb = vertices[ib].Coordinates[temp];
				xc = vertices[ic].Coordinates[icol2];
				yc = vertices[ic].Coordinates[temp];
				sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
				det_capo = -val * iswap * sign2;
			}
			det_bcpo = -o_p[icol4] * sign4;
			if (std::abs(det_bcpo) < eps) {
				mysort.valsort2(p, o, &ia, &ib, &iswap);
				xa = vertices[ia].Coordinates[icol4];
				xb = vertices[ib].Coordinates[icol4];
				sos.sos_minor2_gmp(xa, xb, &val);
				det_bcpo = val * sign4 * iswap;
			}
			det_abcpo = -Mio[icol5] * sign5;
			if (std::abs(det_abcpo) < eps) {
				int temp = 2;
				mysort.valsort3(i, p, o, &ia, &ib, &ic, &iswap);
				xa = vertices[ia].Coordinates[icol5];
				ya = vertices[ia].Coordinates[temp];
				xb = vertices[ib].Coordinates[icol5];
				yb = vertices[ib].Coordinates[temp];
				xc = vertices[ic].Coordinates[icol5];
				yc = vertices[ic].Coordinates[temp];
				sos.sos_minor3_gmp(xa, ya, xb, yb, xc, yc, &val);
				det_abcpo = val * iswap * sign5;
			}
		} else {
			//Second case: o is infinite
			info = o;
			/*
			The three 4x4 determinants become:
			det(i,p) [o,k missing]
			-det(i,p) [o,l missing]
			Const [o,k,l missing]
			The 5x5 determinants become:
			Const*det(i,p) [o,k,l missing]
			 */
			icol1 = inf4_2[k][info];
			sign1 = sign4_2[k][info];
			icol2 = inf4_2[l][info];
			sign2 = sign4_2[l][info];

			missinf_sign(info, k, l, &icol4, &iswap);

			det_abpo = i_p[icol1] * sign1;
			if (std::abs(det_abpo) < eps) {
				mysort.valsort2(i, p, &ia, &ib, &iswap2);
				xa = vertices[ia].Coordinates[icol1];
				xb = vertices[ib].Coordinates[icol1];
				sos.sos_minor2_gmp(xa, xb, &val);
				det_abpo = val * iswap2 * sign1;
			}
			det_capo = -i_p[icol2] * sign2;
			if (std::abs(det_capo) < eps) {
				mysort.valsort2(i, p, &ia, &ib, &iswap2);
				xa = vertices[ia].Coordinates[icol2];
				xb = vertices[ib].Coordinates[icol2];
				sos.sos_minor2_gmp(xa, xb, &val);
				det_capo = -val * iswap2 * sign2;
			}
			det_bcpo = sign4_3[icol4] * iswap;
			det_abcpo = sign5_3[icol4] * iswap * i_p[inf5_3[icol4]];
			if (std::abs(det_abcpo) < eps) {
				mysort.valsort2(i, p, &ia, &ib, &iswap2);
				xa = vertices[ia].Coordinates[inf5_3[icol4]];
				xb = vertices[ib].Coordinates[inf5_3[icol4]];
				sos.sos_minor2_gmp(xa, xb, &val);
				det_abcpo = val * iswap2 * iswap * sign5_3[icol4];
			}
		}
		//if (a,b,c,p,o) regular, no need to flip

		if ((det_abcpo * itest_abcp) < 0) {
			*regular = true;
			return;
		}
		*regular = false;

		//If not regular, we test for convexity

		testc[0] = (det_abpo > 0);
		testc[1] = (det_bcpo > 0);
		testc[2] = (det_capo > 0);
		*test_abpo = testc[ord_rc[idx][0]];
		*test_bcpo = testc[ord_rc[idx][1]];
		*test_capo = testc[ord_rc[idx][2]];

		*convex = false;
		if ((itest_abcp * det_abpo) > 0) return;
		if ((itest_abcp * det_bcpo) > 0) return;
		if ((itest_abcp * det_capo) > 0) return;
		*convex = true;

	} else if (ninf == 3) {
		assert(true);
		//this should not happen
	}
  }

/* ======================================================================
   Sign associated with the missing infinite point
 ====================================================================== */

  void DELCX::missinf_sign(int i, int j, int k, int *l, int *sign) 
  {

	int a, b, c, d;

	*l = 6 - i - j - k;

	a = i;
	b = j;
	c = k;

	*sign = 1;

	if (a > b) {
	    d = a;
	    a = b;
	    b = d;
	    *sign = -*sign;
	}

	if (a > c) {
	    d = a;
	    a = c;
	    c = d;
	    *sign = -*sign;
	}

	if (b > c) {
	    *sign = -*sign;
	}
  }

/* ======================================================================================
	Flip_2_3

	This procedure implements a 2->3 flip in 3D for regular triangulation

	a 2->3 flip is a transformation in which two tetrahedrons are
	flipped into three tetrahedra. The two tetrahedra (abcp) and
	(abco) shares a triangle (abc) which is in the link_facet of the
	current point p added to the triangulation. 
	This flip is only possible if the union of the two tetrahedron is
	convex, and if their shared triangle is not locally regular. 
	We assume here that these tests have been performed and are
	true. 
	Once the flip has been performed, three new tetrahedra are added
	and three new "link facet" are added to the link
	facet queue

	Input:
		- itetra:	index of the tetrahedra (a,b,c,p) considered
		- jtetra:	index of the tetrahedra (a,b,c,o) considered
		- vertices:	the five vertices a,b,c,o,p
		- facei		indices of the vertices a,b,c in (a,b,c,p)
		- facej		indices of the vertices a,b,c in (a,b,c,o)
		- test_abpo:	orientation of the four points a,b,p,o
		- test_bcpo:	orientation of the four points b,c,p,o
		- test_capo:	orientation of the four points c,a,p,o

	Output:
		- nlink_facet:	3 new link facets are added
		- link_facet:	the three faces of the initial tetrahedron
				(a,b,c,o) containing the vertex o are added
				as link facets
		- link_index:   A link_facet is a triangle defined from its
				two neighbouring tetrahedra. I store the position
				of the vertex opposite to the triangle in each
				tetrehedron in the array link_index
		- ierr:		1 if flip was not possible
 ======================================================================================*/

  void DELCX::flip_2_3(std::vector<Tetrahedron>& tetra, int itetra, int jtetra, 
	int *vertices, int *facei, int *facej, bool test_abpo, bool test_bcpo, 
	bool test_capo, int *ierr, int *tetra_last)
  {

	int k,p,o;
	int it, jt, idx, jdx;
	int newtetra;

	std::bitset<8> ikeep, jkeep;

	int face[3];
	int itetra_touch[3], jtetra_touch[3];

	char jtetra_idx[3], itetra_idx[3];

	int idx_list[3][2]={{0,0},{0,1},{1,1}};


	int tests[3], position[3];


	*ierr = 0;

/* ======================================================================================
	If itetra or jtetra are inactive, cannot flip
 ======================================================================================*/

	if(tetra[itetra].info[1] == 0 || tetra[jtetra].info[1] == 0) {
		*ierr = 1;
		return;
	}

/* ======================================================================================

	Define
	- itetra_touch: the three tetrahedra that touches itetra on the
			faces opposite to the 3 vertices a,b,c
	- itetra_idx:	for the three tetrahedra defined by itetra_touch,
			index of the vertex opposite to the face
			common with itetra
	- jtetra_touch: the three tetrahedra that touches jtetra on the
			faces opposite to the 3 vertices a,b,c
	- jtetra_idx:	for the three tetrahedra defined by jtetra_touch,
			index of the vertex opposite to the face
			common with jtetra
 ======================================================================================*/

	for(int i = 0; i < 3; i++)
	{
		itetra_touch[i] = tetra[itetra].Neighbours[facei[i]];
		jtetra_touch[i] = tetra[jtetra].Neighbours[facej[i]];
		itetra_idx[i] = tetra[itetra].nindex[facei[i]];
		jtetra_idx[i] = tetra[jtetra].nindex[facej[i]];
	}

/* ======================================================================================
	First three vertices define face that is removed
 ======================================================================================*/

	face[0] = vertices[0];
	face[1]	= vertices[1];
	face[2] = vertices[2];

	p = vertices[3];
	o = vertices[4];

/* ======================================================================================
	The three new tetrahedra are stored in : 
	- any free space in the tetrahedron list,
	- at the end of the list of known tetrahedra if needed
 ======================================================================================*/

	k = 0;

	while (!free.empty() && k < 3) {
		position[k] = free.top();
		free.pop();
		k++;
	}

	for (int l = k; l < 3; l++) {
		Tetrahedron t;
		t.init();
		position[l] = tetra.size();
		tetra.push_back(t);
	}
        *tetra_last = position[2];

/* ======================================================================================
	Set itetra and jtetra to 0, and add them to kill list
 ======================================================================================*/

	ikeep = tetra[itetra].info;
	jkeep = tetra[jtetra].info;

	tetra[itetra].info[1] = 0;
	tetra[jtetra].info[1] = 0;

	kill.push_back(itetra);
	kill.push_back(jtetra);

/* ======================================================================================
	define the three new tetrahedra: (bcop), (acop) and (abop)
	as well as their neighbours

	tetrahedron bcop : neighbours are acop, abop, neighbour of (abcp)
			   on face bcp, and neighbour of (abco) on face bco
	tetrahedron acop : neighbours are bcop, abop, neighbour of (abcp)
			   on face acp, and neighbour of (abco) on face aco
	tetrahedron abop : neighbours are bcop, acop, neighbour of (abcp)
			   on face abp, and neighbour of (abco) on face abo
 ======================================================================================*/

	tests[0] = 1;
	if(test_bcpo) tests[0] = -1;
	tests[1] = -1;
	if(test_capo) tests[1] = 1;
	tests[2] = 1;
	if(test_abpo) tests[2] = -1;

	for(int i = 0; i < 3; i++)
	{
		newtetra = position[i];

		k = 0;
		for(int j = 0; j < 3; j++)
		{
			if(j==i) continue;
			tetra[newtetra].Vertices[k] = face[j];
			tetra[newtetra].Neighbours[k] = position[j];
			tetra[newtetra].nindex[k]=idx_list[i][k];
			k++;
		}

		tetra[newtetra].Vertices[2] = o;
		it = itetra_touch[i];
		idx = itetra_idx[i];
		tetra[newtetra].Neighbours[2] = it;
		tetra[newtetra].nindex[2] = idx;
		tetra[newtetra].info[5] = ikeep[2+facei[i]];
		if(it !=-1 && idx !=-1) {
			tetra[it].Neighbours[idx]=newtetra;
			tetra[it].nindex[idx] = 2;
		}

		tetra[newtetra].Vertices[3] = p;
		jt = jtetra_touch[i];
		jdx = jtetra_idx[i];
		tetra[newtetra].Neighbours[3] = jt;
		tetra[newtetra].nindex[3] = jdx;
		tetra[newtetra].info[6] = jkeep[2+facej[i]];
		if(jt !=-1 && jdx !=-1) {
			tetra[jt].Neighbours[jdx]=newtetra;
			tetra[jt].nindex[jdx] = 3;
		}

		tetra[newtetra].info[1] = 1;
 
		if(tests[i]==1) {
			tetra[newtetra].info[0] = 1;
		} else {
			tetra[newtetra].info[0] = 0;
		}
	}

/* ======================================================================================

	Now add all three faces of jtetra containing o in the link_facet queue.
	Each link_facet (a triangle) is implicitly defined as the
	intersection of two tetrahedra

	link_facet:	bco	tetrahedra:	bcop and neighbour of (abco)
						on bco
	link_facet:	aco	tetrahedra:	acop and neighbour of (abco)
						on aco
	link_facet:	abo	tetrahedra:	abop and neighbour of (abco)
						on abo
 ======================================================================================*/

	for(int i = 0; i < 3; i++) {
		newtetra = position[i];
		link_facet.push(std::make_pair(newtetra, tetra[newtetra].Neighbours[3]));
		idx = tetra[newtetra].nindex[3];
		link_index.push(std::make_pair(3,idx));
	}

  }

/* ======================================================================================
	This procedure implements a 4->1 flip in 3D for regular triangulation

	a 4->1 flip is a transformation in which a tetrahedron and a single
	vertex included in the tetrahedron are transformed to 4 tetrahedra,
	defined from the 4 four faces of the initial tetrahedron, connected
	to the new point. Each of the faces are then called "link facet",
	and stored on a queue

	Input:
		- ipoint:	index of the point p to be included
		- itetra:	index of the tetrahedra (a,b,c,d) considered

	Output:
		- nlink_facet:	4
		- link_facet:	Add the four faces of the initial tetrahedron
		- link_index:	A link_facet is a triangle defined from its
				two neighboring tetrahedra. I store the position
				of the vertex opposite to the triangle in each
				tetrehedron in the array link_index
 ======================================================================================*/

  void DELCX::flip_1_4(std::vector<Tetrahedron>& tetra, int ipoint, int itetra, 
		int *tetra_last)
  {
	int k, newtetra;
	int jtetra;
	int fact,idx;

	std::bitset<8> ikeep;
	char ival[4];

	int vertex[4], neighbour[4];
	int position[4];
	int idx_list[4][3]={{0,0,0},{0,1,1},{1,1,2},{2,2,2}};

/* ======================================================================================
	Store information about "old" tetrahedron
 ======================================================================================*/

	ikeep = tetra[itetra].info;

	for(int i = 0; i < 4; i++) {
		vertex[i] = tetra[itetra].Vertices[i];
		neighbour[i] = tetra[itetra].Neighbours[i];
		ival[i] = tetra[itetra].nindex[i];
	}
 
	fact = -1;
	if(tetra[itetra].info[0] == 1) fact = 1;

/* ======================================================================================
	The four new tetrahedra are going to be stored in : 
	any free space in the tetrahedron list, and at the end of the list of known tetrahedra
 ======================================================================================*/

	k = 0;
	while (!free.empty() && k < 4) {
		position[k] = free.top();
		free.pop();
		k++;
	}

	for (int l = k; l < 4; l++) {
		Tetrahedron t;
		t.init();
		position[l] = tetra.size();
		tetra.push_back(t);
	}
	*tetra_last = position[3];

/* ======================================================================================
	itetra is set to 0, and added to the "kill" list
 ======================================================================================*/

	tetra[itetra].info[1] = 0;
	kill.push_back(itetra);

/* ======================================================================================
	The tetrahedron is defined as (ijkl); four new tetrahedra are
	created:	jklp, iklp, ijlp, and ijkp, where p is the new
	point to be included

	For each new tetrahedron, define all four neighbours:
	For each neighbour, I store the index of the vertex opposite to 
	the common face in array tetra_nindex

	tetrahedron jklp : neighbours are iklp, ijlp, ijkp and neighbour
			   of (ijkl) on face jkl
	tetrahedron iklp : neighbours are jklp, ijlp, ijkp and neighbour
			   of (ijkl) on face ikl
	tetrahedron ijlp : neighbours are jklp, iklp, ijkp and neighbour
			   of (ijkl) on face ijl
	tetrahedron ijkp : neighbours are jklp, iklp, ijlp and neighbour
			   of (ijkl) on face ijk
 ======================================================================================*/

	for(int i = 0; i < 4; i++)
	{
		newtetra = position[i];

		k = 0;
		for(int j = 0; j < 4; j++) {
			if(j==i) continue;
			tetra[newtetra].Vertices[k] = vertex[j];
			tetra[newtetra].Neighbours[k] = position[j];
			tetra[newtetra].nindex[k] = idx_list[i][k];
			k++;
		}

		jtetra = neighbour[i];
		idx = ival[i];
		tetra[newtetra].Vertices[3] = ipoint;
		tetra[newtetra].Neighbours[3] = jtetra;
		tetra[newtetra].nindex[3] = ival[i];
		tetra[newtetra].info[2+i] = ikeep[2+i];

		if(jtetra !=-1 && idx != -1) {
			tetra[jtetra].Neighbours[idx] = newtetra;
			tetra[jtetra].nindex[idx] = 3;
		}

		tetra[newtetra].info[1] = 1;

		fact = -fact;
		tetra[newtetra].info[0] = 0;
		if(fact==1) tetra[newtetra].info[0] = 1;

	}

/* ======================================================================================
	Now add all fours faces of itetra in the link_facet queue.
	Each link_facet (a triangle) is implicitly defined as the
	intersection of two tetrahedra

	link_facet:	jkl	tetrahedra:	jklp and neighbour of (ijkl)
						on jkl
	link_facet:	ikl	tetrahedra:	iklp and neighbour of (ijkl)
						on ikl
	link_facet:	ijl	tetrahedra:	ijlp and neighbour of (ijkl)
						on ijl
	link_facet:	ijk	tetrahedra:	ijkp and neighbour of (ijkl)
						on ijk
 ======================================================================================*/

	for(int i = 0; i < 4; i++) {
		newtetra = position[i];
		link_facet.push(std::make_pair(newtetra, tetra[newtetra].Neighbours[3]));
		idx = tetra[newtetra].nindex[3];
		link_index.push(std::make_pair(3,idx));
	}

  }

/* ======================================================================================
	flip_3_2 implements a 3->2 flip in 3D for regular triangulation

	a 3->2 flip is a transformation in which three tetrahedrons are
	flipped into two tetrahedra. The two tetrahedra (abpo), (abcp) and
	(abco) shares an edge (ab) which is in the link_facet of the
	current point p added to the triangulation. 
	This flip is only possible if the edge ab is reflex, with degree 3
	We assume here that these tests have been performed and are
	true. 
	Once the flip has been performed, two new tetrahedra are added
	and two new "link facet" are added to the link facet queue

	Input:
		- itetra:	index of the tetrahedra (a,b,c,p) considered
		- jtetra:	index of the tetrahedra (a,b,c,o) considered
		- ktetra:	index of the tetrahedra (a,b,o,p) considered
		- vertices:	the five vertices a,b,c,p,o
		- edgei		indices of a,b in (a,b,c,p)
		- edgej		indices of a,b in (a,b,c,o)
		- edgek		indices of a,b in (a,b,o,p)
		- test_bcpo   : orientation of the four points b,c,p,o
		- test_acpo   : orientation of the four points a,c,p,o

	Output:
		- nlink_facet:	2 new link facets are added
		- link_facet:	the two faces of the initial tetrahedron
				(a,b,o,p) containing the edge op are added
				as link facets
		- link_index:   A link_facet is a triangle defined from its
				two neighbouring tetrahedra. I store the position
				of the vertex opposite to the triangle in each
				tetrehedron in the array link_index
		- ierr:		1 if flip was not possible
 ======================================================================================*/

  void DELCX::flip_3_2(std::vector<Tetrahedron>& tetra, int itetra, int jtetra, 
		int ktetra, int *vertices, int *edgei, int *edgej, int *edgek,
		bool test_bcpo, bool test_acpo, int *ierr, int *tetra_last)
  {

	int k, p, o, c;
	int it, jt, kt, idx, jdx, kdx;
	int newtetra;

	std::bitset<8> ikeep, jkeep, kkeep;

	int edge[2],tests[2];
	int itetra_touch[2],jtetra_touch[2],ktetra_touch[2];
	int position[2];

	char itetra_idx[2],jtetra_idx[2],ktetra_idx[2];

	tests[0] = 1;
	if(test_bcpo) tests[0] = -1;
	tests[1] = 1;
	if(test_acpo) tests[1] = -1;

	*ierr = 0;

/* ======================================================================================
	If itetra, jtetra or ktetra are inactive, cannot flip
 ======================================================================================*/

	if(tetra[itetra].info[1]==0 || tetra[jtetra].info[1]==0 || tetra[ktetra].info[1]==0)
	{
		*ierr = 1;
		return;
	}

/* ======================================================================================
	Store old info
 ======================================================================================*/

	ikeep = tetra[itetra].info;
	jkeep = tetra[jtetra].info;
	kkeep = tetra[ktetra].info;

/* ======================================================================================
	Define
		- itetra_touch:	indices of the two tetrahedra that share the
				faces opposite to a and b in itetra,
				respectively
		- itetra_idx:   for the two tetrahedra defined by itetra_touch,
				index position of the vertex opposite to the face
				common with itetra
		- jtetra_touch:	indices of the two tetrahedra that share the
				faces opposite to a and b in jtetra,
				respectively
		- jtetra_idx:   for the two tetrahedra defined by jtetra_touch,
				index position of the vertex opposite to the face
				common with jtetra
		- ktetra_touch:	indices of the two tetrahedra that share the
				faces opposite to a and b in ktetra,
				respectively
		- ktetra_idx:   for the two tetrahedra defined by ktetra_touch,
				index position of the vertex opposite to the face
				common with ktetra
 ======================================================================================*/

	for(int i = 0; i < 2; i++) {
		itetra_touch[i] = tetra[itetra].Neighbours[edgei[i]];
		jtetra_touch[i] = tetra[jtetra].Neighbours[edgej[i]];
		ktetra_touch[i] = tetra[ktetra].Neighbours[edgek[i]];
		itetra_idx[i]   = tetra[itetra].nindex[edgei[i]];
		jtetra_idx[i]   = tetra[jtetra].nindex[edgej[i]];
		ktetra_idx[i]   = tetra[ktetra].nindex[edgek[i]];
	}

	edge[0] = vertices[0];
	edge[1] = vertices[1];
	c       = vertices[2];
	p       = vertices[3];
	o       = vertices[4];


/* ======================================================================================
	The two new tetrahedra are going to be stored "free" space, or 
	at the end of the list
 ======================================================================================*/

	k = 0;
	while (!free.empty() && k < 2) {
		position[k] = free.top();
		free.pop();
		k++;
	}

	for (int l = k; l < 2; l++) {
		Tetrahedron t;
		t.init();
		position[l] = tetra.size();
		tetra.push_back(t);
	}
	*tetra_last = position[1];


/* ======================================================================================
	itetra, jtetra and ktetra becomes "available"; they are added to the
	"kill" list
 ======================================================================================*/

	tetra[itetra].info[1] = 0;
	tetra[jtetra].info[1] = 0;
	tetra[ktetra].info[1] = 0;

	kill.push_back(itetra);
	kill.push_back(jtetra);
	kill.push_back(ktetra);

/* ======================================================================================
	I need :
		- the two vertices that define their common edge (ab)
		 these vertices are stored in the array edge
		- the vertices c, p and o that form the new triangle
		- for each vertex in the edge (ab), define the opposing
		faces in the three tetrahedra itetra, jtetra and ktetra, and
		the tetrahedron that share these faces with itetra, jtetra and
		ktetra, respectively. This information is stored
		in three arrays, itetra_touch, jtetra_touch and ktetra_touch

	These information are given by the calling program

	For bookkeeping reasons, I always set p to be the last vertex
	of the new tetrahedra

	Now I define the two new tetrahedra: (bcop) and (acop)
	as well as their neighbours

	tetrahedron bcop : neighbours are acop, neighbour of (abop)
			   on face bpo, neighbour of (abcp) on face bcp
			   and neighbour of (abco) on face (bco)
	tetrahedron acop : neighbours are bcop, neighbour of (abop)
			   on face apo, neighbour of (abcp) on face acp
			   and neighbour of (abco) on face (aco)
 ======================================================================================*/

	for(int i = 0; i < 2; i++)
	{

		newtetra = position[i];

		k = 0;
		for(int j = 0; j < 2; j++)
		{
			if(j==i) continue;
			tetra[newtetra].Vertices[k] = edge[j];
			tetra[newtetra].Neighbours[k] = position[j];
			tetra[newtetra].nindex[k] = 0;
			k++;
		}

		tetra[newtetra].Vertices[1] = c;
		kt = ktetra_touch[i];
		kdx = ktetra_idx[i];
		tetra[newtetra].Neighbours[1] = kt;
		tetra[newtetra].nindex[1] = kdx;
		tetra[newtetra].info[3] = kkeep[2+edgek[i]];
		if(kdx != -1 && kt != -1) {
			tetra[kt].Neighbours[kdx] = newtetra;
			tetra[kt].nindex[kdx] = 1;
		}

		tetra[newtetra].Vertices[2] = o;
		it = itetra_touch[i];
		idx = itetra_idx[i];
		tetra[newtetra].Neighbours[2] = it;
		tetra[newtetra].nindex[2] = idx;
		tetra[newtetra].info[4] = ikeep[2+edgek[i]];
		if(idx != -1 && it != -1) {
			tetra[it].Neighbours[idx] = newtetra;
			tetra[it].nindex[idx] = 2;
		}

		tetra[newtetra].Vertices[3] = p;
		jt = jtetra_touch[i];
		jdx = jtetra_idx[i];
		tetra[newtetra].Neighbours[3] = jt;
		tetra[newtetra].nindex[3] = jdx;
		tetra[newtetra].info[5] = jkeep[2+edgej[i]];
		if(jdx != -1 && jt != -1) {
			tetra[jt].Neighbours[jdx] = newtetra;
			tetra[jt].nindex[jdx] = 3;
		}

		tetra[newtetra].info[1] = 1;

		if(tests[i] == 1) {
			tetra[newtetra].info[0] = 1;
		} else {
			tetra[newtetra].info[0] = 0;
		}
	}

/* ======================================================================================

	Now add the two faces of ktetra containing (co) in the link_facet 
	queue.
	Each link_facet (a triangle) is implicitly defined as the
	intersection of two tetrahedra

	link_facet:	bco	tetrahedra:	bcop and neighbour of (abco)
						on bco
	link_facet:	aco	tetrahedra:	acop and neighbour of (abco)
						on aco
 ======================================================================================*/

	for(int i = 0; i < 2; i++) {
		newtetra = position[i];
		link_facet.push(std::make_pair(newtetra, tetra[newtetra].Neighbours[3]));
		link_index.push(std::make_pair(3, tetra[newtetra].nindex[3]));
	}

  }

/* ======================================================================================
	This subroutine implements a 4->1 flip in 3D for regular triangulation

	a 4->1 flip is a transformation in which four tetrahedra are
	flipped into one tetrahedron. The two tetrahedra (abop), (bcop),
	(abcp) and (abco) shares a vertex (b) which is in the link_facet of the
	current point p added to the triangulation. After the flip, b
	is set to redundant. 
	This flip is only possible if the two edges (ab) and (bc)
	are reflex of order 3.
	We assume here that these tests have been performed and are
	true. 
	Once the flip has been performed, one tetrahedron is added
	and one new "link facet" is added to the link facet queue

	Input:
		- itetra:	index of the tetrahedra (a,b,c,p) considered
		- jtetra:	index of the tetrahedra (a,b,c,o) considered
		- ktetra:	index of the tetrahedra (a,b,o,p) considered
		- ltetra:	index of the tetrahedra (b,c,o,p) considered
		- vertices:	index of a,b,c,p,o
		- idp		index of b in (a,b,c,p)
		- jdp		index of b in (a,b,c,o)
		- kdp		index of b in (a,b,o,p)
		- ldp		index of b in (b,c,o,p)
		- test_acpo:	orientation of the 4 points (a,c,p,o)

	Output:
		- nlink_facet:	1 new link facet is added
		- link_facet:	the face of the initial tetrahedron
				(a,b,c,o) opposite to the vertex b is added
				as link facet
		- link_index:   A link_facet is a triangle defined from its
				two neighbouring tetrahedra. I store the position
				of the vertex opposite to the triangle in each
				tetrehedron in the array link_index
		- ierr:		1 if flip was not possible

 ======================================================================================*/

  void DELCX::flip_4_1(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra, 
	int itetra, int jtetra, int ktetra, int ltetra, int *ivertices,
	int idp, int jdp, int kdp, int ldp, bool test_acpo, int *ierr, int *tetra_last)
  {

	int p,o,a,b,c;
	int ishare,jshare,kshare,lshare;
	int idx,jdx,kdx,ldx;
	int test1,newtetra;

	std::bitset<8> ikeep, jkeep, kkeep, lkeep;

	*ierr = 0;

	test1 = 1;
	if(test_acpo) test1 = -1;

/* ======================================================================================
	If itetra, jtetra, ktetra, ltetra are inactive, cannot flip
 ======================================================================================*/

	if(tetra[itetra].info[1]==0 || tetra[jtetra].info[1]==0 ||
	   tetra[ktetra].info[1]==0 || tetra[ltetra].info[1]==0) {
		*ierr = 1;
		return;
	}

/* ======================================================================================
	Store "old" info
 ======================================================================================*/

	ikeep = tetra[itetra].info;
	jkeep = tetra[jtetra].info;
	kkeep = tetra[ktetra].info;
	lkeep = tetra[ltetra].info;
	
/* ======================================================================================
	Define
		- ishare:	index of tetrahedron sharing the face 
				opposite to b in itetra
		- idx		index of the vertex of ishare opposite to the
				face of ishare shared with itetra
		- jshare:	index of tetrahedron sharing the face 
				opposite to b in jtetra
		- jdx		index of the vertex of jshare opposite to the
				face of jshare shared with jtetra
		- kshare:	index of tetrahedron sharing the face 
				opposite to b in ktetra
		- kdx		index of the vertex of kshare opposite to the
				face of kshare shared with ktetra
		- lshare:	index of tetrahedron sharing the face 
				opposite to b in ltetra
		- ldx		index of the vertex of lshare opposite to the
				face of lshare shared with ltetra
 ======================================================================================*/

	ishare = tetra[itetra].Neighbours[idp];
	jshare = tetra[jtetra].Neighbours[jdp];
	kshare = tetra[ktetra].Neighbours[kdp];
	lshare = tetra[ltetra].Neighbours[ldp];

	idx = tetra[itetra].nindex[idp];
	jdx = tetra[jtetra].nindex[jdp];
	kdx = tetra[ktetra].nindex[kdp];
	ldx = tetra[ltetra].nindex[ldp];
 
	if(free.size() > 0) {
		newtetra = free.top();
		free.pop();
	} else {
		newtetra = tetra.size();
		Tetrahedron t;
		t.init();
		tetra.push_back(t);
	}

/* ======================================================================================
	itetra, jtetra, ktetra and ltetra become "available"; they
	are added to the "kill" zone
 ======================================================================================*/

	tetra[itetra].info[1] = 0;
	tetra[jtetra].info[1] = 0;
	tetra[ktetra].info[1] = 0;
	tetra[ltetra].info[1] = 0;
	kill.push_back(itetra);
	kill.push_back(jtetra);
	kill.push_back(ktetra);
	kill.push_back(ltetra);

/* ======================================================================================
	I need :
		- the vertex b that is shared by all 4 tetrahedra
		- the vertices a, c, p and o
		- for each tetrahedron, find neighbour attached to the face
		oposite to b; this information is stored in *share,
		where * can be i, j, k or l
 ======================================================================================*/

	a = ivertices[0];
	b = ivertices[1];
	c = ivertices[2];
	p = ivertices[3];
	o = ivertices[4];

/* ======================================================================================
	For bookkeeping reason, p is set to be the last vertex of the
	new tetrahedron

	Now I define the new tetrahedron: (acop)

	tetrahedron acop : neighbor of (bcop) on face cpo, neighbor of (abop)
			   on face apo, neighbor of (abcp) on face acp
			   and neighbor of (abco) on face aco
 ======================================================================================*/

	vertices[b].info[0] = 0;

	tetra[newtetra].Vertices[0] = a;
	tetra[newtetra].Neighbours[0] = lshare;
	tetra[newtetra].nindex[0] = ldx;
	tetra[newtetra].info[2] = lkeep[2+ldp];
	if(lshare != -1 && ldx != -1) {
		tetra[lshare].Neighbours[ldx] = newtetra;
		tetra[lshare].nindex[ldx] = 0;
	}
	
	tetra[newtetra].Vertices[1] = c;
	tetra[newtetra].Neighbours[1] = kshare;
	tetra[newtetra].nindex[1] = kdx;
	tetra[newtetra].info[3] = kkeep[2+kdp];
	if(kshare != -1 && kdx != -1) {
		tetra[kshare].Neighbours[kdx] = newtetra;
		tetra[kshare].nindex[kdx] = 1;
	}
	
	tetra[newtetra].Vertices[2] = o;
	tetra[newtetra].Neighbours[2] = ishare;
	tetra[newtetra].nindex[2] = idx;
	tetra[newtetra].info[4] = kkeep[2+idp];
	if(ishare != -1 && idx != -1) {
		tetra[ishare].Neighbours[idx] = newtetra;
		tetra[ishare].nindex[idx] = 2;
	}

	tetra[newtetra].Vertices[3] = p;
	tetra[newtetra].Neighbours[3] = jshare;
	tetra[newtetra].nindex[3] = jdx;
	tetra[newtetra].info[5] = kkeep[2+jdp];
	if(jshare != -1 && jdx != -1) {
		tetra[jshare].Neighbours[jdx] = newtetra;
		tetra[jshare].nindex[jdx] = 3;
	}

	tetra[newtetra].info[1] = 1;
	if(test1==1) {
		tetra[newtetra].info[0] = 1;
	} else {
		tetra[newtetra].info[0] = 0;
	}

/* ======================================================================================
	Now add one link facet : 

	link_facet:	aco	tetrahedra:	acop and neighbour of (abco)
						on aco
 ======================================================================================*/

	link_facet.push(std::make_pair(newtetra, jshare));
	link_index.push(std::make_pair(3, jdx));

  }

/* ====================================================================
	This subroutine locates the tetrahedron containing a new
	point to be added in the triangulation

	This implementation of the point location scheme
	uses a "jump-and-walk" technique: first, N active
	tetrahedra are chosen at random. The "distances" between
	these tetrahedra and the point to be added are computed,
	and the tetrahedron closest to the point is chosen as
	a starting point. The program then "walks" from that tetrahedron
	to the point, till we find a tetrahedron that contains
	the point.
	It also checks if the point is redundant in the current
	tetrahedron. It it is, the search terminates.

	Input:

	- ival:	index of the points to be located

	Output:

	- tetra_loc:	tetrahedron containing the point
	- iredundant:	flag for redundancy: 0 is not redundant,
			1 otherwise
 ==================================================================== */

  void DELCX::locate_jw(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra, 
		int ival, int *tetra_loc, int *iredundant)
  {

/* ====================================================================
	Define starting tetrahedron
 ==================================================================== */

	*iredundant = 0;
	int ntetra = tetra.size();

	if(ntetra == 1) {
		*tetra_loc = 0;
		return;
	}

	int itetra=-1;
	if(*tetra_loc < 0) {
		for(int i = ntetra-1; i >=0; i--) {
			if(tetra[i].info[1] == 1) {
				itetra = i;
				break;
			}
		}
	} else {
		itetra = *tetra_loc;
	}

	int a, b, c, d, iorient;
	bool test_in, test_red;
	int idx;

	do {

		a = tetra[itetra].Vertices[0];
		b = tetra[itetra].Vertices[1];
		c = tetra[itetra].Vertices[2];
		d = tetra[itetra].Vertices[3];
		iorient = -1;
		if(tetra[itetra].info[0]==1) iorient = 1;

		inside_tetra(vertices, ival, a, b, c, d, 
		iorient, &test_in, &test_red, &idx);

		if(!test_in) itetra = tetra[itetra].Neighbours[idx];

	} while(!test_in);

	*tetra_loc = itetra;

	if(test_red) *iredundant = 1;

  }

/* ==============================================================================
	After a point has been inserted, this subroutine goes over the 
	link_facet list to restore regularity. When a link_facet is found
	non_regular and "flippable" (see below), the program attempts
	to flip it. If the flip is successful, new link_facets are added
	on the queue.
	The subroutine ends when the link facet is empty
 ============================================================================== */

  void DELCX::flip(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra)
  {
	int ii, ij;
	int idxi, idxj, idxk, idxl;
	int p, o, a, b, c;
	int itetra, jtetra, idx_p, idx_o;
	int itest_abcp;
	int ierr, tetra_last;
	int ireflex, iflip, iorder;
	int ifind, tetra_ab, tetra_ac, tetra_bc; 
	int idx_a = 0, idx_b = 0, idx_c = 0;
	int ia, ib, ic;

	std::pair<int, int> facet, index;
	int facei[3], facej[3], edgei[2], edgej[2], edgek[2];
	int edge_val[3][2];
	int tetra_flip[3], list_flip[3];
	int vert_flip[5];

	bool convex, regular, test_abpo, test_bcpo, test_capo; 
	bool test, test_abpc, test_bcpa, test_acpo, test_acpb;
	bool test_or[3][2];

/* ==============================================================================
	Go over all link facets
 ============================================================================== */

	while(!link_facet.empty()) {

/* ==============================================================================
		First defined the two tetrahedra that contains the link facet as
		itetra and jtetra
 ============================================================================== */

		facet = link_facet.front();
		index = link_index.front();
		link_facet.pop();
		link_index.pop();

		itetra = facet.first;
		jtetra = facet.second;
		idx_p = index.first;
		idx_o = index.second;

/* ==============================================================================
		If the link facet is on the convex hull, discard
 ============================================================================== */

		if(itetra==-1 || jtetra == -1) continue;

/* ==============================================================================
		If these tetrahedra have already been discarded, discard this
		link facet
 ============================================================================== */

		if(tetra[itetra].info[1] == 0) {
			if(tetra[jtetra].info[1] == 0) {
				continue;
			} else {
				itetra = tetra[jtetra].Neighbours[idx_o];
				idx_p = tetra[jtetra].nindex[idx_o];
			}
		}

		if(tetra[jtetra].info[1] == 0) {
			jtetra = tetra[itetra].Neighbours[idx_p];
			idx_o = tetra[itetra].nindex[idx_p];
		}

/* ==============================================================================
		Let us define the vertices of the two tetrahedra:
		itetra:		a,b,c,p
		jtetra:		a,b,c,o
 ============================================================================== */

		a = tetra[itetra].Vertices[0];
		b = tetra[itetra].Vertices[1];
		c = tetra[itetra].Vertices[2];
		p = tetra[itetra].Vertices[3];

		o = tetra[jtetra].Vertices[idx_o];

		itest_abcp = -1;
		if(tetra[itetra].info[0]==1) itest_abcp = 1;

/* ==============================================================================
		Check for local regularity (and convexity, at very little
		extra cost)
 ============================================================================== */

		regular_convex(vertices, a, b, c, p, o, itest_abcp, &regular,
		&convex, &test_abpo, &test_bcpo, &test_capo); 

/* ==============================================================================
		if the link facet is locally regular, discard
 ============================================================================== */

		if(regular) continue;

/* ==============================================================================
		Define neighbors of the facet on itetra and jtetra
 ============================================================================== */

		define_facet(tetra, itetra, jtetra, idx_o, facei, facej);

		test_abpc = (itest_abcp != 1);

/* ==============================================================================
		After discarding the trivial case, we now test if the tetrahedra
		can be flipped. 

		At this stage, I know that the link facet is not locally
		regular. I still don t know if it is "flippable"

		I first check if {itetra} U {jtetra} is convex. If it is, I
		perform a 2-3 flip (this is the convexity test performed
		at the same time as the regularity test)
 ============================================================================== */

		if(convex) {
			vert_flip[0] = a;
			vert_flip[1] = b;
			vert_flip[2] = c;
			vert_flip[3] = p;
			vert_flip[4] = o;
			flip_2_3(tetra, itetra, jtetra, vert_flip, facei, facej,
			test_abpo, test_bcpo, test_capo, &ierr, &tetra_last);

			continue;
		}

/* ==============================================================================
		The union of the two tetrahedra is not convex...
		I now check the edges of the triangle in the link facet, and
		check if they are "reflexes" (see definition in Edelsbrunner and
		Shah, Algorithmica (1996), 15:223-241)
 ============================================================================== */

		ireflex = 0;
		iflip = 0;

/* ==============================================================================
		First check edge (ab): 
		- (ab) is reflex iff o and c lies on opposite sides of
		the hyperplane defined by (abp). We therefore test the
		orientation of (abpo) and (abpc): if they differ (ab)
		is reflex
		- if (ab) is reflex, we test if it is of degree 3.
		(ab) is of degree 3 if it is shared by 3 tetrahedra,
		namely (abcp), (abco) and (abpo). The first two are itetra
		and jtetra, so we only need to check if (abpo) exists.
		since (abpo) contains p, (abp) should then be a link facet
		of p, so we test all tetrahedra that define link facets
 ============================================================================== */

		if(test_abpo != test_abpc) {

			ireflex++;
			find_tetra(tetra, itetra, 2, a, b, o, &ifind,
			&tetra_ab, &idx_a, &idx_b);

			if(ifind==1) {
				tetra_flip[iflip] = tetra_ab;
				list_flip[iflip] = 0;
				edge_val[iflip][0] = idx_a;
				edge_val[iflip][1] = idx_b;
				test_or[iflip][0] = test_bcpo;
				test_or[iflip][1] = !test_capo;
				iflip++;
			}
		}

/* ==============================================================================
		Now check edge (ac): 
		- (ac) is reflex iff o and b lies on opposite sides of
		the hyperplane defined by (acp). We therefore test the
		orientation of (acpo) and (acpb): if they differ (ac)
		is reflex
		- if (ac) is reflex, we test if it is of degree 3.
		(ac) is of degree 3 if it is shared by 3 tetrahedra,
		namely (abcp), (abco) and (acpo). The first two are itetra
		and jtetra, so we only need to check if (acpo) exists.
		since (acpo) contains p, (acp) should then be a link facet
		of p, so we test all tetrahedra that define link facets
 ============================================================================== */

		test_acpo = !test_capo;
		test_acpb = !test_abpc;

		if(test_acpo != test_acpb) {

			ireflex++;
			find_tetra(tetra, itetra, 1, a, c, o, &ifind,
			&tetra_ac, &idx_a, &idx_c);

			if(ifind==1)  {
				tetra_flip[iflip] = tetra_ac;
				list_flip[iflip] = 1;
				edge_val[iflip][0] = idx_a;
				edge_val[iflip][1] = idx_c;
				test_or[iflip][0] = !test_bcpo;
				test_or[iflip][1] = test_abpo;
				iflip++;
			}
		}

/* ==============================================================================
		Now check edge (bc): 
		- (bc) is reflex iff o and a lies on opposite sides of
		the hyperplane defined by (bcp). We therefore test the
		orientation of (bcpo) and (bcpa): if they differ (bc)
		is reflex
		- if (bc) is reflex, we test if it is of degree 3.
		(bc) is of degree 3 if it is shared by 3 tetrahedra,
		namely (abcp), (abco) and (bcpo). The first two are itetra
		and jtetra, so we only need to check if (bcpo) exists.
		since (bcpo) contains p, (bcp) should then be a link facet
		of p, so we test all tetrahedra that define link facets
 ============================================================================== */

		test_bcpa = test_abpc;
		if(test_bcpo != test_bcpa) {

			ireflex++;
			find_tetra(tetra, itetra, 0, b, c, o, &ifind,
			&tetra_bc, &idx_b, &idx_c);

			if(ifind==1)  {
				tetra_flip[iflip] = tetra_bc;
				list_flip[iflip] = 2;
				edge_val[iflip][0] = idx_b;
				edge_val[iflip][1] = idx_c;
				test_or[iflip][0] = test_capo;
				test_or[iflip][1] = !test_abpo;
				iflip++;
			}
		}

		if(ireflex != iflip) continue;

		if(iflip==1) {

/* ==============================================================================
			Only one edge is "flippable": we do a 3-2 flip
 ============================================================================== */

			iorder = list_flip[0];
			ia = table32[iorder][0];
			ib = table32[iorder][1];
			ic = table32[iorder][2];
			vert_flip[ia] = a;
			vert_flip[ib] = b;
			vert_flip[ic] = c;
			vert_flip[3] = p;
			vert_flip[4] = o;
			ia = table32_2[iorder][0];
			ib = table32_2[iorder][1];
			edgei[0] = ia;
			edgei[1] = ib;
			edgej[0] = facej[ia];
			edgej[1] = facej[ib];
			edgek[0] = edge_val[0][0];
			edgek[1] = edge_val[0][1];
			flip_3_2(tetra, itetra, jtetra, tetra_flip[0],
			vert_flip, edgei, edgej, edgek, test_or[0][0],
			test_or[0][1], &ierr, &tetra_last);

		} else if (iflip==2) {

/* ==============================================================================
			In this case, one point is redundant: the point common to
			the two edges that can be flipped. We then perform a 4-1
			flip
 ============================================================================== */

			iorder = list_flip[0] + list_flip[1] -1;
			vert_flip[table41[iorder][0]] = a;
			vert_flip[table41[iorder][1]] = b;
			vert_flip[table41[iorder][2]] = c;
			vert_flip[3] = p;
			vert_flip[4] = o;
			ii = table41_2[iorder][0];
			ij = table41_2[iorder][1];
			idxi = iorder;
			idxj = facej[iorder];
			idxk = edge_val[0][ii];
			idxl = edge_val[1][ij];

			if(iorder==0) {
				test = test_bcpo;
			} else if (iorder==1) {
				test = !test_capo;
			} else {
				test = test_abpo;
			}

			flip_4_1(vertices, tetra, itetra, jtetra, tetra_flip[0],
			tetra_flip[1], vert_flip, idxi, idxj, idxk, idxl, test,
			&ierr, &tetra_last);

		} else {
			std::cout << "Problem.... three edges flippable!!" << std::endl;
			exit(1);
		}

	}

	for(int i = 0; i < kill.size(); i++) {
		free.push(kill[i]);
	}
	kill.clear();

  }

/* ===============================================================================
 Define_facet

	A triangle (or facet) is defined by the intersection of two 
	tetrahedra itetra and jtetra
	If we know the position of its three vertices in the first
	tetrahedron (in fact the first three vertices a,b and c
	of itetra), we need to find the indices of these vertices
	in the second tetrahedron.
	This routine also stores information about the neighbours
	of the two tetrahedra considered

	The vertices are called a,b,c,p, and o, where (abc) is the
	common facet

	Input:
		- itetra:	index of the tetrahedra (a,b,c,p) considered
		- jtetra:	index of the tetrahedra (a,b,c,o) considered
		- idx_o:	position of o in the vertices of jtetra

	Output:
		- jtouch	jtouch(i) is the tetrahedron sharing
				the face opposite to i in tetrahedron jtetra
		- jdx		jdx(i) is the vertex of jtouch(i) opposite
				to the face shared with jtetra
 =============================================================================== */

  void DELCX::define_facet(std::vector<Tetrahedron>& tetra, int itetra, int jtetra, 
		int idx_o, int *facei, int *facej)
  {

	int ia, ib, ie, ig;
	int k;

/* ===============================================================================
	I need to :
		- find the three vertices that define their common face
		 these vertices are stored in the array triangle
		- find the vertices p and o

	To define the common face of the two tetrahedra itetra and jtetra,
	I look at the neighbours of itetra : one of them is jtetra!
	This also provides p. The same procedure is repeated for jtetra,
	to get o
 =============================================================================== */

	for(int i = 0; i < 3; i++) facei[i] = i;

	ia = tetra[itetra].Vertices[0];
	for(int i = 0; i < 3; i++) {
		k = other[idx_o][i];
		ie = tetra[jtetra].Vertices[k];
		if(ia == ie) {
			facej[0] = k;
			break;
		}
	}

	ib = tetra[itetra].Vertices[1];
	ie = other2[idx_o][facej[0]][0];
	ig = other2[idx_o][facej[0]][1];
	if(ib == tetra[jtetra].Vertices[ie]) {
		facej[1] = ie;
		facej[2] = ig;
	} else {
		facej[1] = ig;
		facej[2] = ie;
	}
  }

/* ===============================================================================
	Find_tetra
	This subroutine tests if four given points form an existing
	tetrahedron in the current Delaunay
 =============================================================================== */

  void DELCX::find_tetra(std::vector<Tetrahedron>& tetra, int itetra, 
	int idx_c, int a, int b, int o, int *ifind, int *tetra_loc, 
	int *idx_a, int *idx_b)
  { 

/* ===============================================================================
	We are testing if tetrahedron (abpo) exists. If it exists, it is
	a neighbour of abcp, on the face opposite to vertex c.
	We test that tetrahedron and see if it contains o
 =============================================================================== */

	int ot, otx, otest;

	ot = tetra[itetra].Neighbours[idx_c];
	otx = tetra[itetra].nindex[idx_c];
	otest = tetra[ot].Vertices[otx];

	if(otest == o) {
		*ifind = 1;
		*tetra_loc = ot;

/* ===============================================================================
		We found the tetrahedron, let us define the position
		of a and b in this tetrahedron
 =============================================================================== */

		for(int i = 0; i < 4; i++) {
			if(tetra[*tetra_loc].Vertices[i] == a) {
				*idx_a = i;
			} else if(tetra[*tetra_loc].Vertices[i] == b) {
				*idx_b = i;
			}
		}
	} else {
		*ifind = 0;
	}

  }

/* ===============================================================================
	Remove_inf
	This subroutine sets to 0 the status of tetrahedron that
	contains infinite points
 =============================================================================== */

  void DELCX::remove_inf(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra)
  {
	int a, b, c, d;
	int ntetra = tetra.size();

	for(int i = 0; i < ntetra; i++)
	{
		if(tetra[i].info[1]==0) continue;

		a = tetra[i].Vertices[0];
		b = tetra[i].Vertices[1];
		c = tetra[i].Vertices[2];
		d = tetra[i].Vertices[3];

		if(a<4 || b<4 || c<4 || d<4) {
			tetra[i].info[2] = 1;
			tetra[i].info[1] = 0;
			if(a < 4) mark_zero(tetra, i, 0);
			if(b < 4) mark_zero(tetra, i, 1);
			if(c < 4) mark_zero(tetra, i, 2);
			if(d < 4) mark_zero(tetra, i, 3);
		}
	}

	for(int i = 0; i < 4; i++) {
		vertices[i].info[0] = 0;
	}
  }
/* ===============================================================================
	This subroutine marks the tetrahedron that touches
	a tetrahedron with infinite point as part of the
	convex hull (i.e. one of its neighbor is 0)
 =============================================================================== */

  void DELCX::mark_zero(std::vector<Tetrahedron>& tetra, int itetra, int ivertex)
  {
	int jtetra, jvertex;

	jtetra = tetra[itetra].Neighbours[ivertex];

	if(jtetra != -1) {
		jvertex = tetra[itetra].nindex[ivertex];
		tetra[jtetra].Neighbours[jvertex] = -1;
	}
  }

/* ====================================================================
	Peel
	This subroutine removes the flat tetrahedra at the boundary 
	of the DT
 ==================================================================== */

  void DELCX::peel(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra)
  {

	int ia, ib, ic, id;
	int k, l;
	int res;
	int ntetra = tetra.size();
	double coorda[3], coordb[3], coordc[3], coordd[3];

	for(int i = 0; i < ntetra; i++) {

		if(tetra[i].info[1]==0) continue;

		bool itest = false;
		for(int j = 0; j < 4; j++) {
			if(tetra[i].Neighbours[j]==-1) itest = true;
		}

		if(!itest) continue; 

/* ====================================================================
		This is a tetrahedron at the boundary: we test
		if it is flat, i.e. if its volume is 0
 ==================================================================== */

		ia = tetra[i].Vertices[0];
		ib = tetra[i].Vertices[1];
		ic = tetra[i].Vertices[2];
		id = tetra[i].Vertices[3];

		for(int j = 0; j < 3; j++) {
			coorda[j] = vertices[ia].Coordinates[j];
			coordb[j] = vertices[ib].Coordinates[j];
			coordc[j] = vertices[ic].Coordinates[j];
			coordd[j] = vertices[id].Coordinates[j];
		}

		double vol = tetra_vol(coorda, coordb, coordc, coordd);

		if(std::abs(vol) < eps) {
			sos.minor4_gmp(coorda, coordb, coordc, coordd, &res);
			if(res == 0) {
				tetra[i].info[2] = 1;
			}
		}
	}

/* ====================================================================
	Now we remove those flat tetrahedra, and update the links
	to their neighbours
 ==================================================================== */

	for(int i = 0; i < ntetra; i++) {
		if(tetra[i].info[2]==1) {
			if(tetra[i].info[1]==1) {
				tetra[i].info[1] = 0;
				for(int j = 0; j < 4; j++) {
					k = tetra[i].Neighbours[j];
					if(k != -1) {
						l = tetra[i].nindex[j];
						tetra[k].Neighbours[l] = -1;
					}
				}
			}
		}
	}

  }

/* ====================================================================
	Computes the volume of a tetrahedron
 ==================================================================== */

   double DELCX::tetra_vol(double *a, double *b, double *c, double *d)
  {
	double vol;
	double ad[3], bd[3], cd[3];
	double Sbcd[3];

/* ====================================================================
	The volume of the tetrahedron is proportional to:

	vol = det | a(1)  a(2)  a(3)  1|
		  | b(1)  b(2)  b(3)  1|
		  | c(1)  c(2)  c(3)  1|
		  | d(1)  d(2)  d(3)  1|

	After substracting the last row from the first 3 rows, and
	developping with respect to the last column, we obtain:

	vol = det | ad(1)  ad(2)  ad(3) |
		  | bd(1)  bd(2)  bd(3) |
		  | cd(1)  cd(2)  cd(3) |

	where ad(i) = a(i) - d(i), ...
 ==================================================================== */

	for(int i = 0; i < 3; i++) {
		ad[i] = a[i] - d[i];
		bd[i] = b[i] - d[i];
		cd[i] = c[i] - d[i];
	}

	Sbcd[2] = bd[0]*cd[1] - cd[0]*bd[1];
	Sbcd[1] = bd[0]*cd[2] - cd[0]*bd[2];
	Sbcd[0] = bd[1]*cd[2] - cd[1]*bd[2];

	vol = ad[0]*Sbcd[0] - ad[1]*Sbcd[1] + ad[2]*Sbcd[2];

	return vol;
  }

/* ====================================================================
	reorder_tetra

	reorders the vertices of a list of tetrahedron,
	such that now the indices are in increasing order
 ==================================================================== */

   void DELCX::reorder_tetra(std::vector<Tetrahedron>& tetra)
   {

	int ntetra = tetra.size();
	int vert[4], idx[4], neighbor[4];
	int n = 4;
	int nswap;
	char nidx[4];
	bool nsurf[4];

	for(int i = 0; i < ntetra; i++)
	{
		if(tetra[i].info[1]==0) continue;

		for(int j = 0; j < 4; j++) {
			vert[j] = tetra[i].Vertices[j];
		}
		mysort.sort4_sign(vert, idx, &nswap, n);

		for(int j = 0; j < 4; j++) {
			neighbor[j] = tetra[i].Neighbours[idx[j]];
			nidx[j] = tetra[i].nindex[idx[j]];
			std::string s = tetra[i].info.to_string();
			nsurf[j] = (s[2+idx[j]]=='1');
			if(neighbor[j] != -1) {
				tetra[neighbor[j]].nindex[nidx[j]]=j;
			}
		}

		for(int j = 0; j < 4; j++) {
			tetra[i].Vertices[j] = vert[j];
			tetra[i].Neighbours[j] = neighbor[j];
			tetra[i].nindex[j] = nidx[j];
			tetra[i].info.set(2+j,nsurf[j]);
		}

		if(nswap==-1) {
			if(tetra[i].info[0] == 0) {
				tetra[i].info[0] = 1;
			} else {
				tetra[i].info[0] = 0;
			}
		}
	}
  }
/*===========================================================================================
	DelaunayEdges

	This procedure generates the list of edges in the Delaunay triangulation
 
 ========================================================================================== */

  void DELCX::delaunayEdges(std::vector<Tetrahedron>& tetra, std::vector<std::pair<int, int> >& edges)
  {
 
	int face_info[6][2] = {
		{ 0, 1},
		{ 0, 2},
		{ 0, 3},
		{ 1, 2},
		{ 1, 3},
		{ 2, 3},
	};

	int pair[6][2] = {
		{ 2, 3},
		{ 1, 3},
		{ 1, 2},
		{ 0, 3},
		{ 0, 2},
		{ 0, 1},
	};

/* ============================================================================================
	define mask on edges of tetrahedra, to check if already seen
 ============================================================================================ */

	int ntetra = tetra.size();
	std::bitset<6> *tetra_mask = new std::bitset<6>[ntetra];

	std::bitset<6> zero(std::string("000000"));
	for(int i = 0; i < ntetra; i++) {
		tetra_mask[i] = zero;
	}

/* ============================================================================================
	define mask on edges of tetrahedra, to check if already seen
 	loop over all tetrahedron: if it belongs to the Delaunay triangulation,
 	check its edges; include in edge list if not seen before
 ============================================================================================ */
 
	edges.clear();

	int i, j;
	int trig1, trig2;
	int jtetra, ktetra;
	int npass;
	bool done;
	int trig_in, trig_out;
	int triga, trigb;
	int ipair;
	std::pair<int, int> p;

	for(int idx = 0; idx < ntetra; idx++) {

		if(tetra[idx].info[1]==0) continue;
 
/* ============================================================================================
 	    Check all six edges
 ============================================================================================ */
 					
		for(int iedge = 0; iedge < 6; iedge++) {
 
/* ============================================================================================
 		If this edge has already been considered (from another tetrahedron),
 		discard
 ============================================================================================ */
 
			if(tetra_mask[idx][iedge] == 1) continue;
 
/* ============================================================================================
 			iedge is the edge number in the tetrahedron idx, with:
 			iedge = 1		(c,d)
 			iedge = 2		(b,d)
 			iedge = 3		(b,c)
 			iedge = 4		(a,d)
 			iedge = 5		(a,c)
 			iedge = 6		(a,b)
 		
 			Define indices of the two vertices of the edge
 ============================================================================================ */
 
			i = tetra[idx].Vertices[pair[iedge][0]];
			j = tetra[idx].Vertices[pair[iedge][1]];
			p = std::make_pair(i,j);
 
			edges.push_back(p);
 
/* ============================================================================================
 			trig1 and trig2 are the two faces of idx that share
 			iedge
 ============================================================================================ */
 
			trig1 = face_info[iedge][0];
			trig2 = face_info[iedge][1];
 
/* ============================================================================================
 			Now we look at the star of the edge:
 ============================================================================================ */
 
			ktetra = idx;
			npass = 0;
			trig_out = trig1;
			jtetra = tetra[ktetra].Neighbours[trig1];
			done = false;

			while(!done) {

/* ============================================================================================
 				Leave this side of the star if we hit the convex hull
 				in this case, the edge is not buried
 ============================================================================================ */

				if(jtetra==-1) {
					if(npass==1) done = true;
					npass++;
					ktetra = idx;
					trig_out = trig2;
					jtetra = tetra[ktetra].Neighbours[trig_out];
				} else {
					if(jtetra==idx) done = true;
					ipair = findEdge(tetra[jtetra], i, j);
					tetra_mask[jtetra][ipair] = 1;
					trig_in = tetra[ktetra].nindex[trig_out];
					triga = face_info[ipair][0];
					trigb = face_info[ipair][1];
					trig_out = triga;
					if(trig_in == triga) trig_out = trigb;
					ktetra = jtetra;
					jtetra = tetra[ktetra].Neighbours[trig_out];
				}
			}

		}
	}	
 
/* ============================================================================================
 	Sort list of all edges in increasing order
 ============================================================================================ */
 
	std::sort(edges.begin(), edges.end());

	delete [] tetra_mask;

 }

/* ============================================================================================
 	Given two vertices of a tetrahedron, find the index of the edge they form
 ============================================================================================ */
 
  int DELCX::findEdge(Tetrahedron t, int i1, int j1)
  {

	int ipair;

	if(i1==t.Vertices[0]) {
		if(j1==t.Vertices[1]) {
			ipair = 5;
		} else if(j1==t.Vertices[2]) {
			ipair = 4;
		} else {
			ipair = 3;
		}
	} else if(i1==t.Vertices[1]) {
		if(j1==t.Vertices[2]) {
			ipair = 2;
		} else {
			ipair = 1;
		}
	} else {
		ipair = 0;
	}

	return ipair;

   }
/* ============================================================================================
  Add bogus points as needed so that we have at least 4 points
 ============================================================================================ */
 
  void DELCX::addBogus(int npoints, double *coord, double *radii, double *bcoord, double *brad)
  {

	if(npoints > 3) return;

	int np = 4 - npoints;
	memset(bcoord, 0, 3*np*sizeof(double));
	Vector c, c1, c2, c3;
	Vector u1, v1, w1;
	double Rmax, d, d1, d2, d3;

	if(npoints==1) {
		Rmax = radii[0];
		for(int i = 0; i < np; i++) {
			bcoord[3*i+i] = coord[i] + 3*Rmax;
			brad[i] = Rmax/20;
		}
	} else if(npoints==2) {
		Rmax = std::max(radii[0], radii[1]);
		c1[0] = coord[0]; c1[1] = coord[1]; c1[2] = coord[2];
		c2[0] = coord[3]; c2[1] = coord[4]; c2[2] = coord[5];
		c = 0.5*(c1+c2);
		u1 = c2-c1; 
		if((u1[2]!=0) || (u1[0]!=-u1[1])) {
			v1[0] = u1[2]; v1[1] = u1[2]; v1[2] = -u1[0]-u1[2];
		} else {
			v1[0] = -u1[1]-u1[2]; v1[1] = u1[0]; v1[2] = u1[0];
		}
		w1 = u1^v1;
		d = u1.norm();
		for(int i = 0; i < 3; i++) {
			bcoord[i] = c[i] + (2*d+3*Rmax)*v1[i];
			bcoord[i+3] = c[i] + (2*d+3*Rmax)*w1[i];
		}
		brad[0] = Rmax/20; brad[1] = Rmax/20;
	} else {
		Rmax = std::max(std::max(radii[0], radii[1]), radii[2]);
		c1[0] = coord[0]; c1[1] = coord[1]; c1[2] = coord[2];
		c2[0] = coord[3]; c2[1] = coord[4]; c2[2] = coord[5];
		c3[0] = coord[6]; c3[1] = coord[7]; c3[2] = coord[8];
		c = (c1+c2+c3)/3;
		u1 = c2-c1; v1 = c3-c1;
		w1 = u1^v1;
		d1 = u1.norm(); d2 = v1.norm(); d3 = (c3-c2).norm(); d = std::max(std::max(d1,d2),d3);
		for(int i = 0; i < 3; i++) {
			bcoord[i] = c[i] + (2*d+3*Rmax)*w1[i];
		}
		brad[0] = Rmax/20;
	}

  }

#endif
