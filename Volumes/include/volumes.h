/* ====================================================================

This file contains a series of procedures used to compute the
intrinsic volumes of a union of balls, i.e. its surface area, volume,
mean integrated curvature, and Gaussian curvature, and optionally their derivatives
with respect to the coordinates of the centers of the ball.

 ==================================================================== */

#ifndef VOLUMES_H
#define VOLUMES_H

#include <vector>
#include "tetra.h"
#include "Vector.h"
#include "Tetrahedron.h"
#include "gauss_corners.h"
#include "Edge.h"
#include "Face.h"

TETRAGEOM tetrageom;
GAUSSCORNER gauss;
 
/* ====================================================================
   class
 ==================================================================== */

  class VOLUMES {

	public:

		void ball_dvolumes(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra,
			std::vector<Edge>& edges, std::vector<Face>& faces,
			double *WSurf, double *WVol, double *WMean, double *WGauss,
			double *Surf, double *Vol, double *Mean, double *Gauss,
			double *ballwsurf, double *ballwvol, double *ballwmean, double *ballwgauss,
			double *dsurf_coord, double *dvol_coord, double *dmean_coord,
			double *dgauss_coord, int option);

	private:

		double distance2(std::vector<Vertex>& vertices, int n1, int n2);

		void twosphere_info(double ra, double ra2, double rb, double rb2,
			double rab, double rab2, double *surfa, double *surfb,
			double *vola, double *volb, double *r, double *phi, double *l);

		void twosphere_dinfo(double ra, double ra2, double rb, double rb2,
			double rab, double rab2, double *surfa, double *surfb,
			double *vola, double *volb, double *r, double *phi, double *l,
			double *dsurfa, double *dsurfb, double *dvola, double *dvolb, 
			double *dr, double *dphi, double *dl, int option);

		void threesphere_dvol(double ra, double rb,double rc, double ra2,
			double rb2, double rc2, double rab, double rac, double rbc,
			double rab2, double rac2, double rbc2, double *angle, double deriv[6][3],
			double *surfa, double *surfb, double *surfc, double *vola, double *valb,
			double *volc, double *dsurfa, double *dsurfb, double *dsurfc, 
			double *dvola, double *dvolb, double *dvolc, int option);

		double plane_dist(double ra2, double rb2, double rab2);

		double pi = M_PI;
		double twopi = 2.0*pi;

		double coef_E = 1.;
		double coef_F = 2.;

  };

/* ====================================================================
	ball_dvolumes
 ==================================================================== */

   void VOLUMES::ball_dvolumes(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra,
	std::vector<Edge>& edges, std::vector<Face>& faces,
	double *WSurf, double *WVol, double *WMean, double *WGauss,
	double *Surf, double *Vol, double *Mean, double *Gauss,
	double *ballwsurf, double *ballwvol, double *ballwmean, double *ballwgauss,
	double *dsurf_coord, double *dvol_coord, double *dmean_coord, double *dgauss_coord, int option)
   {

	int ia, ib, ic, id;
	int e1, e2, e3;

	double  ra, ra2, rb, rb2, rc, rc2, rd, rd2;
	double	rab, rac, rad, rbc, rbd, rcd;
	double	rab2, rac2, rad2, rbc2, rbd2, rcd2;
	double	val, val1S, val2S, val3S, val4S;
	double  val1M1, val2M1, val3M1, val4M1;
	double  val1G1;
	double  d1, d2, d3, d4;
	double  val1V, val2V, val3V, val4V;
	double	coefval, coefvalS;
	double	surfa, surfb, surfc;
	double	gaussa, gaussb, gaussc;
	double	vola, volb, volc, vold;
	double  r, phi, l, dr, dphi, dl;
	double  coefaS, coefbS, coefcS, coefdS;
	double  coefaV, coefbV, coefcV, coefdV;
	double  coefaM, coefbM, coefcM, coefdM;
	double  coefaG, coefbG, coefcG, coefdG;
	double	dsurfa2, dsurfb2;
	double	dvola2, dvolb2;

	double	u[3];
	double  angle[6], cosine[6], sine[6];
	double	deriv[6][6], deriv2[6][3], dg[3][3];
	double	dsurfa3[3], dsurfb3[3], dsurfc3[3];
	double	dvola3[3], dvolb3[3], dvolc3[3];
	double	dvola[6], dvolb[6], dvolc[6], dvold[6];

	int edge_list[6];

	int nedges = edges.size();
	int nvertices = vertices.size();
	int nfaces = faces.size();
	int ntetra = tetra.size();

/* ====================================================================
	Initialize results arrays
 ==================================================================== */

	*WSurf = 0; *Surf  = 0; *WVol  = 0; *Vol   = 0;
	*WMean = 0; *Mean  = 0; *WGauss= 0; *Gauss = 0;

	memset(ballwsurf,  0, nvertices*sizeof(double));
	memset(ballwvol,   0, nvertices*sizeof(double));
	memset(ballwmean,  0, nvertices*sizeof(double));
	memset(ballwgauss, 0, nvertices*sizeof(double));

/* ====================================================================
	Initialize edge and vertex info
 ==================================================================== */

	for(int i = 0; i < nedges; i++) {

		edges[i].gamma = 1.;
		edges[i].sigma = 1.;

		ia = edges[i].Vertices[0];
		ib = edges[i].Vertices[1];

		if(vertices[ia].status==0 || vertices[ib].status==0) continue;

		ra = vertices[ia].Radius; ra2 = ra*ra;
		rb = vertices[ib].Radius; rb2 = rb*rb;

		coefaS = vertices[ia].CoefS; coefbS = vertices[ib].CoefS;
		coefaV = vertices[ia].CoefV; coefbV = vertices[ib].CoefV;
		coefaM = vertices[ia].CoefM; coefbM = vertices[ib].CoefM;
		coefaG = vertices[ia].CoefG; coefbG = vertices[ib].CoefG;

		rab2 = distance2(vertices, ia, ib); rab = std::sqrt(rab2);

		twosphere_info(ra, ra2, rb, rb2, rab, rab2, &surfa, &surfb,
		&vola, &volb, &r, &phi, &l);

		edges[i].Length = rab;
		edges[i].Surf   = (coefaS*surfa + coefbS*surfb)/twopi;
		edges[i].Vol    = (coefaV*vola + coefbV*volb)/twopi;
		edges[i].CoefM1 = (coefaM*surfa/ra + coefbM*surfb/rb)/twopi;
		edges[i].CoefM2  = coef_E*(coefaM+coefbM)*r*phi/2;
		edges[i].CoefG1 = (coefaG*surfa/ra2 + coefbG*surfb/rb2)/twopi;
		edges[i].CoefG2  = coef_E*(coefaG+coefbG)*l/2;
		edges[i].dsurf  = 0;
		edges[i].dvol   = 0;
		edges[i].dmean  = 0;
		edges[i].dgauss = 0;

	}

	for(int i = 0; i < nvertices; i++) {
		vertices[i].gamma = 1.;
	}

/* ====================================================================
	Contributions of four overlapping spheres:

	We are using the weighted inclusion-exclusion formula:
	Each tetrahedron in the Alpha Complex only contributes to the weight of each
	its edges and each of its vertices
 ==================================================================== */


	for(int idx = 0; idx < ntetra; idx++) {

		if(tetra[idx].info[6]==0) continue;

		ia = tetra[idx].Vertices[0];
		ib = tetra[idx].Vertices[1];
		ic = tetra[idx].Vertices[2];
		id = tetra[idx].Vertices[3];

		if(vertices[ia].status==0 || vertices[ib].status==0
		|| vertices[ic].status==0 || vertices[id].status==0) continue;

		ra = vertices[ia].Radius; ra2 = ra*ra;
		rb = vertices[ib].Radius; rb2 = rb*rb;
		rc = vertices[ic].Radius; rc2 = rc*rc;
		rd = vertices[id].Radius; rd2 = rd*rd;

		coefaS = vertices[ia].CoefS; coefaV = vertices[ia].CoefV; 
		coefaM = vertices[ia].CoefM; coefaG = vertices[ia].CoefG;
		coefbS = vertices[ib].CoefS; coefbV = vertices[ib].CoefV; 
		coefbM = vertices[ib].CoefM; coefbG = vertices[ib].CoefG;
		coefcS = vertices[ic].CoefS; coefcV = vertices[ic].CoefV; 
		coefcM = vertices[ic].CoefM; coefcG = vertices[ic].CoefG;
		coefdS = vertices[id].CoefS; coefdV = vertices[id].CoefV; 
		coefdM = vertices[id].CoefM; coefdG = vertices[id].CoefG;

		for(int iedge = 0; iedge < 6; iedge++) {

/* 		====================================================================
		iedge is the edge number in the tetrahedron idx, with:
		iedge = 1		(c,d)
		iedge = 2		(b,d)
		iedge = 3		(b,c)
		iedge = 4		(a,d)
		iedge = 5		(a,c)
		iedge = 6		(a,b)
 		==================================================================== */

			edge_list[5-iedge] = tetra[idx].info_edge[iedge];

		}

		rab = edges[edge_list[0]].Length; rab2 = rab*rab;
		rac = edges[edge_list[1]].Length; rac2 = rac*rac;
		rad = edges[edge_list[2]].Length; rad2 = rad*rad;
		rbc = edges[edge_list[3]].Length; rbc2 = rbc*rbc;
		rbd = edges[edge_list[4]].Length; rbd2 = rbd*rbd;
		rcd = edges[edge_list[5]].Length; rcd2 = rcd*rcd;

/* 		====================================================================
	   	Characterize tetrahedron (A,B,C,D)
 		==================================================================== */

		if(option==0) {
			tetrageom.tetra_dihed(rab2, rac2, rad2, rbc2, rbd2, rcd2,
			angle, cosine, sine);
		} else {
			tetrageom.tetra_dihed_der(rab2, rac2, rad2, rbc2, rbd2, rcd2,
			angle, cosine, sine, deriv);
		}

/* 		====================================================================
           	Add fraction of tetrahedron that "belongs" to each ball
 		==================================================================== */

		tetrageom.tetra_Voronoi_der(ra2, rb2, rc2, rd2, rab, rac, rad, rbc,
		rbd, rcd, rab2, rac2, rad2, rbc2, rbd2, rcd2, cosine, sine,
		deriv, &vola, &volb, &volc, &vold, dvola, dvolb, dvolc, dvold, option);

		ballwvol[ia] += vola; ballwvol[ib] += volb;
		ballwvol[ic] += volc; ballwvol[id] += vold;

		if(option==1) {
			for(int iedge = 0; iedge < 6; iedge++) {
				int i1 = edge_list[iedge];
				edges[i1].dvol += coefaV*dvola[iedge]
				+coefbV*dvolb[iedge]+coefcV*dvolc[iedge]
				+coefdV*dvold[iedge];
			}
		}


/* 		====================================================================
	    	Weights on each vertex: fraction of solid angle
 		==================================================================== */

		vertices[ia].gamma -= (angle[0]+angle[1]+angle[2])/2 -0.250;
		vertices[ib].gamma -= (angle[0]+angle[3]+angle[4])/2 -0.250;
		vertices[ic].gamma -= (angle[1]+angle[3]+angle[5])/2 -0.250;
		vertices[id].gamma -= (angle[2]+angle[4]+angle[5])/2 -0.250;

/* 		====================================================================
		Weights on each edge: fraction of dihedral angle
 		==================================================================== */

		for(int iedge = 0; iedge < 6; iedge++) {
			int i1 = edge_list[iedge];
			if(edges[i1].gamma != 0) {
				edges[i1].gamma -= angle[iedge];
				edges[i1].sigma -= angle[iedge];
			}
		}

		if(option==1) {

/* 			====================================================================
			Derivative: take into account the derivatives of the edge weight in weighted 
			inclusion-exclusion formula
 			==================================================================== */

			for(int iedge = 0; iedge < 6; iedge++) {
				int i1 = edge_list[iedge];
				val1S = edges[i1].Surf;
				val1V = edges[i1].Vol;
				val1M1 = edges[i1].CoefM1 + edges[i1].CoefM2;
				val1G1 = edges[i1].CoefG1 + edges[i1].CoefG2;
				for(int ie = 0; ie < 6; ie++) {
					int je = edge_list[ie];
					edges[je].dsurf  += val1S*deriv[iedge][ie];
					edges[je].dvol   += val1V*deriv[iedge][ie];
					edges[je].dmean  += val1M1*deriv[iedge][ie];
					edges[je].dgauss += val1G1*deriv[iedge][ie];
				}
			}

/* 			====================================================================
			Derivative: take into account the derivatives of the vertex weight in weighted 
			inclusion-exclusion formula
 			==================================================================== */

			val1S = ra2*coefaS; val2S = rb2*coefbS;
			val3S = rc2*coefcS; val4S = rd2*coefdS;

			val1V = ra2*ra*coefaV/3; val2V = rb2*rb*coefbV/3;
			val3V = rc2*rc*coefcV/3; val4V = rd2*rd*coefdV/3;

			val1M1 = ra*coefaM; val2M1 = rb*coefbM;
			val3M1 = rc*coefcM; val4M1 = rd*coefdM;

			for(int ie = 0; ie < 6; ie++) {

				int je = edge_list[ie];
				d1 = deriv[0][ie]+deriv[1][ie]+deriv[2][ie];
				d2 = deriv[0][ie]+deriv[3][ie]+deriv[4][ie];
				d3 = deriv[1][ie]+deriv[3][ie]+deriv[5][ie];
				d4 = deriv[2][ie]+deriv[4][ie]+deriv[5][ie];

				val = val1S*d1 + val2S*d2 + val3S*d3 + val4S*d4;
				edges[je].dsurf -= val;

				val = val1V*d1 + val2V*d2 + val3V*d3 + val4V*d4;
				edges[je].dvol -= val;

				val = val1M1*d1 + val2M1*d2 + val3M1*d3 + val4M1*d4;
				edges[je].dmean -= val;

				val = coefaG*d1 + coefbG*d2 + coefcG*d3 + coefdG*d4;
				edges[je].dgauss -= val;

			}
		}

	}

/* ====================================================================
	Contribution of 3-balls (i.e. triangles of the alpha complex)
 ==================================================================== */

	for(int idx = 0; idx < nfaces; idx++) {

		coefval = faces[idx].gamma;
		if(coefval==0) continue;

		ia = faces[idx].Vertices[0];
		ib = faces[idx].Vertices[1];
		ic = faces[idx].Vertices[2];

		if(vertices[ia].status==0 || vertices[ib].status==0
		|| vertices[ic].status==0 ) continue;


		coefaS = vertices[ia].CoefS; coefaV = vertices[ia].CoefV; 
		coefaM = vertices[ia].CoefM; coefaG = vertices[ia].CoefG;
		coefbS = vertices[ib].CoefS; coefbV = vertices[ib].CoefV; 
		coefbM = vertices[ib].CoefM; coefbG = vertices[ib].CoefG;
		coefcS = vertices[ic].CoefS; coefcV = vertices[ic].CoefV; 
		coefcM = vertices[ic].CoefM; coefcG = vertices[ic].CoefG;

		e1 = faces[idx].Edges[0];
		e2 = faces[idx].Edges[1];
		e3 = faces[idx].Edges[2];

		ra = vertices[ia].Radius; ra2 = ra*ra;
		rb = vertices[ib].Radius; rb2 = rb*rb;
		rc = vertices[ic].Radius; rc2 = rc*rc;

		rab = edges[e1].Length; rab2=rab*rab;
		rac = edges[e2].Length; rac2=rac*rac;
		rbc = edges[e3].Length; rbc2=rbc*rbc;

		threesphere_dvol(ra, rb, rc, ra2, rb2, rc2, rab, rac, rbc, rab2, rac2, rbc2,
		angle, deriv2, &surfa, &surfb, &surfc, &vola, &volb, &volc, 
		dsurfa3, dsurfb3, dsurfc3, dvola3, dvolb3, dvolc3, option);

		gauss.threesphere_dgauss(ra, rb, rc, ra2, rb2, rc2, rab, rac, rbc,
		rab2, rac2, rbc2, &gaussa, &gaussb, &gaussc, dg, option);

		ballwsurf[ia] += coefval*surfa;
		ballwsurf[ib] += coefval*surfb;
		ballwsurf[ic] += coefval*surfc;

		ballwvol[ia] += coefval*vola;
		ballwvol[ib] += coefval*volb;
		ballwvol[ic] += coefval*volc;

		ballwgauss[ia] += coef_F*coefval*gaussa;
		ballwgauss[ib] += coef_F*coefval*gaussb;
		ballwgauss[ic] += coef_F*coefval*gaussc;

		edges[e1].sigma -= 2*coefval*angle[0];
		edges[e2].sigma -= 2*coefval*angle[1];
		edges[e3].sigma -= 2*coefval*angle[3];

		if(option==1) {
			edges[e1].dsurf += coefval*(coefaS*dsurfa3[0]+coefbS*dsurfb3[0]+coefcS*dsurfc3[0]);
			edges[e2].dsurf += coefval*(coefaS*dsurfa3[1]+coefbS*dsurfb3[1]+coefcS*dsurfc3[1]);
			edges[e3].dsurf += coefval*(coefaS*dsurfa3[2]+coefbS*dsurfb3[2]+coefcS*dsurfc3[2]);

			edges[e1].dvol += coefval*(coefaV*dvola3[0]+coefbV*dvolb3[0]+coefcV*dvolc3[0]);
			edges[e2].dvol += coefval*(coefaV*dvola3[1]+coefbV*dvolb3[1]+coefcV*dvolc3[1]);
			edges[e3].dvol += coefval*(coefaV*dvola3[2]+coefbV*dvolb3[2]+coefcV*dvolc3[2]);

			edges[e1].dmean += coefval*(coefaM*dsurfa3[0]/ra+coefbM*dsurfb3[0]/rb+
					coefcM*dsurfc3[0]/rc);
			edges[e2].dmean += coefval*(coefaM*dsurfa3[1]/ra+coefbM*dsurfb3[1]/rb+
					coefcM*dsurfc3[1]/rc);
			edges[e3].dmean += coefval*(coefaM*dsurfa3[2]/ra+coefbM*dsurfb3[2]/rb+
					coefcM*dsurfc3[2]/rc);

			edges[e1].dmean += 2*coefval*(edges[e1].CoefM2*deriv2[0][0]+
					edges[e2].CoefM2*deriv2[1][0]+edges[e3].CoefM2*deriv2[3][0]);
			edges[e2].dmean += 2*coefval*(edges[e1].CoefM2*deriv2[0][1]+
					edges[e2].CoefM2*deriv2[1][1]+edges[e3].CoefM2*deriv2[3][1]);
			edges[e3].dmean += 2*coefval*(edges[e1].CoefM2*deriv2[0][2]+
					edges[e2].CoefM2*deriv2[1][2]+edges[e3].CoefM2*deriv2[3][2]);

			edges[e1].dgauss += coefval*(coefaG*dsurfa3[0]/ra2+coefbG*dsurfb3[0]/rb2+
					coefcG*dsurfc3[0]/rc2);
			edges[e2].dgauss += coefval*(coefaG*dsurfa3[1]/ra2+coefbG*dsurfb3[1]/rb2+
					coefcG*dsurfc3[1]/rc2);
			edges[e3].dgauss += coefval*(coefaG*dsurfa3[2]/ra2+coefbG*dsurfb3[2]/rb2+
					coefcG*dsurfc3[2]/rc2);

			edges[e1].dgauss += 2*coefval*(edges[e1].CoefG2*deriv2[0][0]+
					edges[e2].CoefG2*deriv2[1][0]+edges[e3].CoefG2*deriv2[3][0]);
			edges[e2].dgauss += 2*coefval*(edges[e1].CoefG2*deriv2[0][1]+
					edges[e2].CoefG2*deriv2[1][1]+edges[e3].CoefG2*deriv2[3][1]);
			edges[e3].dgauss += 2*coefval*(edges[e1].CoefG2*deriv2[0][2]+
					edges[e2].CoefG2*deriv2[1][2]+edges[e3].CoefG2*deriv2[3][2]);

			edges[e1].dgauss += coef_F*coefval*(coefaG*dg[0][0]+coefbG*dg[1][0]+coefcG*dg[2][0]);
			edges[e2].dgauss += coef_F*coefval*(coefaG*dg[0][1]+coefbG*dg[1][1]+coefcG*dg[2][1]);
			edges[e3].dgauss += coef_F*coefval*(coefaG*dg[0][2]+coefbG*dg[1][2]+coefcG*dg[2][2]);

		}

	}

/* ====================================================================
	Now add contribution of two-sphere
 ==================================================================== */
					
	double eps = 1.e-10;
	for(int iedge = 0; iedge < nedges; iedge++) {

		coefval = edges[iedge].gamma;
		coefvalS = edges[iedge].sigma;


		if(std::abs(coefval) < eps) continue;

		ia = edges[iedge].Vertices[0];
		ib = edges[iedge].Vertices[1];

		if(vertices[ia].status==0 || vertices[ib].status==0) continue;

		coefaS = vertices[ia].CoefS; coefaV = vertices[ia].CoefV; 
		coefaM = vertices[ia].CoefM; coefaG = vertices[ia].CoefG;
		coefbS = vertices[ib].CoefS; coefbV = vertices[ib].CoefV; 
		coefbM = vertices[ib].CoefM; coefbG = vertices[ib].CoefG;

		ra = vertices[ia].Radius; ra2 = ra*ra;
		rb = vertices[ib].Radius; rb2 = rb*rb;

		rab = edges[iedge].Length; rab2 = rab*rab;

		twosphere_dinfo(ra, ra2, rb, rb2, rab, rab2, &surfa, &surfb,
		&vola, &volb, &r, &phi, &l, &dsurfa2, &dsurfb2, &dvola2, &dvolb2, 
		&dr, &dphi, &dl, option);

		ballwsurf[ia] -= coefval*surfa; 
		ballwsurf[ib] -= coefval*surfb; 
		ballwvol[ia]  -= coefval*vola; 
		ballwvol[ib]  -= coefval*volb; 

		val = coef_E*pi*coefvalS*r*phi;
		ballwmean[ia] -= val;
		ballwmean[ib] -= val;

		val = coef_E*pi*coefvalS*l;
		ballwgauss[ia] -= val;
		ballwgauss[ib] -= val;

		if(option==1) {
			edges[iedge].dsurf  -= coefval* (coefaS*dsurfa2 + coefbS*dsurfb2);
			edges[iedge].dvol   -= coefval* (coefaV*dvola2 + coefbV*dvolb2);
			edges[iedge].dmean  -= coefval* (coefaM*dsurfa2/ra + coefbM*dsurfb2/rb);
			edges[iedge].dmean  -= coef_E*coefvalS*pi*(r*dphi+phi*dr)*(coefaM+coefbM);
			edges[iedge].dgauss -= coefval* (coefaG*dsurfa2/ra2 + coefbG*dsurfb2/rb2);
			edges[iedge].dgauss -= coef_E*coefvalS*pi*dl*(coefaG+coefbG);
		}
	}

/* ====================================================================
	Now loop over vertices
 ==================================================================== */
					
	for(int i = 4; i < nvertices; i++) {

		coefval = vertices[i].gamma;
		if(vertices[i].info[0]==0) continue;
		if(vertices[i].info[7]==0) continue;
		if(coefval==0) continue;
		if(vertices[i].status==0) continue;

		ra = vertices[i].Radius; ra2 = ra*ra;
		surfa = 4*pi*ra*ra;
		vola  = surfa*ra/3;
		ballwsurf[i]  += coefval*surfa;
		ballwvol[i]   += coefval*vola;
		ballwmean[i]  += ballwsurf[i]/ra;
		ballwgauss[i] += ballwsurf[i]/ra2;
	}

/* ====================================================================
	Compute total surface, volume (weighted, and unweighted)
 ==================================================================== */
					
	for(int i = 4; i < nvertices; i++) {

		if(vertices[i].info[0]==0) continue;
		if(vertices[i].status==0) continue;

		coefaS = vertices[i].CoefS; coefaV = vertices[i].CoefV;
		coefaM = vertices[i].CoefM; coefaG = vertices[i].CoefG;

		*Surf          += ballwsurf[i];
		ballwsurf[i]    = ballwsurf[i]*coefaS;
		*WSurf         += ballwsurf[i];

		*Vol          += ballwvol[i];
		ballwvol[i]    = ballwvol[i]*coefaV;
		*WVol         += ballwvol[i];

		*Mean         += ballwmean[i];
		ballwmean[i]   = ballwmean[i]*coefaM;
		*WMean        += ballwmean[i];

		*Gauss        += ballwgauss[i];
		ballwgauss[i]  = ballwgauss[i]*coefaG;
		*WGauss       += ballwgauss[i];

	}

/* ====================================================================
	Shift as 4 first vertices are pseudo atoms
 ==================================================================== */

	int nballs = 0;
	for(int i = 0; i < nvertices; i++) if(vertices[i].status==1) nballs++;

	for(int i = 0; i < nballs; i++) {
		ballwsurf[i]   = ballwsurf[i+4];
		ballwvol[i]    = ballwvol[i+4];
		ballwmean[i]   = ballwmean[i+4];
		ballwgauss[i]  = ballwgauss[i+4];
	}

	if(option==0) return;

/* ====================================================================
	Convert derivatives wrt to distance to derivatives wrt to coordinates
 ==================================================================== */

	memset(dsurf_coord, 0, 3*nvertices*sizeof(double));
	memset(dvol_coord, 0, 3*nvertices*sizeof(double));
	memset(dmean_coord, 0, 3*nvertices*sizeof(double));
	memset(dgauss_coord, 0, 3*nvertices*sizeof(double));
					
	for(int iedge = 0; iedge < nedges; iedge++) {
		
		ia = edges[iedge].Vertices[0];
		ib = edges[iedge].Vertices[1];

		for(int i = 0; i < 3; i++) {
			u[i] = vertices[ia].Coordinates[i] - vertices[ib].Coordinates[i];
		}

		rab  = edges[iedge].Length;
		val1S  = edges[iedge].dsurf/rab;
		val1V  = edges[iedge].dvol/rab;
		val1M1 = edges[iedge].dmean/rab;
		val1G1 = edges[iedge].dgauss/rab;

		for(int j = 0; j < 3; j++) {
			dsurf_coord[3*ia+j]  += u[j]*val1S;
			dsurf_coord[3*ib+j]  -= u[j]*val1S;
			dvol_coord[3*ia+j]   += u[j]*val1V;
			dvol_coord[3*ib+j]   -= u[j]*val1V;
			dmean_coord[3*ia+j]  += u[j]*val1M1;
			dmean_coord[3*ib+j]  -= u[j]*val1M1;
			dgauss_coord[3*ia+j] += u[j]*val1G1;
			dgauss_coord[3*ib+j] -= u[j]*val1G1;
		}
	}

/* ====================================================================
	Shift as 4 first vertices are pseudo atoms
 ==================================================================== */

	for(int i = 0; i < 3*nballs; i++) {
		dsurf_coord[i]   = dsurf_coord[i+12];
		dvol_coord[i]    = dvol_coord[i+12];
		dmean_coord[i]   = dmean_coord[i+12];
		dgauss_coord[i]  = dgauss_coord[i+12];
	}

  }

/* ====================================================================
  distance2 computes the square of the distance between two
	sphere centers
 ==================================================================== */

  double VOLUMES::distance2(std::vector<Vertex>& vertices, int n1, int n2)
  {

	double x;
	double dist = 0;
	for(int i = 0; i < 3; i++) {
		x = vertices[n1].Coordinates[i] - vertices[n2].Coordinates[i];
		dist += x*x;
	}

	return dist;
  }

/* ====================================================================
  plane_dist computes the fraction of the distance between the center of sphere B
	and the Voronoi plane between this sphere and another sphere A.
 ==================================================================== */

  double VOLUMES::plane_dist(double ra2, double rb2, double rab2)
  {

	double lambda = 0.50 - (ra2-rb2)/(2*rab2);

	return lambda;
  }

/* ====================================================================
	twosphere_info calculates the volume of the intersection
	of two balls and the surface area  of the
	intersection of two corresponding spheres; it is only called when the
	intersection exists
	
	Input:
			rab	: distance between the centers of the 2 spheres
			rab2	: distance between the centers of the 2 spheres
				  (squared)
			ra,rb	: radii of sphere A and B, respectively
			ra2 and rb2 are the squared of the quantities
			above)
	Output
			surfa	: partial contribution of A to the total
				  surface of the intersection
			surfb	: partial contribution of B to the total
				  surface of the intersection
			vola	: partial contribution of A to the total
				  volume of the intersection
			volb	: partial contribution of B to the total
				  volume of the intersection
 ==================================================================== */

  void VOLUMES::twosphere_info(double ra, double ra2, double rb, double rb2,
		double rab, double rab2, double *surfa, double *surfb,
		double *vola, double *volb, double *r, double *phi, double *l)
  {

	double cosine, vala, valb, lambda, ha, hb;
	double Aab, sa, ca, sb, cb;

/* ====================================================================
	Get distance between center of sphere A and Voronoi plane
	between A and B
 ==================================================================== */

	lambda = plane_dist(ra2, rb2, rab2);
	valb = lambda*rab;
	vala = rab-valb;

/* ====================================================================
	Get height of the cap of sphere A occluded by sphere B
 ==================================================================== */

	ha = ra - vala;

/* ====================================================================
	same for sphere B ...
 ==================================================================== */

	hb = rb - valb;

/* ====================================================================
	Get surfaces of intersection
 ==================================================================== */

	*surfa = twopi*ra*ha;
	*surfb = twopi*rb*hb;

/* ====================================================================
	Now get volume
 ==================================================================== */

	Aab = pi*(ra2-vala*vala);

	sa = ra*(*surfa);
	ca = vala*Aab;

	*vola = (sa-ca)/3;

	sb = rb*(*surfb);
	cb = valb*Aab;

	*volb = (sb-cb)/3;

/* ====================================================================
	Get radius of the circle of intersection between the two spheres
 ==================================================================== */

	*r = std::sqrt(ra2 - vala*vala);

/* ====================================================================
	Get angle between normals of the sphere at a point on this circle
 ==================================================================== */

	cosine = (ra2+rb2-rab2)/(2.0*ra*rb);
	*phi = std::acos(cosine);

	*l = vala/ra + valb/rb;

  }

/* ====================================================================
	twosphere_vol calculates the volume of the intersection
	of two balls and the surface area  of the
	intersection of two corresponding spheres; it is only called when the
	intersection exists

	It also computes the derivatives of the surface area
	and volume with respect to the distance between the two centers
	
	Input:
			rab	: distance between the centers of the 2 spheres
			rab2	: distance between the centers of the 2 spheres
				  (squared)
			ra,rb	: radii of sphere A and B, respectively
			ra2 and rb2 are the squared of the quantities
			above)
			option  : flag; if 1 computes derivatives
	Output
			surfa	: partial contribution of A to the total
				  surface of the intersection
			surfb	: partial contribution of B to the total
				  surface of the intersection
			vola	: partial contribution of A to the total
				  volume of the intersection
			volb	: partial contribution of B to the total
				  volume of the intersection
			dsurfa	: derivative of surfa with respect to rab
			dsurfb	: derivative of surfb with respect to rab
			dvola	: derivative of vola with respect to rab
			dvolb	: derivative of volb with respect to rab

 ==================================================================== */

  void VOLUMES::twosphere_dinfo(double ra, double ra2, double rb, double rb2,
		double rab, double rab2, double *surfa, double *surfb,
		double *vola, double *volb, double *r, double *phi, double *l,
		double *dsurfa, double *dsurfb, double *dvola, double *dvolb, 
		double *dr, double *dphi, double *dl, int option)
  {

	double cosine, vala, valb, lambda, ha, hb;
	double Aab, sa, ca, sb, cb;
	double	dera, derb;

/* ====================================================================
	Get distance between center of sphere A and Voronoi plane
	between A and B
 ==================================================================== */

	lambda = plane_dist(ra2, rb2, rab2);
	valb = lambda*rab;
	vala = rab-valb;

/* ====================================================================
	Get height of the cap of sphere A occluded by sphere B
 ==================================================================== */

	ha = ra - vala;

/* ====================================================================
	same for sphere B ...
 ==================================================================== */

	hb = rb - valb;

/* ====================================================================
	Get surfaces of intersection
 ==================================================================== */

	*surfa = twopi*ra*ha;
	*surfb = twopi*rb*hb;

/* ====================================================================
	Now get volume
 ==================================================================== */

	Aab = pi*(ra2-vala*vala);

	sa = ra*(*surfa);
	ca = vala*Aab;

	*vola = (sa-ca)/3;

	sb = rb*(*surfb);
	cb = valb*Aab;

	*volb = (sb-cb)/3;

/* ====================================================================
	Get radius of the circle of intersection between the two spheres
 ==================================================================== */

	*r = std::sqrt(ra2 - vala*vala);

/* ====================================================================
	Get angle between normals of the sphere at a point on this circle
 ==================================================================== */

	cosine = (ra2+rb2-rab2)/(2.0*ra*rb);
	*phi = std::acos(cosine);
	*l = vala/ra + valb/rb;

	if(option==0) return;

	dera = - lambda;
	derb = lambda - 1;

	*dsurfa = twopi*ra*dera;
	*dsurfb = twopi*rb*derb;

	*dvola = -Aab*lambda;
	*dvolb = -(*dvola) - Aab;

	*dr   = -vala*lambda/(*r);
	*dphi = rab/(ra*rb*std::sqrt(1-cosine*cosine));
	*dl   = lambda/ra + (1.0-lambda)/rb;

  }

/* ====================================================================
	threesphere_dvol

	This procedure computes the surface area and volume of the intersection 
	of three spheres; it is only called when the intersection exists

	It also computes the derivatives of the surface areas and volumes with respect
	to the three distances rAB, rAC and rBC

	Input:
			ra,rb,rc  : radii of sphere A, B and C, respectively
			ra2,rb2,rc2: radii squared
			rab,rab2: distance between the centers of sphere A and B
			rac,rac2: distance between the centers of sphere A and C
			rbc,rbc2: distance between the centers of sphere B and C
			option:   flag; if set to 1, compute derivatives
	Output
			surfa,surfb,surfc : contribution of A, B and C to
			the total surface of the intersection of A,B,C
			dsurfa		  : derivatives of surfa wrt rAB, rAC and rBC
			dsurfb		  : derivatives of surfb wrt rAB, rAC and rBC
			dsurfc		  : derivatives of surfc wrt rAB, rAC and rBC
			dvola		  : derivatives of vola wrt rAB, rAC and rBC
			dvolb		  : derivatives of volb wrt rAB, rAC and rBC
			dvolc		  : derivatives of volc wrt rAB, rAC and rBC

 ==================================================================== */

   void VOLUMES::threesphere_dvol(double ra, double rb,double rc, double ra2,
	double rb2, double rc2, double rab, double rac, double rbc,
	double rab2, double rac2, double rbc2, double *angle, double deriv[6][3],
	double *surfa, double *surfb, double *surfc, double *vola, double *volb, double *volc,
	double *dsurfa, double *dsurfb, double *dsurfc,
	double *dvola, double *dvolb, double *dvolc, int option)
  {

	double	a1, a2, a3, s2, c1, c2;
	double	seg_ang_ab, seg_ang_ac, seg_ang_bc;
	double	ang_dih_ap, ang_dih_bp, ang_dih_cp;
	double	val1, val2, val3, l1, l2, l3;
	double	val1b, val2b, val3b;
	double	ang_abc, ang_acb, ang_bca;
	double	cos_abc, cos_acb, cos_bca;
	double	sin_abc, sin_acb, sin_bca;
	double	s_abc, s_acb, s_bca;
	double	rho_ab2, rho_ac2, rho_bc2;
	double  drho_ab2, drho_ac2, drho_bc2;
	double	val_abc, val_acb, val_bca;
	double  val2_abc, val2_acb, val2_bca;
	double	der_val1b, der_val1, der_val2b, der_val2, der_val3b, der_val3;

	double cosine[6], sine[6];

	l1 = plane_dist(ra2, rb2, rab2);
	l2 = plane_dist(ra2, rc2, rac2);
	l3 = plane_dist(rb2, rc2, rbc2);

	val1 = l1*rab; val2 = l2*rac; val3 = l3*rbc;
	val1b = rab - val1; val2b = rac - val2; val3b = rbc - val3;

/* ====================================================================
	We consider the tetrahedron (A,B,C,P) where P is the
	point of intersection of the three spheres such that (A,B,C,P) is ccw.

	The edge lengths in this tetrahedron are: rab, rac, rAP=ra, rbc, rBP=rb, rCP=rc
 ==================================================================== */

	tetrageom.tetra_dihed_der3(rab2, rac2, ra2, rbc2, rb2, rc2, angle, 
		cosine, sine, deriv, option);

/* ====================================================================
	The seg_ang_ are the dihedral angles around the three edges AB, AC and BC
 ==================================================================== */

	seg_ang_ab = angle[0];
	seg_ang_ac = angle[1];
	seg_ang_bc = angle[3];

/* ====================================================================
	The ang_dih_ are the dihedral angles around the three edges AP, BP and CP
 ==================================================================== */

	ang_dih_ap = angle[2];
	ang_dih_bp = angle[4];
	ang_dih_cp = angle[5];

	a1 = ra*(1-2*ang_dih_ap);
	a2 = 2*seg_ang_ab*val1b;
	a3 = 2*seg_ang_ac*val2b;

	*surfa = twopi*ra*(a1 - a2 - a3);

	a1 = rb*(1-2*ang_dih_bp);
	a2 = 2*seg_ang_ab*val1;
	a3 = 2*seg_ang_bc*val3b;

	*surfb = twopi*rb*(a1 - a2 - a3);

	a1 = rc*(1-2*ang_dih_cp);
	a2 = 2*seg_ang_ac*val2;
	a3 = 2*seg_ang_bc*val3;

	*surfc = twopi*rc*(a1 - a2 - a3);

/* ====================================================================
	compute volumes of the three caps
 ==================================================================== */

	ang_abc = twopi*seg_ang_ab;
	ang_acb = twopi*seg_ang_ac;
	ang_bca = twopi*seg_ang_bc;

	cos_abc = cosine[0];
	sin_abc = sine[0];
	cos_acb = cosine[1];
	sin_acb = sine[1];
	cos_bca = cosine[3];
	sin_bca = sine[3];

	rho_ab2 = ra2 - val1b*val1b;
	rho_ac2 = ra2 - val2b*val2b;
	rho_bc2 = rb2 - val3b*val3b;

	val_abc = ang_abc - sin_abc*cos_abc; s_abc = rho_ab2*val_abc;
	val_acb = ang_acb - sin_acb*cos_acb; s_acb = rho_ac2*val_acb;
	val_bca = ang_bca - sin_bca*cos_bca; s_bca = rho_bc2*val_bca;

	s2 = ra*(*surfa);
	c1 = val1b*s_abc;
	c2 = val2b*s_acb;

	*vola = (s2 - c1 - c2)/3;

	s2 = rb*(*surfb);
	c1 = val1*s_abc;
	c2 = val3b*s_bca;

	*volb = (s2 - c1 - c2)/3;

	s2 = rc*(*surfc);
	c1 = val2*s_acb;
	c2 = val3*s_bca;

	*volc = (s2 - c1 - c2)/3;

	if(option==0) return;

	der_val1b = l1; der_val1  = 1-l1;
	der_val2b = l2; der_val2  = 1-l2;
	der_val3b = l3; der_val3  = 1-l3;

	dsurfa[0] = -2*ra*(
		    twopi*seg_ang_ab*der_val1b +
		    (ra*deriv[2][0] +
			val1b*deriv[0][0] +val2b*deriv[1][0]));
	dsurfa[1] = -2*ra*(
		    twopi*seg_ang_ac*der_val2b +
		    (ra*deriv[2][1] +
			val1b*deriv[0][1] +val2b*deriv[1][1]));
	dsurfa[2] = -2*ra*( ra*deriv[2][2] +
			val1b*deriv[0][2]+val2b*deriv[1][2]);

	dsurfb[0] = -2*rb*(
			twopi*seg_ang_ab*der_val1
			+(rb*deriv[4][0]+
			val1*deriv[0][0]+val3b*deriv[3][0]));
	dsurfb[1] = -2*rb*(rb*deriv[4][1]+
			val1*deriv[0][1]+val3b*deriv[3][1]);
	dsurfb[2] = -2*rb*(
			twopi*seg_ang_bc*der_val3b
			+(rb*deriv[4][2]+
			val1*deriv[0][2]+val3b*deriv[3][2]));

	dsurfc[0] = -2*rc*(rc*deriv[5][0]+
			val2*deriv[1][0]+val3*deriv[3][0]);
	dsurfc[1] = -2*rc*(
			twopi*seg_ang_ac*der_val2
			+(rc*deriv[5][1]+
			val2*deriv[1][1]+val3*deriv[3][1]));
	dsurfc[2] = -2*rc*(
			twopi*seg_ang_bc*der_val3
			+(rc*deriv[5][2]+
			val2*deriv[1][2]+val3*deriv[3][2]));

	drho_ab2 = -2*der_val1b*val1b;
	drho_ac2 = -2*der_val2b*val2b;
	drho_bc2 = -2*der_val3b*val3b;

	val2_abc = rho_ab2*(1 - cos_abc*cos_abc + sin_abc*sin_abc);
	val2_acb = rho_ac2*(1 - cos_acb*cos_acb + sin_acb*sin_acb);
	val2_bca = rho_bc2*(1 - cos_bca*cos_bca + sin_bca*sin_bca);

	dvola[0] = ra*dsurfa[0] - der_val1b*s_abc - 
		(val1b*deriv[0][0]*val2_abc + val2b*deriv[1][0]*val2_acb)
		- val1b*drho_ab2*val_abc;
	dvola[0] = dvola[0]/3;
	dvola[1] = ra*dsurfa[1] - der_val2b*s_acb - 
		(val1b*deriv[0][1]*val2_abc + val2b*deriv[1][1]*val2_acb)
		- val2b*drho_ac2*val_acb;
	dvola[1] = dvola[1]/3;
	dvola[2] = ra*dsurfa[2] - 
		(val1b*deriv[0][2]*val2_abc + val2b*deriv[1][2]*val2_acb);
	dvola[2] = dvola[2]/3;

	dvolb[0] = rb*dsurfb[0] - der_val1*s_abc - 
		(val1*deriv[0][0]*val2_abc + val3b*deriv[3][0]*val2_bca)
		- val1*drho_ab2*val_abc;
	dvolb[0] = dvolb[0]/3;
	dvolb[1] = rb*dsurfb[1] - 
		(val1*deriv[0][1]*val2_abc + val3b*deriv[3][1]*val2_bca);
	dvolb[1] = dvolb[1]/3;
	dvolb[2] = rb*dsurfb[2] - der_val3b*s_bca - 
		(val1*deriv[0][2]*val2_abc + val3b*deriv[3][2]*val2_bca)
		- val3b*drho_bc2*val_bca;
	dvolb[2] = dvolb[2]/3;

	dvolc[0] = rc*dsurfc[0] - 
		(val2*deriv[1][0]*val2_acb + val3*deriv[3][0]*val2_bca);
	dvolc[0] = dvolc[0]/3;
	dvolc[1] = rc*dsurfc[1] - der_val2*s_acb - 
		(val2*deriv[1][1]*val2_acb + val3*deriv[3][1]*val2_bca)
		- val2*drho_ac2*val_acb;
	dvolc[1] = dvolc[1]/3;
	dvolc[2] = rc*dsurfc[2] - der_val3*s_bca - 
		(val2*deriv[1][2]*val2_acb + val3*deriv[3][2]*val2_bca)
		- val3*drho_bc2*val_bca;
	dvolc[2] = dvolc[2]/3;

  }

#endif
