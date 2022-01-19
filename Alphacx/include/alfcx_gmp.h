/* =======================================================================================
	alfx_gmp

        Performs all predicates for alfcx using gmp
                                                                        
 ======================================================================================= */

#ifndef ALFCX_GMP_H
#define ALFCX_GMP_H

/* Includes :								  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* =======================================================================================
	the class
 ======================================================================================= */

  class ALFCX_GMP {

	public:

		void tetra_radius_gmp(double *a, double *b, double *c, double *d, double ra,
		double rb, double rc, double rd, int *test, double alpha);

		void vertex_attach_gmp(double *a, double *b, double ra, double rb, int *testa, int *testb);

		void edge_attach_gmp(double *a, double *b, double *c, double ra,
		double rb, double rc, int *test, int *memory);

		void edge_radius_gmp(double *a, double *b, double ra, double rb,
		int *test, double alpha, int *memory);

		void triangle_attach_gmp(double *a, double *b, double *c, double *d,
        	double ra, double rb, double rc, double rd, int *test, int *memory);

		void triangle_radius_gmp(double *a, double *b, double *c, double ra,
		double rb, double rc, int *test, double alpha, int *memory);

		void set_alf_gmp();
		void clear_alf_gmp();

	private:

		void set_edge(double *a, double *b, double ra, double rb);

		void set_triangle(double *a, double *b, double *c, double ra, double rb, double rc);

		void real_to_gmp(double *coord, int idx, mpz_t val);

		void build_weight(mpz_t ax, mpz_t ay, mpz_t az, mpz_t r, mpz_t w);

		void scalar_to_gmp(double coord, mpz_t val);

		mpz_t temp1, temp2, temp3;

		mpz_t ra2,rb2,dist2,dtest, num, den;
		mpz_t r_11, r_22, r_33, r_14, r_313, r_212,diff, det0, det1, det2, det3, det4;
		mpz_t Dabc, Dabd, Dacd, Dbcd, Dabcd;
		mpz_t wa,wb,wc,wd;
		
		mpz_t ra_mp,rb_mp, rc_mp, rd_mp;
		mpz_t alp;
		
		mpz_t res[4][5], res2_c[4][5];
		mpz_t a_mp[5], b_mp[5], c_mp[5], d_mp[5];
		mpz_t Tab[4], Sab[4], Dab[5];
		mpz_t Sac[4], Sad[4], Sbc[4], Sbd[4], Scd[4];
		mpz_t Sa[4], Sb[4], Sd[4];
		mpz_t Sam1[4], Sbm1[4], Scm1[4], Sdm1[4];
		mpz_t Deter[4];
		mpz_t Tc[4],Sc[4];
		mpz_t Mab[4][5], Mac[4][5], Mbc[4][5], S[4][5], T[4][5];

		double scale = 1.e8;
    };

/* =======================================================================================
	scalar_to_gmp: converts one coordinate of one point into mpz_t type
 ======================================================================================= */

void ALFCX_GMP::scalar_to_gmp(double coord, mpz_t val)
{
        double x;
        int ivalue;

        mpz_set_d(temp3, scale);

        ivalue= (int) coord;
        mpz_set_si(temp1,ivalue);
        mpz_mul(temp1,temp1,temp3);
        x = (coord-ivalue)*(scale);
        ivalue = (int) round(x);
        mpz_set_si(temp2,ivalue);
        mpz_add(val,temp1,temp2);
}

/* =======================================================================================
	real_to_gmp: converts one coordinate of one point into mpz_t type
 ======================================================================================= */

void ALFCX_GMP::real_to_gmp(double *coord, int i, mpz_t val)
{
        double x;
        int ivalue;

        mpz_set_d(temp3, scale);

        ivalue= (int) coord[i];
        mpz_set_si(temp1,ivalue);
        mpz_mul(temp1,temp1,temp3);
        x = (coord[i]-ivalue)*(scale);
        ivalue = (int) round(x);
        mpz_set_si(temp2,ivalue);
        mpz_add(val,temp1,temp2);
}

/* =======================================================================================
	build_weight: builds the weight of a point:  w = x**2 + y**2 + z**2 - ra**2
 ======================================================================================= */

void ALFCX_GMP::build_weight(mpz_t ax, mpz_t ay, mpz_t az, mpz_t r, mpz_t w)
{
	mpz_mul(temp1,r,r);
	mpz_mul(temp2,ax,ax), mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,ay,ay), mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,az,az), mpz_add(w,temp2,temp1);
}

/* =======================================================================================
 	set_triangle: sets all common arrays for gmp calculation over edges
 ======================================================================================= */

void ALFCX_GMP::set_triangle(double *a, double *b, double *c, double ra, double rb, double rc)
{
	int i,j;

//       0. define coordinates

	for (i=0; i<3; i++)
	{
		real_to_gmp(a, i, a_mp[i+1]);
		real_to_gmp(b, i, b_mp[i+1]);
		real_to_gmp(c, i, c_mp[i+1]);
	}
	
	scalar_to_gmp(ra, ra_mp);
	scalar_to_gmp(rb, rb_mp);
	scalar_to_gmp(rc, rc_mp);

	build_weight(a_mp[1],a_mp[2],a_mp[3],ra_mp,a_mp[4]);
	build_weight(b_mp[1],b_mp[2],b_mp[3],rb_mp,b_mp[4]);
	build_weight(c_mp[1],c_mp[2],c_mp[3],rc_mp,c_mp[4]);

/*
       1. Computes all Minors Mab(i,j)= M(a,b,i,j)   = Det | a(i)  a(j) |
                                                           | b(i)  b(j) |
*/
	for (i=1;  i<4; i++)
	{
		for (j=i+1; j<5 ; j++)
		{
			mpz_mul(temp1, a_mp[j], b_mp[i]); 
			mpz_mul(temp2, a_mp[i], b_mp[j]);
			mpz_sub(Mab[i][j], temp2, temp1);
			mpz_mul(temp1, a_mp[j], c_mp[i]); 
			mpz_mul(temp2, a_mp[i], c_mp[j]);
			mpz_sub(Mac[i][j], temp2, temp1);
			mpz_mul(temp1, b_mp[j], c_mp[i]); 
			mpz_mul(temp2, b_mp[i], c_mp[j]);
			mpz_sub(Mbc[i][j], temp2, temp1);
		}
	}

/*
       Now compute all Minors
               S(i,j) = M(a,b,c,i,j,0)    = Det | a(i) a(j) 1 |
                                                | b(i) b(j) 1 |
                                                | c(i) c(j) 1 |

       a,b,c are the 3 vertices of the triangle, i and j correspond
       to two of the coordinates of the vertices

       for all i in [1,3] and all j in [i+1,4]
*/

	for (i=1;  i<4; i++)
	{
		for (j=i+1; j<5 ; j++)
		{
			mpz_sub(temp1,Mbc[i][j],Mac[i][j]);
			mpz_add(S[i][j],temp1,Mab[i][j]);
		}
	}

/*
       Now compute all Minors
               T(i,j) = M(a,b,c,i,j,4)    = Det | a(i) a(j) a(4) |
                                                | b(i) b(j) b(4) |
                                                | c(i) c(j) c(4) |

       for all i in [1,2] and all j in [i+1,3]
*/

	for (i=1;  i<3; i++)
	{
		for (j=i+1; j<4 ; j++)
		{
			mpz_mul(temp1,a_mp[4],Mbc[i][j]);
			mpz_mul(temp2,b_mp[4],Mac[i][j]);
			mpz_sub(temp1,temp1,temp2);
			mpz_mul(temp2,c_mp[4],Mab[i][j]);
			mpz_add(T[i][j],temp1,temp2);
		}
	}
/*
       Finally,  need Dabc = M(a,b,c,1,2,3) Det | a(1) a(2) a(3) |
                                                | b(2) b(2) b(3) |
                                                | c(3) c(2) c(3) |
*/
	mpz_mul(temp1,a_mp[1],Mbc[2][3]); mpz_mul(temp2,b_mp[1],Mac[2][3]); 
	mpz_sub(temp1,temp1,temp2);
	mpz_mul(temp2,c_mp[1],Mab[2][3]); mpz_add(Dabc,temp1,temp2);

}

/* =======================================================================================
 	set_edge: sets all common arrays for gmp calculation over edges
 ======================================================================================= */

void ALFCX_GMP::set_edge(double *a, double *b, double ra, double rb)
{
	int i,j,k;

//       0. define coordinates


	for (i=0; i<3; i++)
	{
		real_to_gmp(a, i, a_mp[i+1]);
		real_to_gmp(b, i, b_mp[i+1]);
	}
	
	scalar_to_gmp(ra,ra_mp);
	scalar_to_gmp(rb,rb_mp);

	build_weight(a_mp[1],a_mp[2],a_mp[3],ra_mp,a_mp[4]);
	build_weight(b_mp[1],b_mp[2],b_mp[3],rb_mp,b_mp[4]);

/*
       1. Compute all Minors Dab(i) = M(a,b,i,0) = Det | a(i) 1 |
                                                       | b(i) 1 |
*/
	for (i=1; i<5; i++)
	{
		mpz_sub(Dab[i],a_mp[i],b_mp[i]);
	}

/*
       2. Computes all Minors Sab(i,j)= M(a,b,i,j)   = Det | a(i)  a(j) |
                                                           | b(i)  b(j) |
*/
	for (i=1;  i<3; i++)
	{
		for (j=i+1; j<4 ; j++)
		{
			k=i+j-2;
			mpz_mul(temp1,a_mp[j],b_mp[i]); 
			mpz_mul(temp2,a_mp[i],b_mp[j]);
			mpz_sub(Sab[k],temp2,temp1);
		}
	}

/*
      3. Computes all Minors Tab(i)= M(a,b,i,4)   = Det | a(i)  a(4) |
                                                         | b(i)  b(4) |
*/
	for (i=1; i<4; i++)
	{
		mpz_mul(temp1,a_mp[i],b_mp[4]); 
		mpz_mul(temp2,a_mp[4],b_mp[i]);
		mpz_sub(Tab[i],temp1,temp2);
	}
}

/* =======================================================================================
	tetra_radius_gmp : computes the radius R of the circumsphere containing
 	a tetrahedron [A,B,C,D], as well as check if any fourth point L 
 	(in A, B, C, D) of the tetrahedron is "hidden" by its opposite 
 	face [I,J,K], i.e. is interior to the circumsphere of [I,J,K]

	The array "coord" contains the coordinates of A,B,C,D in that order
	The array "radius" contains the radii of A,B,C and D
 
 	Since we are only interested at how R compares to Alpha, we don't
 	output R, rather the result of the comparison
 
 	Computation extends to all four faces of the tetrahedron, as
 	well as to all six edges of the tetrahedron
 
	This procedure works with Multiple Precision Integer Arithmetics  (MPIA)
 	The package GMP is used for MPIA (with a C wrapper)
 ======================================================================================= */

void ALFCX_GMP::tetra_radius_gmp( double *a, double *b, double *c, double *d, double ra,
		double rb, double rc, double rd, int *testr, double alpha)
{
	int i,j,k,coef;
	int ivalue;
	double value;

/*	Transfer data in multiple precision */

	for (i=0; i<3; i++)
	{
		real_to_gmp(a, i, a_mp[i+1]);
		real_to_gmp(b, i, b_mp[i+1]);
		real_to_gmp(c, i, c_mp[i+1]);
		real_to_gmp(d, i, d_mp[i+1]);
	}
	
	scalar_to_gmp(ra, ra_mp);
	scalar_to_gmp(rb, rb_mp);
	scalar_to_gmp(rc, rc_mp);
	scalar_to_gmp(rd, rd_mp);

	build_weight(a_mp[1], a_mp[2], a_mp[3], ra_mp, wa);
	build_weight(b_mp[1], b_mp[2], b_mp[3], rb_mp, wb);
	build_weight(c_mp[1], c_mp[2], c_mp[3], rc_mp, wc);
	build_weight(d_mp[1], d_mp[2], d_mp[3], rd_mp, wd);

	value = alpha*scale; ivalue = (int) floor(value); 
	mpz_set_si(alp,ivalue);

/*	1. Computes all Minors Smn(i+j-2)= M(m,n,i,j) = Det | m(i)  m(j) |
						            | n(i)  n(j) |
	for all i in [1,2] and all j in [i+1,3]                         */

	for (i=1;  i<3; i++)
	{
		for (j=i+1; j<4 ; j++)
		{
			k=i+j-2;
			mpz_mul(temp1, a_mp[j], b_mp[i]); 
			mpz_mul(temp2, a_mp[i], b_mp[j]);
			mpz_sub(Sab[k], temp2, temp1);
			mpz_mul(temp1, a_mp[j], c_mp[i]); 
			mpz_mul(temp2, a_mp[i], c_mp[j]);
			mpz_sub(Sac[k], temp2, temp1);
			mpz_mul(temp1, a_mp[j], d_mp[i]); 
			mpz_mul(temp2, a_mp[i], d_mp[j]);
			mpz_sub(Sad[k], temp2, temp1);
			mpz_mul(temp1, b_mp[j], c_mp[i]); 
			mpz_mul(temp2, b_mp[i], c_mp[j]);
			mpz_sub(Sbc[k], temp2, temp1);
			mpz_mul(temp1, b_mp[j], d_mp[i]); 
			mpz_mul(temp2, b_mp[i], d_mp[j]);
			mpz_sub(Sbd[k], temp2, temp1);
			mpz_mul(temp1, c_mp[j], d_mp[i]); 
			mpz_mul(temp2, c_mp[i], d_mp[j]);
			mpz_sub(Scd[k], temp2, temp1);
		}
	}

/*	Now compute all Minors 
		Sq(i+j-2) = M(m,n,p,i,j,0) = Det | m(i) m(j) 1 |
		       			         | n(i) n(j) 1 |
						 | p(i) p(j) 1 |

	and all Minors
		Det(i+j-2) = M(m,n,p,q,i,j,4,0) = Det | m(i) m(j) m(4) 1 |
						      | n(i) n(j) n(4) 1 |
						      | p(i) p(j) p(4) 1 |
						      | q(i) q(j) q(4) 1 |

	m,n,p,q are the four vertices of the tetrahedron, i and j correspond
	to two of the coordinates of the vertices, and m(4) refers to the
	"weight" of vertices m                                           */
 
	for (i=1; i<4; i++)
	{
		mpz_sub(temp1, Scd[i], Sbd[i]); mpz_add(Sa[i], temp1, Sbc[i]);
		mpz_mul(temp2, Sa[i], wa);
		mpz_sub(temp1, Scd[i], Sad[i]); mpz_add(Sb[i], temp1, Sac[i]);
		mpz_mul(temp3, Sb[i], wb); mpz_sub(temp2, temp2, temp3);
		mpz_sub(temp1, Sbd[i], Sad[i]); mpz_add(Sc[i], temp1, Sab[i]);
		mpz_mul(temp3, Sc[i], wc); mpz_add(temp2, temp2, temp3);
		mpz_sub(temp1, Sbc[i], Sac[i]); mpz_add(Sd[i], temp1, Sab[i]);
		mpz_mul(temp3, Sd[i], wd); mpz_sub(Deter[i], temp2, temp3);
		mpz_neg(Sam1[i], Sa[i]); mpz_neg(Sbm1[i], Sb[i]);
		mpz_neg(Scm1[i], Sc[i]); mpz_neg(Sdm1[i], Sd[i]);
	}
 
/*
	Now compute the determinant needed to compute the radius of the
	circumsphere of the tetrahedron :

		Det1 = Minor(a,b,c,d,4,2,3,0)
		Det2 = Minor(a,b,c,d,1,3,4,0)
		Det3 = Minor(a,b,c,d,1,2,4,0)
		Det4 = Minor(a,b,c,d,1,2,3,0)
									*/

	mpz_set(det1, Deter[3]);
	mpz_set(det2, Deter[2]);
	mpz_set(det3, Deter[1]);

	mpz_mul(temp1, a_mp[1], Sa[3]);mpz_mul(temp2, b_mp[1], Sb[3]);
	mpz_sub(temp3, temp1, temp2);
	mpz_mul(temp1, c_mp[1], Sc[3]);mpz_mul(temp2, d_mp[1], Sd[3]);
	mpz_sub(temp1, temp1, temp2);
	mpz_add(det4, temp1, temp3);

/*
	Now compute all minors:
		Dmnp = Minor(m, n, p, 1, 2, 3) = Det | m(1) m(2) m(3) |
						| n(1) n(2) n(3) |
						| p(1) p(2) p(3) |
									*/

	mpz_mul(temp1, a_mp[1], Sbc[3]); mpz_mul(temp2, b_mp[1], Sac[3]);
	mpz_sub(temp3, temp1, temp2);
	mpz_mul(temp1, c_mp[1], Sab[3]);mpz_add(Dabc, temp3, temp1);

	mpz_mul(temp1, a_mp[1], Sbd[3]); mpz_mul(temp2, b_mp[1], Sad[3]);
	mpz_sub(temp3, temp1, temp2);
	mpz_mul(temp1, d_mp[1], Sab[3]);mpz_add(Dabd, temp3, temp1);

	mpz_mul(temp1, a_mp[1], Scd[3]); mpz_mul(temp2, c_mp[1], Sad[3]);
	mpz_sub(temp3, temp1, temp2);
	mpz_mul(temp1, d_mp[1], Sac[3]);mpz_add(Dacd, temp3, temp1);

	mpz_mul(temp1, b_mp[1], Scd[3]); mpz_mul(temp2, c_mp[1], Sbd[3]);
	mpz_sub(temp3, temp1, temp2);
	mpz_mul(temp1, d_mp[1], Sbc[3]);mpz_add(Dbcd, temp3, temp1);

/*
	We also need :
		Det = Det | m(1) m(2) m(3) m(4) |
			  | n(1) n(2) n(3) n(4) |
			  | p(1) p(2) p(3) p(4) |
			  | q(1) q(2) q(3) q(4) |
								*/

	mpz_mul(temp1, wa, Dbcd); mpz_mul(temp2, wb, Dacd);
	mpz_sub(temp3, temp2, temp1);
	mpz_mul(temp1, wc, Dabd); mpz_mul(temp2, wd, Dabc);
	mpz_sub(temp1, temp2, temp1); mpz_add(Dabcd, temp3, temp1);

/*
	The radius of the circumsphere of the weighted tetrahedron is then:
	r_t = (Det1*Det1 + Det2*Det2 + Det3*Det3 + 4*Det4*Dabcd)/(4*Det4*Det4)
								*/

	mpz_mul(temp1, det4, det4); coef=4; mpz_mul_si(den, temp1, coef);

	mpz_mul(temp1, det1, det1); mpz_mul(temp2, det2, det2);
	mpz_add(temp1, temp1, temp2); mpz_mul(temp2, det3, det3);
	mpz_add(temp1, temp1, temp2); mpz_mul(temp2, det4, Dabcd);
	mpz_mul_si(temp2, temp2, coef); mpz_add(num, temp1, temp2);
	
	mpz_mul(temp1, den, alp); mpz_sub(temp2, num, temp1);
	
/* 
 	If tetrahedron is part of the alpha shape, then the 4 triangles,
 	the 6 edges and the four vertices are also part of the alpha
 	complex
 								*/
	if(!(mpz_sgn(temp2) > 0)) 
	{
		(*testr)=1;
	}
	else
	{
		(*testr)=0;
	}
}

/* =======================================================================================

 vertex_attach_gmp: checks if a vertex is attached to another vertex 
	Input:
		a, b	: coordinates of the two points
		ra,rb	: radii of the two points
	Output:
		testa	: flag equal to 1 if a is attached to b
		testb	: flag equal to 1 if b is attached to a
 ======================================================================================= */

void ALFCX_GMP::vertex_attach_gmp(double *a, double *b, double ra, double rb, 
	int *testa, int *testb)

{
	int i;

	(*testa = 0);
	(*testb = 0);

	for (i=0; i<4; i++)
	{
		real_to_gmp(a, i, temp1);
		real_to_gmp(b, i, temp2);
		mpz_sub(Dab[i], temp1, temp2);
	}

	scalar_to_gmp(ra, ra_mp);
	scalar_to_gmp(rb, rb_mp);

	mpz_mul(temp1, Dab[0], Dab[0]);
	mpz_mul(temp2, Dab[1], Dab[1]);
	mpz_mul(temp3, Dab[2], Dab[2]);
	mpz_add(temp1, temp1, temp2); mpz_add(dist2, temp1, temp3);

	mpz_mul(ra2, ra_mp, ra_mp); mpz_mul(rb2, rb_mp, rb_mp);

	mpz_add(dtest, dist2, ra2);
	mpz_sub(dtest, dtest, rb2);
	if(mpz_sgn(dtest) < 0) (*testa = 1);

	mpz_sub(dtest, dist2, ra2);
	mpz_add(dtest, dtest, rb2);
	if(mpz_sgn(dtest) < 0) (*testb = 1);

}

/* =======================================================================================
   edge_attach_gmp: checks if an edge of the regular triangulation is "attached"
	to another vertex (i.e.  if the vertex belongs to the smallest circumsphere
	of the edge). 
 ======================================================================================= */

void ALFCX_GMP::edge_attach_gmp(double *a, double *b, double *c, double ra, 
			double rb, double rc, int *testa, int *memory)

{
	int i,j,k,coef;

/* Set up calculation, if not already done) */

	if((*memory) != 1) set_edge(a, b, ra, rb);

	for (i=0; i<3; i++)
	{
		real_to_gmp(c, i, c_mp[i+1]);
	}

	scalar_to_gmp(rc, rc_mp);
	build_weight(c_mp[1], c_mp[2], c_mp[3], rc_mp, c_mp[4]);

/* 
       Need to compute:
       Sc      : minor(a,b,c,i,j,0) for i=1,2 and j = i+1,3
       Tc      : minor(a,b,c,i,4,0) for i = 1,2,3
*/

	for (i=1; i<3 ; i++)
	{
		for (j=i+1; j<4; j++)
		{
			k=i+j-2;
			mpz_mul(temp1, c_mp[i], Dab[j]);
			mpz_mul(temp2, c_mp[j], Dab[i]);
			mpz_sub(temp1, temp1, temp2);
			mpz_add(Sc[k], temp1, Sab[k]);
		}
	}

	for (i=1; i<4; i++)
	{
		mpz_mul(temp1, c_mp[i], Dab[4]);
		mpz_mul(temp2, c_mp[4], Dab[i]);
		mpz_sub(temp1, temp1, temp2);
		mpz_add(Tc[i], temp1, Tab[i]);
	}

/*	This is the "hidden1" part */

	(*testa) = 0;

	if( mpz_cmp(a_mp[1], b_mp[1]) != 0) 
	{
		for (i = 1; i < 4 ; i++)
		{
			mpz_set(res[0][i], Dab[i]);
			mpz_set(res2_c[i][4], Tc[i]);
		}
		mpz_set(res[1][2], Sab[1]); mpz_set(res[1][3], Sab[2]);
		mpz_set(res[2][3], Sab[3]);
		mpz_set(res2_c[1][2], Sc[1]); mpz_set(res2_c[1][3], Sc[2]);
		mpz_set(res2_c[2][3], Sc[3]);
	}
	else if ( mpz_cmp(a_mp[2], b_mp[2]) != 0)
	{
		mpz_set(res[0][1], Dab[2]); mpz_set(res[0][2], Dab[3]);
		mpz_set(res[0][3], Dab[1]);
		mpz_set(res[1][2], Sab[3]);
		mpz_neg(res[1][3], Sab[1]); mpz_neg(res[2][3], Sab[2]);
		mpz_set(res2_c[1][2], Sc[3]);
		mpz_neg(res2_c[1][3], Sc[1]); mpz_neg(res2_c[2][3], Sc[2]);
		mpz_set(res2_c[1][4], Tc[2]); mpz_set(res2_c[2][4], Tc[3]);
		mpz_set(res2_c[3][4], Tc[1]);
	}
	else if (  mpz_cmp(a_mp[3], b_mp[3]) != 0)
	{
		mpz_set(res[0][1], Dab[3]); mpz_set(res[0][2], Dab[1]);
		mpz_set(res[0][3], Dab[2]);
		mpz_neg(res[1][2], Sab[2]);
		mpz_neg(res[1][3], Sab[3]); mpz_set(res[2][3], Sab[1]);
		mpz_neg(res2_c[1][2], Sc[2]);
		mpz_neg(res2_c[1][3], Sc[3]); mpz_set(res2_c[2][3], Sc[1]);
		mpz_set(res2_c[1][4], Tc[3]); mpz_set(res2_c[2][4], Tc[1]);
		mpz_set(res2_c[3][4], Tc[2]);
	}
	else
	{
		exit(1);
	}

	mpz_mul(r_11, res[0][1], res[0][1]);
	mpz_mul(r_22, res[0][2], res[0][2]);
	mpz_mul(r_33, res[0][3], res[0][3]);
	mpz_mul(temp1, res[0][3], res[1][2]); 
	mpz_mul(temp2, res[0][2], res[1][3]);
	mpz_sub(diff, temp1, temp2);

/* Compute det0 */

	mpz_add(temp1, r_22, r_33); mpz_add(temp1, temp1, r_11); 
	mpz_mul(temp1, temp1, res[0][1]); 
	coef = -2;
	mpz_mul_si(det0, temp1, coef);


/* Now check if edge (ab) is attached to c */

	mpz_mul(temp1, res[1][2], res2_c[1][2]);
	mpz_mul(temp2, res[1][3], res2_c[1][3]);
	mpz_add(temp1, temp1, temp2);
	coef = -2;
	mpz_mul_si(temp1, temp1, coef);

	mpz_set_si(temp2, 0);
	for (i=1; i<4; i++)
	{
		mpz_mul(temp3, res[0][i], res2_c[i][4]);
		mpz_add(temp2, temp2, temp3);
	}
	mpz_add(temp1, temp2, temp1);
	mpz_mul(temp1, temp1, res[0][1]);

	mpz_mul(temp2, res2_c[2][3], diff);
	mpz_mul_si(temp2, temp2, coef);
	mpz_sub(temp3, temp1, temp2);
	mpz_mul(dtest, temp3, det0);

	if(mpz_sgn(dtest) < 0) (*testa = 1);

}

/* =======================================================================================
 edge_radius_gmp: checks if the radius of the smallest circumsphere of an edge
	of the regular triangulation is smaller than the value of alpha
 ======================================================================================= */

void ALFCX_GMP::edge_radius_gmp(double *a, double *b, double ra, double rb, 
			int *testr, double alpha, int *memory)

{
	int i, coef, ivalue;
	double value;

	value = alpha*scale; ivalue = (int) floor(value); 
	mpz_set_si(alp, ivalue);

	if((*memory) != 1) set_edge(a, b, ra, rb);

/*	This is the "hidden1" part */

	(*testr) = 0;
	mpz_set(res[0][4],Dab[4]);


	if( mpz_cmp(a_mp[1],b_mp[1]) != 0) 
	{
		for (i = 1; i < 4 ; i++)
		{
			mpz_set(res[0][i], Dab[i]);
			mpz_mul(temp1, b_mp[i], a_mp[4]);
			mpz_mul(temp2, a_mp[i], b_mp[4]);
			mpz_sub(res[i][4], temp2, temp1);
		}
		mpz_set(res[1][2], Sab[1]); mpz_set(res[1][3], Sab[2]);
		mpz_set(res[2][3], Sab[3]);
	}
	else if ( mpz_cmp(a_mp[2], b_mp[2]) != 0)
	{
		mpz_set(res[0][1], Dab[2]); mpz_set(res[0][2], Dab[3]);
		mpz_set(res[0][3], Dab[1]);
		mpz_set(res[1][2], Sab[3]);
		mpz_neg(res[1][3], Sab[1]); mpz_neg(res[2][3], Sab[2]);
		mpz_mul(temp1, a_mp[2], b_mp[4]); mpz_mul(temp2, b_mp[2], a_mp[4]);
		mpz_sub(res[1][4], temp1, temp2);
		mpz_mul(temp1, a_mp[3], b_mp[4]); mpz_mul(temp2, b_mp[3], a_mp[4]);
		mpz_sub(res[2][4], temp1, temp2);
		mpz_mul(temp1, a_mp[1], b_mp[4]); mpz_mul(temp1, b_mp[1], a_mp[4]);
		mpz_sub(res[3][4], temp1, temp2);
	}
	else if (  mpz_cmp(a_mp[3], b_mp[3]) != 0)
	{
		mpz_set(res[0][1], Dab[3]); mpz_set(res[0][2], Dab[1]);
		mpz_set(res[0][3], Dab[2]);
		mpz_neg(res[1][2], Sab[2]);
		mpz_neg(res[1][3], Sab[3]); mpz_set(res[2][3], Sab[1]);
		mpz_mul(temp1, a_mp[3], b_mp[4]); mpz_mul(temp2, b_mp[3], a_mp[4]);
		mpz_sub(res[1][4], temp1, temp2);
		mpz_mul(temp1, a_mp[1], b_mp[4]); mpz_mul(temp2, b_mp[1], a_mp[4]);
		mpz_sub(res[2][4], temp1, temp2);
		mpz_mul(temp1, a_mp[2], b_mp[4]); mpz_mul(temp1, b_mp[2], a_mp[4]);
		mpz_sub(res[3][4], temp1, temp2);
	}
	else
	{
		exit(1);
	}

	mpz_mul(r_11, res[0][1], res[0][1]);
	mpz_mul(r_22, res[0][2], res[0][2]);
	mpz_mul(r_33, res[0][3], res[0][3]);
	mpz_mul(r_14, res[0][1], res[0][4]);
	mpz_mul(r_313, res[0][3], res[1][3]);
	mpz_mul(r_212, res[0][2], res[1][2]);
	mpz_mul(temp1, res[0][3], res[1][2]); 
	mpz_mul(temp2, res[0][2], res[1][3]);
	mpz_sub(diff, temp1, temp2);

/* Compute det0 */

	mpz_add(temp1, r_22, r_33); mpz_add(temp1, temp1, r_11); 
	mpz_mul(temp1, temp1, res[0][1]); 
	coef = -2;
	mpz_mul_si(det0, temp1, coef);

/* Compute det1 */

	mpz_add(temp1, r_313, r_212);
	coef = 2;
	mpz_mul_si(temp1, temp1, coef);
	mpz_sub(temp1, temp1, r_14);
	mpz_mul(det1, res[0][1], temp1);


/* Compute det2 */

	mpz_add(temp1, r_11, r_33);
	mpz_mul(temp1, temp1, res[1][2]);
	coef = -2;
	mpz_mul_si(temp1, temp1, coef);
	mpz_mul_si(temp2, r_313, coef);
	mpz_add(temp2, temp2, r_14);
	mpz_mul(temp2, temp2, res[0][2]);
	mpz_sub(det2, temp1, temp2);

/* Compute det3 */

	mpz_add(temp1, r_11, r_22);
	mpz_mul(temp1, temp1, res[1][3]);
	mpz_mul_si(temp1, temp1, coef);
	mpz_mul_si(temp2, r_212, coef);
	mpz_add(temp2, temp2, r_14);
	mpz_mul(temp2, temp2, res[0][3]);
	mpz_sub(det3, temp1, temp2);

/* Compute det4 */

	mpz_mul(temp1, res[0][3], res[3][4]);
	mpz_mul(temp2, res[0][2], res[2][4]);
	mpz_mul(temp3, res[0][1], res[1][4]);
	mpz_add(temp1, temp1, temp2); mpz_add(temp1, temp3, temp1);
	coef = 2;
	mpz_mul_si(temp1, temp1, coef);
	mpz_mul(temp1, temp1, res[0][1]);
	mpz_mul(temp2, res[1][3], res[1][3]);
	mpz_mul(temp3, res[1][2], res[1][2]);
	mpz_add(temp2, temp3, temp2);
	mpz_mul(temp2, temp2, res[0][1]);
	mpz_mul(temp3, res[2][3], diff);
	mpz_sub(temp2, temp3, temp2);
	coef = 4;
	mpz_mul_si(temp2, temp2, coef);
	mpz_add(det4, temp1, temp2);

/* Compute numerator of the radius of the smallest circumsphere of the edge */

	mpz_mul(temp1, det0, det4);
	mpz_mul(temp2, det3, det3);
	mpz_sub(temp2, temp2, temp1);
	mpz_mul(temp1, det2, det2);
	mpz_add(temp2, temp2, temp1);
	mpz_mul(temp1, det1, det1);
	mpz_add(num, temp1, temp2);

/* Compute denominator of the radius of the smallest circumsphere of the edge */

	mpz_mul(den, det0, det0);

/* check if radius is lower than ALPHA         */

	mpz_mul(temp1, den, alp);
	mpz_sub(temp2, num, temp1);

	if(mpz_sgn(temp2) < 0) (*testr)=1;

}

/* =======================================================================================
 triangle_radius_gmp: checks if the radius of the circumsphere of a facet of the
   regular triangulation is smaller than alpha

       For the three points a,b,c that form the triangles, the program
       needs as input the following determinants:

       S(i,j)   = Minor(a,b,c,i,j,0)= det | a(i)  a(j)  1 |
                                          | b(i)  b(j)  1 |
                                          | c(i)  c(j)  1 |

       for i in [1,3] and j in [i+1,4]

       and:

       T(i,j) = Minor(a,b,c,i,j,4)=det | a(i) a(j) a(4) |
                                       | b(i) b(j) b(4) |
                                       | c(i) c(j) c(4) |

       and

       Dabc  = Minor(a,b,c,1,2,3)

	Output:

	testr	: flag set to 1 if ALPHA is larger than rho, the radius
		  of the circumsphere of the triangle
 ======================================================================================= */

void ALFCX_GMP::triangle_radius_gmp(double *a, double *b, double *c, double ra, 
		double rb, double rc, int *testr, double alpha, int *memory)
{
	int i, j, coef,  ivalue;
	double value;

	value = alpha*scale; ivalue = (int) floor(value); 
	mpz_set_si(alp,  ivalue);

	if(*memory !=1) set_triangle(a,  b,  c,  ra,  rb,  rc);

	(*testr) = 0;

	mpz_set_si(temp1, 0);
	for (i=1; i<3; i++)
	{
		for (j=i+1; j<4; j++)
		{
			mpz_mul(temp2, S[i][j], S[i][j]);
			mpz_add(temp1, temp1, temp2);
		}
	}

/* Compute det0 */

	coef = 4;
	mpz_mul_si(det0, temp1, coef);

/* Compute det1 */

	mpz_mul(temp1, Dabc, S[2][3]);
	coef = -2;
	mpz_mul_si(temp1, temp1, coef);
	mpz_mul(temp2, S[1][2], S[2][4]);
	mpz_add(temp1, temp2, temp1);
	mpz_mul(temp2, S[1][3], S[3][4]);
	mpz_add(temp1, temp2, temp1);
	mpz_mul_si(det1, temp1, coef);

/* Compute det2 */

	coef = 2;
	mpz_mul(temp1, Dabc, S[1][3]);
	mpz_mul_si(temp1, temp1, coef);
	mpz_mul(temp2, S[2][3], S[3][4]);
	mpz_add(temp1, temp2, temp1);
	mpz_mul(temp2, S[1][2], S[1][4]);
	mpz_sub(temp1, temp2, temp1);
	mpz_mul_si(det2, temp1, coef);

/* Compute det3 */

	mpz_mul(temp1, Dabc, S[1][2]);
	mpz_mul_si(temp1, temp1, coef);
	mpz_mul(temp2, S[1][3], S[1][4]);
	mpz_add(temp1, temp2, temp1);
	mpz_mul(temp2, S[2][3], S[2][4]);
	mpz_add(temp1, temp2, temp1);
	mpz_mul_si(det3, temp1, coef);

/* Compute det4 */

	mpz_mul(temp1, Dabc, Dabc);
	coef = -2;
	mpz_mul_si(temp1, temp1, coef);

	for (i=1; i<3; i++)
	{
		for (j=i+1; j<4; j++)
		{
			mpz_mul(temp2, S[i][j], T[i][j]);
			mpz_add(temp1, temp1, temp2);
		}
	}
	coef = -4;
	mpz_mul_si(det4, temp1, coef);

/* Now compute numerator of the radius of the circumsphere of the triangle */

	mpz_mul(temp1, det0, det4);
	mpz_mul(temp2, det3, det3);
	mpz_sub(temp2, temp2, temp1);
	mpz_mul(temp1, det2, det2);
	mpz_add(temp2, temp2, temp1);
	mpz_mul(temp1, det1, det1);
	mpz_add(num, temp1, temp2);

/* Now compute denominator of the radius of the circumsphere of the triangle */

	mpz_mul(den, det0, det0);

/* Check if radius is lower than ALPHA */

	mpz_mul(temp1, den, alp);
	mpz_sub(temp2, num, temp1);

	if(mpz_sgn(temp2) < 0) (*testr) = 1;

}

/* =======================================================================================
 triangle_attach_gmp: checks if a facet is attached to a vertex

	Input:

	For the three points a,b,c that form the triangles, the program
	needs as input the following determinants:

         S(i,j) = Minor(a,b,c,i,j,0)= det | a(i)  a(j)  1 |
                                          | b(i)  b(j)  1 |
                                          | c(i)  c(j)  1 |
       for all i in [1,3], j in [i+1,4]

       T(i,j) = M(a,b,c,i,j,4)    = Det | a(i) a(j) a(4) |
                                        | b(i) b(j) b(4) |
                                        | c(i) c(j) c(4) |

       for all i in [1,2] and all j in [i+1,3]

       Dabc = Det | a(1) a(2) a(3) |
                  | b(1) b(2) b(3) |
                  | c(1) c(2) c(3) |

       and the coordinates of the fourth vertex d


	testa	: flag set to 1 if the fourth point d is inside the
		  circumsphere of (a,b,c)
 ======================================================================================= */

void ALFCX_GMP::triangle_attach_gmp(double *a, double *b, double *c, double *d,
        double ra, double rb, double rc, double rd, int *testa, int *memory)

{
	int i,coef;

	if(*memory !=1) set_triangle(a, b, c, ra, rb, rc);

	for (i=0; i<3; i++)
	{
		real_to_gmp(d, i, d_mp[i+1]);
	}
	scalar_to_gmp(rd, rd_mp);
	build_weight(d_mp[1], d_mp[2], d_mp[3], rd_mp, d_mp[4]);

/*
       We need to compute:

       det1 = Minor(a,b,c,d,2,3,4,0)
       det2 = Minor(a,b,c,d,1,3,4,0)
       det3 = Minor(a,b,c,d,1,2,4,0)
       det4 = Minor(a,b,c,d,1,2,3,0)
*/
	mpz_mul(temp1, d_mp[2], S[3][4]); mpz_mul(temp2, d_mp[3], S[2][4]); mpz_sub(temp1, temp2, temp1);
	mpz_mul(temp2, d_mp[4], S[2][3]); mpz_sub(temp2, T[2][3], temp2); mpz_add(det1, temp2, temp1);
	mpz_mul(temp1, d_mp[1], S[3][4]); mpz_mul(temp2, d_mp[3], S[1][4]); mpz_sub(temp1, temp2, temp1);
	mpz_mul(temp2, d_mp[4], S[1][3]); mpz_sub(temp2, T[1][3], temp2); mpz_add(det2, temp2, temp1);
	mpz_mul(temp1, d_mp[1], S[2][4]); mpz_mul(temp2, d_mp[2], S[1][4]); mpz_sub(temp1, temp2, temp1);
	mpz_mul(temp2, d_mp[4], S[1][2]); mpz_sub(temp2, T[1][2], temp2); mpz_add(det3, temp2, temp1);
	mpz_mul(temp1, d_mp[1], S[2][3]); mpz_mul(temp2, d_mp[2], S[1][3]); mpz_sub(temp1, temp2, temp1);
	mpz_mul(temp2, d_mp[3], S[1][2]); mpz_sub(temp2, Dabc, temp2); mpz_add(det4, temp2, temp1);

/* check if triangle is attached                               */

	(*testa) = 0;

	mpz_set_si(temp1, 1);
/*	for (i=1; i<4; i++)
	{
		mpz_mul(temp2, S[i], S[i]);
		mpz_add(temp1, temp1, temp2);
	}
*/
	mpz_mul(temp2, det4, Dabc);
	coef = -2;
	mpz_mul_si(temp2, temp2, coef);
	mpz_mul(temp3, det3, S[1][2]);
	mpz_add(temp2, temp3, temp2);
	mpz_mul(temp3, det2, S[1][3]);
	mpz_add(temp2, temp3, temp2);
	mpz_mul(temp3, det1, S[2][3]);
	mpz_add(temp2,temp3,temp2);
	mpz_mul(dtest,temp1,temp2);

	if(mpz_sgn(dtest) > 0) (*testa)=1;

/* 	Clear local GMP variables */

}
/* =======================================================================================
 set_alf_gmp: initialises all gmp variables that can be used for computing the dual complex
 ======================================================================================= */

void ALFCX_GMP::set_alf_gmp()
{
/*	Initialise local GMP variables */

	int i,j;

	mpz_init(temp1); mpz_init(temp2), mpz_init(temp3);
	mpz_init(ra2); mpz_init(rb2); mpz_init(dist2);
	mpz_init(dtest);
	mpz_init (r_11); mpz_init (r_22); mpz_init (r_33);
	mpz_init (diff);
	mpz_init (det0);

	for (i= 0; i < 4; i++)
	{
		for (j=0; j < 5; j++)
		{
			mpz_init(res[i][j]);
			mpz_init(res2_c[i][j]);
			mpz_init(Mab[i][j]);
			mpz_init(Mac[i][j]);
			mpz_init(Mbc[i][j]);
			mpz_init(S[i][j]);
			mpz_init(T[i][j]);
		}
	}
	mpz_init (r_14); mpz_init (r_313); mpz_init (r_212);
	mpz_init (det1); mpz_init (det2); mpz_init (det3); mpz_init (det4); 
	mpz_init (wa); mpz_init(wb); mpz_init(wc); mpz_init(wd);

	for (i = 0; i < 5; i++) 
	{
		mpz_init(a_mp[i]);
		mpz_init(b_mp[i]);
		mpz_init(c_mp[i]);
		mpz_init(d_mp[i]);
		mpz_init(Dab[i]);
	}
	mpz_init(ra_mp);mpz_init(rb_mp);
	mpz_init(rc_mp);mpz_init(rd_mp);

	mpz_init (num); mpz_init (den);

	for (i=0; i < 4; i++)
	{
		mpz_init(Sab[i]);
		mpz_init(Sac[i]);
		mpz_init(Sad[i]);
		mpz_init(Sbc[i]);
		mpz_init(Sbd[i]);
		mpz_init(Scd[i]);
		mpz_init(Tab[i]);
		mpz_init(Sa[i]);
		mpz_init(Sb[i]);
		mpz_init(Sc[i]);
		mpz_init(Sd[i]);
		mpz_init(Sam1[i]);
		mpz_init(Sbm1[i]);
		mpz_init(Scm1[i]);
		mpz_init(Sdm1[i]);
		mpz_init(Tc[i]);
		mpz_init(Deter[i]);
	}

	mpz_init(Dabc); mpz_init(Dabd); mpz_init(Dacd); mpz_init(Dbcd);
	mpz_init(alp);mpz_init(Dabcd);

}

/* =======================================================================================
 clear_alf_gmp:	clears all gmp variables that were used for computing the dual complex
 ======================================================================================= */

void ALFCX_GMP::clear_alf_gmp()
{

	int i,j;

	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3);
	mpz_clear(ra2); mpz_clear(rb2); mpz_clear(dist2);
	mpz_clear(dtest);
	mpz_clear (r_11); mpz_clear (r_22); mpz_clear (r_33);
	mpz_clear (diff);
	mpz_clear (det0);

	for (i= 0; i < 4; i++)
	{
		for (j=0; j < 5; j++)
		{
			mpz_clear(res[i][j]);
			mpz_clear(res2_c[i][j]);
			mpz_clear(Mab[i][j]);
			mpz_clear(Mac[i][j]);
			mpz_clear(Mbc[i][j]);
			mpz_clear(S[i][j]);
			mpz_clear(T[i][j]);
		}
	}
	
	mpz_clear (r_14); mpz_clear (r_313); mpz_clear (r_212);
	mpz_clear (det1); mpz_clear (det2); mpz_clear (det3); mpz_clear (det4); 
	mpz_clear (wa); mpz_clear(wb); mpz_clear(wc); mpz_clear(wd);

	for (i = 0; i < 5; i++) 
	{
		mpz_clear(a_mp[i]);
		mpz_clear(b_mp[i]);
		mpz_clear(c_mp[i]);
		mpz_clear(d_mp[i]);
		mpz_clear(Dab[i]);
	}
	mpz_clear(ra_mp);mpz_clear(rb_mp);
	mpz_clear(rc_mp);mpz_clear(rd_mp);

	mpz_clear (num); mpz_clear (den);

	for (i=0; i < 4; i++)
	{
		mpz_clear(Sab[i]);
		mpz_clear(Sac[i]);
		mpz_clear(Sad[i]);
		mpz_clear(Sbc[i]);
		mpz_clear(Sbd[i]);
		mpz_clear(Scd[i]);
		mpz_clear(Tab[i]);
		mpz_clear(Sa[i]);
		mpz_clear(Sb[i]);
		mpz_clear(Sc[i]);
		mpz_clear(Sd[i]);
		mpz_clear(Sam1[i]);
		mpz_clear(Sbm1[i]);
		mpz_clear(Scm1[i]);
		mpz_clear(Sdm1[i]);
		mpz_clear(Tc[i]);
		mpz_clear(Deter[i]);
	}

	mpz_clear(Dabc);mpz_clear(Dabd); mpz_clear(Dacd); mpz_clear(Dbcd);
	mpz_clear(Dabcd);
	mpz_clear(alp);

}

#endif
