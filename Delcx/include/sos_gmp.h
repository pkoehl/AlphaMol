/* ==============================================================================
 *	Sos_gmp
 *									  
 *  Performs all operations   
 *  with multi precision arithmetics, using the package GMP               
 *									  
   ============================================================================== */

#ifndef SOS_GMP_H
#define SOS_GMP_H

  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include "gmp.h"

  #define round(x) ((x)>=0?(int)((x)+0.5):(int)((x)-0.5))

/* ==============================================================================
   The SOS class
   ============================================================================== */

  class SOS {
  
  public:

	void real_to_gmp(double coord, mpz_t val);

	void build_weight(mpz_t ax, mpz_t ay, mpz_t az, mpz_t r, mpz_t w);

	void sos_minor2_gmp(double xa, double xb, int *res);

	void sos_minor3_gmp(double xa, double ya, double xb, double yb, double xc, double yc, int *res);

	void sos_minor4_gmp(double *coord_a, double *coord_b, double *coord_c, double *coord_d, int *res);

	void sos_minor5_gmp(double *coord_a, double ra, double *coord_b, double rb, 
	double *coord_c, double rc, double *coord_d, double rd, double *coord_e, double re,
	int *res) ;

	void minor4_gmp(double *coord_a, double *coord_b, double *coord_c, double *coord_d, int *res);

	void init_sos_gmp();

	void clear_sos_gmp();

  private:

	void deter2_gmp(mpz_t deter, mpz_t a, mpz_t b);

	void deter3_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a21, mpz_t a22,
		mpz_t a31, mpz_t a32);

	void deter4_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a13, mpz_t a21,
		mpz_t a22, mpz_t a23, mpz_t a31, mpz_t a32,
		mpz_t a33, mpz_t a41, mpz_t a42, mpz_t a43);

	void deter5_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a13, mpz_t a14,
		mpz_t a21, mpz_t a22, mpz_t a23, mpz_t a24,
		mpz_t a31, mpz_t a32, mpz_t a33, mpz_t a34,
		mpz_t a41, mpz_t a42, mpz_t a43, mpz_t a44,
		mpz_t a51, mpz_t a52, mpz_t a53, mpz_t a54);

	mpz_t a11_mp,a12_mp,a13_mp,a14_mp;
	mpz_t a21_mp,a22_mp,a23_mp,a24_mp;
	mpz_t a31_mp,a32_mp,a33_mp,a34_mp;
	mpz_t a41_mp,a42_mp,a43_mp,a44_mp;
	mpz_t a51_mp,a52_mp,a53_mp,a54_mp;
	mpz_t r1_mp,r2_mp, r3_mp, r4_mp, r5_mp;

	mpz_t temp1,temp2,temp3,temp4;
	mpz_t val1,val2,val3;

	mpz_t c11,c12,c13,c14,c21,c22,c23,c24,c31,c32,c33,c34,c41,c42,c43,c44;
	mpz_t d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3;

	double scale;

  };

/* ==============================================================================
   init_sos_gmp:						  
   Initilize all private gmp variables needed by delcx
   ============================================================================== */

  void SOS::init_sos_gmp()
  {

	mpz_init(a11_mp);mpz_init(a12_mp); mpz_init(a13_mp); mpz_init(a14_mp);
	mpz_init(a21_mp);mpz_init(a22_mp); mpz_init(a23_mp); mpz_init(a24_mp);
	mpz_init(a31_mp);mpz_init(a32_mp); mpz_init(a33_mp); mpz_init(a34_mp);
	mpz_init(a41_mp);mpz_init(a42_mp); mpz_init(a43_mp); mpz_init(a44_mp);
	mpz_init(a51_mp);mpz_init(a52_mp); mpz_init(a53_mp); mpz_init(a54_mp);

	mpz_init(r1_mp);mpz_init(r2_mp); mpz_init(r3_mp); mpz_init(r4_mp); mpz_init(r5_mp);

	mpz_init(temp1); mpz_init(temp2); mpz_init(temp3); mpz_init(temp4);
	mpz_init(val1);mpz_init(val2); mpz_init(val3);

	mpz_init(c11); mpz_init(c12); mpz_init(c13); mpz_init(c14);
	mpz_init(c21); mpz_init(c22); mpz_init(c23); mpz_init(c24);
	mpz_init(c31); mpz_init(c32); mpz_init(c33); mpz_init(c34);
	mpz_init(c41); mpz_init(c42); mpz_init(c43); mpz_init(c44);
	mpz_init(d1); mpz_init(d2); mpz_init(d3);
	mpz_init(e1); mpz_init(e2); mpz_init(e3);
	mpz_init(f1); mpz_init(f2); mpz_init(f3);
	mpz_init(g1); mpz_init(g2); mpz_init(g3);

	scale = 1.e8;
  }

/* ==============================================================================
   clear_sos_gmp:						  
   clear all private gmp variables used by delcx
   ============================================================================== */

  void SOS::clear_sos_gmp()
  {

	mpz_clear(a11_mp);mpz_clear(a12_mp); mpz_clear(a13_mp); mpz_clear(a14_mp);
	mpz_clear(a21_mp);mpz_clear(a22_mp); mpz_clear(a23_mp); mpz_clear(a24_mp);
	mpz_clear(a31_mp);mpz_clear(a32_mp); mpz_clear(a33_mp); mpz_clear(a34_mp);
	mpz_clear(a41_mp);mpz_clear(a42_mp); mpz_clear(a43_mp); mpz_clear(a44_mp);
	mpz_clear(a51_mp);mpz_clear(a52_mp); mpz_clear(a53_mp); mpz_clear(a54_mp);

	mpz_clear(r1_mp);mpz_clear(r2_mp); mpz_clear(r3_mp); mpz_clear(r4_mp); mpz_clear(r5_mp);

	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3); mpz_clear(temp4);
	mpz_clear(val1);mpz_clear(val2); mpz_clear(val3);

	mpz_clear(c11); mpz_clear(c12); mpz_clear(c13); mpz_clear(c14);
	mpz_clear(c21); mpz_clear(c22); mpz_clear(c23); mpz_clear(c24);
	mpz_clear(c31); mpz_clear(c32); mpz_clear(c33); mpz_clear(c34);
	mpz_clear(c41); mpz_clear(c42); mpz_clear(c43); mpz_clear(c44);
	mpz_clear(d1); mpz_clear(d2); mpz_clear(d3);
	mpz_clear(e1); mpz_clear(e2); mpz_clear(e3);
	mpz_clear(f1); mpz_clear(f2); mpz_clear(f3);
	mpz_clear(g1); mpz_clear(g2); mpz_clear(g3);

  }
/* ==============================================================================
   real_to_gmp:						  
   converts one coordinate of one point into mpz_t type
   ============================================================================== */

  void SOS::real_to_gmp(double coord, mpz_t val)
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
 
/* ==============================================================================
 build_weight:
	builds the weight of a point:
	w = x**2 + y**2 + z**2 - ra**2
   ============================================================================== */

  void SOS::build_weight(mpz_t ax, mpz_t ay, mpz_t az, mpz_t r, mpz_t w)
  {
	mpz_mul(temp1,r,r);
	mpz_mul(temp2,ax,ax), mpz_sub(temp1,temp2,temp1);
	mpz_mul(temp2,ay,ay), mpz_add(temp1,temp2,temp1);
	mpz_mul(temp2,az,az), mpz_add(w,temp2,temp1);
  }

/* ==============================================================================
   deter2_gmp:								  
	This subroutine evaluates the determinant:
		D = | b11 1 |
		    | b21 1 |

	Input:
		b11, b21
	Output:
		deter
   ============================================================================== */

  void SOS::deter2_gmp(mpz_t deter, mpz_t b11, mpz_t b21)
  {
	mpz_sub(deter,b11,b21);
  }

/* ==============================================================================
   deter3_gmp:								  
	This subroutine evaluates the determinant:
		D = | b11 b12 1 |
		    | b21 b22 1 |
		    | b31 b32 1 |

	Input:
		b11, b12, b21, b22, b31, b32
	Output:
		deter3_gmp
   ============================================================================== */

  void SOS::deter3_gmp(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b21, 
		mpz_t b22, mpz_t b31, mpz_t b32)
  {

	mpz_sub(temp1,b21,b11);
	mpz_sub(temp2,b22,b12);
	mpz_sub(temp3,b31,b11);
	mpz_sub(temp4,b32,b12);

	mpz_mul(val1,temp1,temp4);
	mpz_mul(val2,temp2,temp3);

	mpz_sub(deter,val1,val2);

  }

/* ==============================================================================
   deter4_gmp:								  
	This subroutine evaluates the determinant:
		D = | b11 b12 b13 1 |
		    | b21 b22 b23 1 |
		    | b31 b32 b33 1 |
		    | b41 b42 b43 1 |

	Input:
		b11, b12, b13, b21, b22, b23, b31, b32, b33
		b41, b42, b43
	Output:
		deter4_gmp
   ============================================================================== */

  void SOS::deter4_gmp(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13, mpz_t b21, 
		mpz_t b22, mpz_t b23, mpz_t b31, mpz_t b32, mpz_t b33, 
		mpz_t b41, mpz_t b42, mpz_t b43)
  {
	mpz_sub(c11,b21,b11);mpz_sub(c12,b22,b12);mpz_sub(c13,b23,b13);
	mpz_sub(c21,b31,b11);mpz_sub(c22,b32,b12);mpz_sub(c23,b33,b13);
	mpz_sub(c31,b41,b11);mpz_sub(c32,b42,b12);mpz_sub(c33,b43,b13);

	mpz_mul(temp1,c22,c33);mpz_mul(temp2,c32,c23);mpz_sub(val1,temp1,temp2);
	mpz_mul(temp1,c12,c33);mpz_mul(temp2,c32,c13);mpz_sub(val2,temp1,temp2);
	mpz_mul(temp1,c12,c23);mpz_mul(temp2,c22,c13);mpz_sub(val3,temp1,temp2);

	mpz_mul(temp1,c21,val2);mpz_mul(temp2,c11,val1);mpz_mul(temp3,c31,val3);

	mpz_add(val1,temp2,temp3);
	mpz_sub(deter,temp1,val1);

  }

/* ==============================================================================
   deter5_gmp:								  
	This subroutine evaluates the determinant:
		D = | b11 b12 b13 b14 1 |
		    | b21 b22 b23 b24 1 |
		    | b31 b32 b33 b34 1 |
		    | b41 b42 b43 b44 1 |
		    | b51 b52 b53 b54 1 |

	Input:
		b11, b12, b13, b14, b21, b22, b23, b24, b31, b32, b33, b34
		b41, b42, b43, b44, b51, b52, b53, b54
	Output:
		deter5_gmp
   ============================================================================== */

  void SOS::deter5_gmp(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13, mpz_t b14, 
		 mpz_t b21, mpz_t b22, mpz_t b23, mpz_t b24,
		 mpz_t b31, mpz_t b32, mpz_t b33, mpz_t b34,
		 mpz_t b41, mpz_t b42, mpz_t b43, mpz_t b44,
		 mpz_t b51, mpz_t b52, mpz_t b53, mpz_t b54)
  {

	mpz_sub(c11,b21,b11); mpz_sub(c12,b22,b12); mpz_sub(c13,b23,b13);
	mpz_sub(c14,b24,b14);
	mpz_sub(c21,b31,b11); mpz_sub(c22,b32,b12); mpz_sub(c23,b33,b13);
	mpz_sub(c24,b34,b14);
	mpz_sub(c31,b41,b11); mpz_sub(c32,b42,b12); mpz_sub(c33,b43,b13);
	mpz_sub(c34,b44,b14);
	mpz_sub(c41,b51,b11); mpz_sub(c42,b52,b12); mpz_sub(c43,b53,b13);
	mpz_sub(c44,b54,b14);

	mpz_mul(temp1,c32,c43); mpz_mul(temp2,c42,c33); mpz_sub(d1,temp1,temp2);
	mpz_mul(temp1,c32,c44); mpz_mul(temp2,c42,c34); mpz_sub(d2,temp1,temp2);
	mpz_mul(temp1,c33,c44); mpz_mul(temp2,c43,c34); mpz_sub(d3,temp1,temp2);

	mpz_mul(temp1,c12,c23); mpz_mul(temp2,c22,c13); mpz_sub(e1,temp1,temp2);
	mpz_mul(temp1,c12,c24); mpz_mul(temp2,c22,c14); mpz_sub(e2,temp1,temp2);
	mpz_mul(temp1,c13,c24); mpz_mul(temp2,c23,c14); mpz_sub(e3,temp1,temp2);

	mpz_mul(temp1,c11,c24); mpz_mul(temp2,c21,c14); mpz_sub(f1,temp1,temp2);
	mpz_mul(temp1,c11,c23); mpz_mul(temp2,c21,c13); mpz_sub(f2,temp1,temp2);
	mpz_mul(temp1,c11,c22); mpz_mul(temp2,c21,c12); mpz_sub(f3,temp1,temp2);

	mpz_mul(temp1,c31,c44); mpz_mul(temp2,c41,c34); mpz_sub(g1,temp1,temp2);
	mpz_mul(temp1,c31,c43); mpz_mul(temp2,c41,c33); mpz_sub(g2,temp1,temp2);
	mpz_mul(temp1,c31,c42); mpz_mul(temp2,c41,c32); mpz_sub(g3,temp1,temp2);
 
	mpz_mul(temp1,e3,g3); mpz_mul(temp2,e2,g2); mpz_sub(temp3,temp1,temp2);
	mpz_mul(temp1,e1,g1); mpz_add(temp3,temp3,temp1);
	mpz_mul(temp1,d3,f3); mpz_add(temp3,temp3,temp1);
	mpz_mul(temp1,d2,f2); mpz_sub(temp3,temp3,temp1);
	mpz_mul(temp1,d1,f1); mpz_add(deter,temp3,temp1);

  }


/* ==============================================================================
 sos_minor2_gmp:							  
	This subroutine tests the sign of the determinant
		D = | a11 1 |
		    | a21 1 |
	If the determinant is found to be 0, then the SoS procedure is used:
	a development of the determinant with respect to a perturbation EPS
	applied to the coordinates in the determinant is computed, and
	the sign of the first non zero term defines the sign of the 
	determinant.
	In the case of a 2x2 determinant, the first term in the expansion
	is the coefficient 1 ...					
   ============================================================================== */

  void SOS::sos_minor2_gmp(double xa, double xb, int *res)
  {
	int icomp;

/* Get coordinates */

	real_to_gmp(xa,a11_mp);
	real_to_gmp(xb,a21_mp);

/* Compute determinant */

	deter2_gmp(temp1,a11_mp,a21_mp);

	icomp = mpz_sgn(temp1);

	if (icomp != 0) {
		*res = icomp;
	}
	else {
		*res = 1;
	}

  }

/* ==============================================================================
   sos_minor3_gmp:							  
	This subroutine tests the sign of the determinant
		D = | a11 a12 1 |
		    | a21 a22 1 |
		    | a31 a32 1 |
	If the determinant is found to be 0, then the SoS procedure is used:
	a development of the determinant with respect to a perturbation EPS
	applied to the coordinates in the determinant is computed, and
	the sign of the first non zero term defines the sign of the 
	determinant.
	In the case of a 3x3 determinant, the maximum number of terms to be
	checked is 4 ...					
   ============================================================================== */

  void SOS::sos_minor3_gmp(double xa, double ya, double xb, double yb, double xc, double yc, int *res)
  {
	int icomp;

/* Transfer coordinates to GMP */

	real_to_gmp(xa,a11_mp);
	real_to_gmp(ya,a12_mp);
	real_to_gmp(xb,a21_mp);
	real_to_gmp(yb,a22_mp);
	real_to_gmp(xc,a31_mp);
	real_to_gmp(yc,a32_mp);

/* Compute determinant */

	deter3_gmp(temp1,a11_mp,a12_mp,a21_mp,a22_mp,a31_mp,a32_mp);

	icomp = mpz_sgn(temp1);

/* if major determinant is non 0, return its sign */

	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Look now at each term in the expansion of the determinant with
   respect to EPS                                               
	The initial determinant is:
		Minor3(i,j,k,1,2,0)				*/

/* Term 1: - Minor2(j,k,1,0) */

	deter2_gmp(temp1,a21_mp,a31_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = - icomp;
		return;
	}

/* Term 2: Minor2(j,k,2,0)  */

	deter2_gmp(temp1,a22_mp,a32_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Term 3: Minor2(i,k,1,0)  */

	deter2_gmp(temp1,a11_mp,a31_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Term 4: 1 */

	*res = 1;

}
/* ==============================================================================
   sos_minor4_:								  
	This subroutine tests the sign of the determinant
		D = | a11 a12 a13 1 |
		    | a21 a22 a23 1 |
		    | a31 a32 a33 1 |
		    | a41 a42 a43 1 |
	If the determinant is found to be 0, then the SoS procedure is used:
	a development of the determinant with *respect to a perturbation EPS
	applied to the coordinates in the determinant is computed, and
	the sign of the first non zero term defines the sign of the 
	determinant.
	In the case of a 4x4 determinant, the maximum number of terms to be
	checked is 14 ...					
   ============================================================================== */

  void SOS::sos_minor4_gmp(double *coord_a, double *coord_b, double *coord_c, double *coord_d, int *res)
  {
	int icomp;

/* Transfer coordinates to gmp */

	real_to_gmp(coord_a[0],a11_mp);
	real_to_gmp(coord_a[1],a12_mp);
	real_to_gmp(coord_a[2],a13_mp);
	real_to_gmp(coord_b[0],a21_mp);
	real_to_gmp(coord_b[1],a22_mp);
	real_to_gmp(coord_b[2],a23_mp);
	real_to_gmp(coord_c[0],a31_mp);
	real_to_gmp(coord_c[1],a32_mp);
	real_to_gmp(coord_c[2],a33_mp);
	real_to_gmp(coord_d[0],a41_mp);
	real_to_gmp(coord_d[1],a42_mp);
	real_to_gmp(coord_d[2],a43_mp);

/* Compute determinant */

	deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
			   a31_mp,a32_mp,a33_mp,a41_mp,a42_mp,a43_mp);

//	printf("Deter4 = %s\n",mpz_get_str(NULL,base,temp1)); 

	icomp = mpz_sgn(temp1);

/* if major determinant is non 0, return its sign 
   (don't forget to clear GMP variables !!)	*/

	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Look now at each term in the expansion of the determinant with
   *respect to EPS                                                
	The initial determinant is:
		Minor4(i,j,k,l,1,2,3,0)				*/

/* Term 1:	Minor3(j,k,l,1,2,0)  */

	deter3_gmp(temp1,a21_mp,a22_mp,a31_mp,a32_mp,a41_mp,a42_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Term 2:	-Minor3(j,k,l,1,3,0) */

	deter3_gmp(temp1,a21_mp,a23_mp,a31_mp,a33_mp,a41_mp,a43_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = - icomp;
		return;
	}

/* Term 3:	Minor3(j,k,l,2,3,0) */

	deter3_gmp(temp1,a22_mp,a23_mp,a32_mp,a33_mp,a42_mp,a43_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Term 4:	- Minor3(i,k,l,1,2,0) */

	deter3_gmp(temp1,a11_mp,a12_mp,a31_mp,a32_mp,a41_mp,a42_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = - icomp;
		return;
	}

/* Term 5:	Minor2(k,l,1,0) */

	deter2_gmp(temp1,a31_mp,a41_mp);	
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Term 6:	-Minor2(k,l,2,0) */

	deter2_gmp(temp1,a32_mp,a42_mp);	
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = - icomp;
		return;
	}

/* Term 7:	Minor3(i,k,l,1,3,0) */

	deter3_gmp(temp1,a11_mp,a13_mp,a31_mp,a33_mp,a41_mp,a43_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 8:	Minor2(k,l,3,0) */

	deter2_gmp(temp1,a33_mp,a43_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 9:	- Minor3(i,k,l,2,3,0) */

	deter3_gmp(temp1,a12_mp,a13_mp,a32_mp,a33_mp,a42_mp,a43_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  - icomp;
		return;
	}

/* Term 10:	Minor3(i,j,l,1,2,0) */

	deter3_gmp(temp1,a11_mp,a12_mp,a21_mp,a22_mp,a41_mp,a42_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 11: 	- Minor2(j,l,1,0) */

	deter2_gmp(temp1,a21_mp,a41_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  - icomp;
		return;
	}

/* Term 12:	Minor2(j,l,2,0) */

	deter2_gmp(temp1,a22_mp,a42_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 13:	Minor2(i,l,1,0) */

	deter2_gmp(temp1,a11_mp,a41_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 14:	1 */

	*res = 1;
}

/* ==============================================================================
   sos_minor5_gmp:							  
	This subroutine tests the sign of the determinant
		D = | a11 a12 a13 a14 1 |
		    | a21 a22 a23 a24 1 |
		    | a31 a32 a33 a34 1 |
		    | a41 a42 a43 a44 1 |
		    | a51 a52 a53 a54 1 |
	If the determinant is found to be 0, then the SoS procedure is used:
	a development of the determinant with *respect to a perturbation EPS
	applied to the coordinates in the determinant is computed, and
	the sign of the first non zero term defines the sign of the 
	determinant.
	In the case of a 5x5 determinant, the maximum number of terms to be
	checked is 49 ...					
   ============================================================================== */

  void SOS::sos_minor5_gmp(double *coord_a, double ra, double *coord_b, double rb, 
	double *coord_c, double rc, double *coord_d, double rd, double *coord_e, double re,
	int *res) 
{
	int icomp;

/*	Initialise local GMP variables */

	real_to_gmp(coord_a[0],a11_mp);
	real_to_gmp(coord_a[1],a12_mp);
	real_to_gmp(coord_a[2],a13_mp);
	real_to_gmp(coord_b[0],a21_mp);
	real_to_gmp(coord_b[1],a22_mp);
	real_to_gmp(coord_b[2],a23_mp);
	real_to_gmp(coord_c[0],a31_mp);
	real_to_gmp(coord_c[1],a32_mp);
	real_to_gmp(coord_c[2],a33_mp);
	real_to_gmp(coord_d[0],a41_mp);
	real_to_gmp(coord_d[1],a42_mp);
	real_to_gmp(coord_d[2],a43_mp);
	real_to_gmp(coord_e[0],a51_mp);
	real_to_gmp(coord_e[1],a52_mp);
	real_to_gmp(coord_e[2],a53_mp);

	real_to_gmp(ra,r1_mp);
	real_to_gmp(rb,r2_mp);
	real_to_gmp(rc,r3_mp);
	real_to_gmp(rd,r4_mp);
	real_to_gmp(re,r5_mp);

	build_weight(a11_mp,a12_mp,a13_mp,r1_mp,a14_mp);
	build_weight(a21_mp,a22_mp,a23_mp,r2_mp,a24_mp);
	build_weight(a31_mp,a32_mp,a33_mp,r3_mp,a34_mp);
	build_weight(a41_mp,a42_mp,a43_mp,r4_mp,a44_mp);
	build_weight(a51_mp,a52_mp,a53_mp,r5_mp,a54_mp);

/* Compute determinant */

	deter5_gmp(temp1,a11_mp,a12_mp,a13_mp,a14_mp,a21_mp,a22_mp,
	  a23_mp,a24_mp,a31_mp,a32_mp,a33_mp,a34_mp,a41_mp,a42_mp,
	   a43_mp,a44_mp,a51_mp,a52_mp,a53_mp,a54_mp);

	icomp = mpz_sgn(temp1);

/* if major determinant is non 0, return its sign */

	if (icomp != 0) {
		*res = icomp;
		return;
	}

/* Look now at each term in the expansion of the determinant with
   *respect to EPS                                                
	The initial determinant is:
	Minor5(i,j,k,l,m,1,2,3,4,0)			*/

/* Term 1: 	-Minor4(j,k,l,m,1,2,3,0) */

	deter4_gmp(temp1,a21_mp,a22_mp,a23_mp,a31_mp,a32_mp,a33_mp,
		a41_mp,a42_mp,a43_mp,a51_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res = - icomp;
		return;
	}

/* Term 2:	Minor4(j,k,l,m,1,2,4,0) */
	
	deter4_gmp(temp1,a21_mp,a22_mp,a24_mp,a31_mp,a32_mp,a34_mp,
		a41_mp,a42_mp,a44_mp,a51_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 3:	- Minor4(j,k,l,m,1,3,4,0) */
	
	deter4_gmp(temp1,a21_mp,a23_mp,a24_mp,a31_mp,a33_mp,a34_mp,
		a41_mp,a43_mp,a44_mp,a51_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  - icomp;
		return;
	}

/* Term 4:	Minor4(j,k,l,m,2,3,4,0) */
	
	deter4_gmp(temp1,a22_mp,a23_mp,a24_mp,a32_mp,a33_mp,a34_mp,
		a42_mp,a43_mp,a44_mp,a52_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 5:	Minor4(i,k,l,m,1,2,3,0) */
	
	deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a31_mp,a32_mp,a33_mp,
		a41_mp,a42_mp,a43_mp,a51_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 6:	Minor3(k,l,m,1,2,0) */
	
	deter3_gmp(temp1,a31_mp,a32_mp,a41_mp,a42_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 7:	-Minor3(k,l,m,1,3,0) */
	
	deter3_gmp(temp1,a31_mp,a33_mp,a41_mp,a43_mp,a51_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 8:	Minor3(k,l,m,2,3,0) */
	
	deter3_gmp(temp1,a32_mp,a33_mp,a42_mp,a43_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 9:	-Minor4(i,k,l,m,1,2,4,0) */
	
	deter4_gmp(temp1,a11_mp,a12_mp,a14_mp,a31_mp,a32_mp,a34_mp,
		a41_mp,a42_mp,a44_mp,a51_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 10:	Minor3(k,l,m,1,4,0) */
	
	deter3_gmp(temp1,a31_mp,a34_mp,a41_mp,a44_mp,a51_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 11:	-Minor3(k,l,m,2,4,0) */
	
	deter3_gmp(temp1,a32_mp,a34_mp,a42_mp,a44_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 12:	Minor4(i,k,l,m,1,3,4,0) */
	
	deter4_gmp(temp1,a11_mp,a13_mp,a14_mp,a31_mp,a33_mp,a34_mp,
		a41_mp,a43_mp,a44_mp,a51_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 13:	Minor3(k,l,m,3,4,0) */
	
	deter3_gmp(temp1,a33_mp,a34_mp,a43_mp,a44_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 14:	-Minor4(i,k,l,m,2,3,4,0) */
	
	deter4_gmp(temp1,a12_mp,a13_mp,a14_mp,a32_mp,a33_mp,a34_mp,
		a42_mp,a43_mp,a44_mp,a52_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 15:	-Minor4(i,j,l,m,1,2,3,0) */
	
	deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
		a41_mp,a42_mp,a43_mp,a51_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 16:	-Minor3(j,l,m,1,2,0) */
	
	deter3_gmp(temp1,a21_mp,a22_mp,a41_mp,a42_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 17:	Minor3(j,l,m,1,3,0) */
	
	deter3_gmp(temp1,a21_mp,a23_mp,a41_mp,a43_mp,a51_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 18:	-Minor3(j,l,m,2,3,0) */
	
	deter3_gmp(temp1,a22_mp,a23_mp,a42_mp,a43_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 19:	Minor3(i,l,m,1,2,0) */
	
	deter3_gmp(temp1,a11_mp,a12_mp,a41_mp,a42_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 20:	-Minor2(l,m,1,0) */

	deter2_gmp(temp1,a41_mp,a51_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 21:	Minor2(l,m,2,0) */

	deter2_gmp(temp1,a42_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 22:	-Minor3(i,l,m,1,3,0) */
	
	deter3_gmp(temp1,a11_mp,a13_mp,a41_mp,a43_mp,a51_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 23:	-Minor2(l,m,3,0) */

	deter2_gmp(temp1,a43_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 24:	Minor3(i,l,m,2,3,0) */
	
	deter3_gmp(temp1,a12_mp,a13_mp,a42_mp,a43_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 25:	Minor4(i,j,l,m,1,2,4,0) */
	
	deter4_gmp(temp1,a11_mp,a12_mp,a14_mp,a21_mp,a22_mp,a24_mp,
		a41_mp,a42_mp,a44_mp,a51_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 26:	-Minor3(j,l,m,1,4,0) */
	
	deter3_gmp(temp1,a21_mp,a24_mp,a41_mp,a44_mp,a51_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 27:	Minor3(j,l,m,2,4,0) */
	
	deter3_gmp(temp1,a22_mp,a24_mp,a42_mp,a44_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 28:	Minor3(i,l,m,1,4,0) */
	
	deter3_gmp(temp1,a11_mp,a14_mp,a41_mp,a44_mp,a51_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 29:	Minor2(l,m,4,0) */

	deter2_gmp(temp1,a44_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 30:	-Minor3(i,l,m,2,4,0) */
	
	deter3_gmp(temp1,a12_mp,a14_mp,a42_mp,a44_mp,a52_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 31:	-Minor4(i,j,l,m,1,3,4,0) */
	
	deter4_gmp(temp1,a11_mp,a13_mp,a14_mp,a21_mp,a23_mp,a24_mp,
		a41_mp,a43_mp,a44_mp,a51_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 32:	-Minor3(j,l,m,3,4,0) */
	
	deter3_gmp(temp1,a23_mp,a24_mp,a43_mp,a44_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 33:	Minor3(i,l,m,3,4,0) */
	
	deter3_gmp(temp1,a13_mp,a14_mp,a43_mp,a44_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 34:	Minor4(i,j,l,m,2,3,4,0) */
	
	deter4_gmp(temp1,a12_mp,a13_mp,a14_mp,a22_mp,a23_mp,a24_mp,
		a42_mp,a43_mp,a44_mp,a52_mp,a53_mp,a54_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 35:	Minor4(i,j,k,m,1,2,3,0) */
	
	deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
		a31_mp,a32_mp,a33_mp,a51_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 36:	Minor3(j,k,m,1,2,0) */
	
	deter3_gmp(temp1,a21_mp,a22_mp,a31_mp,a32_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 37:	-Minor3(j,k,m,1,3,0) */
	
	deter3_gmp(temp1,a21_mp,a23_mp,a31_mp,a33_mp,a51_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 38:	Minor3(j,k,m,2,3,0) */
	
	deter3_gmp(temp1,a22_mp,a23_mp,a32_mp,a33_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 39:	-Minor3(i,k,m,1,2,0) */
	
	deter3_gmp(temp1,a11_mp,a12_mp,a31_mp,a32_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 40:	Minor2(k,m,1,0) */

	deter2_gmp(temp1,a31_mp,a51_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 41:	-Minor2(k,m,2,0) */

	deter2_gmp(temp1,a32_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 42:	Minor3(i,k,m,1,3,0) */
	
	deter3_gmp(temp1,a11_mp,a13_mp,a31_mp,a33_mp,a51_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 43:	Minor2(k,m,3,0) */

	deter2_gmp(temp1,a33_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 44:	-Minor3(i,k,m,2,3,0) */
	
	deter3_gmp(temp1,a12_mp,a13_mp,a32_mp,a33_mp,a52_mp,a53_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 45:	Minor3(i,j,m,1,2,0) */
	
	deter3_gmp(temp1,a11_mp,a12_mp,a21_mp,a22_mp,a51_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 46:	-Minor2(j,m,1,0) */

	deter2_gmp(temp1,a21_mp,a51_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  -icomp;
		return;
	}

/* Term 47:	Minor2(j,m,2,0) */

	deter2_gmp(temp1,a22_mp,a52_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 48:	Minor2(i,m,1,0) */

	deter2_gmp(temp1,a11_mp,a51_mp);
	icomp = mpz_sgn(temp1);
	if (icomp != 0) {
		*res =  icomp;
		return;
	}

/* Term 49:	1 */

	*res = 1;

  }

/* ==============================================================================
   minor4_:								  
	This subroutine tests the sign of the determinant
		D = | a11 a12 a13 1 |
		    | a21 a22 a23 1 |
		    | a31 a32 a33 1 |
		    | a41 a42 a43 1 |
	and return 1 if positive, -1 if negative, 0 otherwise 
   ============================================================================== */

  void SOS::minor4_gmp(double *coord_a, double *coord_b, double *coord_c, double *coord_d, int *res)
  {
	int icomp;

/* Transfer coordinates to gmp */

	real_to_gmp(coord_a[0],a11_mp);
	real_to_gmp(coord_a[1],a12_mp);
	real_to_gmp(coord_a[2],a13_mp);
	real_to_gmp(coord_b[0],a21_mp);
	real_to_gmp(coord_b[1],a22_mp);
	real_to_gmp(coord_b[2],a23_mp);
	real_to_gmp(coord_c[0],a31_mp);
	real_to_gmp(coord_c[1],a32_mp);
	real_to_gmp(coord_c[2],a33_mp);
	real_to_gmp(coord_d[0],a41_mp);
	real_to_gmp(coord_d[1],a42_mp);
	real_to_gmp(coord_d[2],a43_mp);

/* Compute determinant */

	deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
			   a31_mp,a32_mp,a33_mp,a41_mp,a42_mp,a43_mp);

	icomp = mpz_sgn(temp1);

	*res = icomp;

  }

#endif
