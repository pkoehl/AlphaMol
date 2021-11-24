/* ====================================================================
	sets of tools to compute the contributions of the corners
	of the space-filling diagram to the Gaussian curvature
 ==================================================================== */

#ifndef GAUSSCORNER_H
#define GAUSSCORNER_H

/* ====================================================================
   class
 ==================================================================== */

  class GAUSSCORNER {

	public:

		void threesphere_dgauss(double ra, double rb, double rc, 
		double ra2, double rb2, double rc2,
		double rab, double rac, double rbc,
		double rab2, double rac2, double rbc2,
		double *areaA, double *areaB, double *areaC, double darea[3][3], int option);

	private:

		double trig_darea(double a, double b, double c, double *der_S, int option);

		double trig_dradius(double a, double b, double c, double *der_r, int option);

		double sign(double a, double b, double c);

  };

/* ====================================================================
   If R is the radius of the circumcircle to 3 points on the surface of the
   unit sphere, compute r = cos^2(R)

  The 3 points considered are A, B, C; phi_AB is the arc length on the surface
  of the unit sphere between A and B, with similar definitions for phi_AC
  and phi_BC

  Input:
	a = cos^2(phi_AB/2)
	b = cos^2(phi_AC/2)
	c = cos^2(phi_BC/2)
	option: flag: 0: do not compute derivatives, 1: compute derivatives
  Output:
	r = cos^2(R)
	dr: array containing dr/da, dr/db, dr/dc
		
 ==================================================================== */

  double GAUSSCORNER::trig_dradius(double a, double b, double c, double *der_r, int option)
  {

	double u = 4*a*b*c - (a+b+c-1)*(a+b+c-1);
	double sqr_u = std::sqrt(u);

	double v = (a-1)*(a-1) + (b-1)*(b-1) + (c-1)*(c-1) - (a-b)*(a-b) - (a-c)*(a-c) - (b-c)*(b-c);
	double sqr_v = std::sqrt(v);

	double r = 0.5 + 0.5* sqr_u/sqr_v;

	if(option==0) return r;

	double du_da = 4*b*c - 2*(a+b+c-1);
	double du_db = 4*a*c - 2*(a+b+c-1);
	double du_dc = 4*a*b - 2*(a+b+c-1);

	double dv_da = 2*(b+c-a-1);
	double dv_db = 2*(a+c-b-1);
	double dv_dc = 2*(a+b-c-1);

	der_r[0] = 0.5*(r-0.5)*(du_da/u -dv_da/v);
	der_r[1] = 0.5*(r-0.5)*(du_db/u -dv_db/v);
	der_r[2] = 0.5*(r-0.5)*(du_dc/u -dv_dc/v);

	return r;
  }

/* ====================================================================
   Computes the surface area S of a spherical triangle and its derivatives

  The 3 points considered are A, B, C; phi_AB is the arc length on the surface
  of the unit sphere between A and B, with similar definitions for phi_AC
  and phi_BC

  Input:
	a = cos^2(phi_AB/2)
	b = cos^2(phi_AC/2)
	c = cos^2(phi_BC/2)
	option: flag: 0: do not compute derivatives, 1: compute derivatives
  Output:
	S: surface area
	dS: array containing dS/da, dS/db, dS/dc
		
 ==================================================================== */

   double GAUSSCORNER::trig_darea(double a, double b, double c, double *der_S, int option)
   {

	double tol = 1.e-14;

	double u = 4*a*b*c - (a+b+c-1)*(a+b+c-1);
	double v = 4*a*b*c;
	if(std::abs(u) < tol) u = 0.0;
	double w = std::sqrt(std::abs(u));
	double t = std::sqrt(v);

	double S = 2*std::asin(w/t);

	if(option==0) return S;

	if(w>0) {
		der_S[0] = (b+c-a-1)/(a*w);
		der_S[1] = (a+c-b-1)/(b*w);
		der_S[2] = (a+b-c-1)/(c*w);
	} else {
		der_S[0] = 0;
		der_S[1] = 0;
		der_S[2] = 0;
	}

	return S;

   }

/* ====================================================================
   Sign function: defines if the surface area of a spherical triangle
   is positive or negative
 ==================================================================== */

  double GAUSSCORNER::sign(double a, double b, double c)
  {

	double s = -1;
	if( a + b <= 1 + c) {
		s = 1;
	} else {
		s = -1.;
	}

	return s;

  }

/* ====================================================================
   Computes the relative contribution of the 3 vertices of a triangle
   in the alpha complex to the Gaussian curvature
 ==================================================================== */

  void GAUSSCORNER::threesphere_dgauss(double ra, double rb, double rc, 
	double ra2, double rb2, double rc2, double rab, double rac, double rbc,
	double rab2, double rac2, double rbc2, double *areaA, double *areaB, 
	double *areaC, double darea[3][3], int option)
  {

	double cos_ab, cos_ac, cos_bc;
	double da_ab, db_bc, dc_ac;
	double a, b, c;
	double r;
	double sign_a, sign_b, sign_c;
	double S, Sa, Sb, Sc;
	double der_r[3], der_S[3], der_Sa[3], der_Sb[3], der_Sc[3];
	double dSa[3], dSb[3], dSc[3];
	double dSa_dab, dSa_dac, dSa_dbc;
	double dSb_dab, dSb_dac, dSb_dbc;
	double dSc_dab, dSc_dac, dSc_dbc;

	cos_ab = (ra2+rb2-rab2)/(2.0*ra*rb);
	a = 0.5*(cos_ab+1);

	cos_bc = (rb2+rc2-rbc2)/(2.0*rb*rc);
	b = 0.5*(cos_bc+1);

	cos_ac = (ra2+rc2-rac2)/(2.0*ra*rc);
	c = 0.5*(cos_ac+1);

	sign_a = sign(a, c, b);
	sign_b = sign(a, b, c);
	sign_c = sign(c, b, a);

	der_r[0] = 0; der_r[1] = 0; der_r[2] = 0;
	r = trig_dradius(a, b, c, der_r, option);

	der_S[0] = 0; der_S[1] = 0; der_S[2] = 0;
	der_Sa[0] = 0; der_Sa[1] = 0; der_Sa[2] = 0;
	der_Sb[0] = 0; der_Sb[1] = 0; der_Sb[2] = 0;
	der_Sc[0] = 0; der_Sc[1] = 0; der_Sc[2] = 0;

	S  = trig_darea(a, b, c, der_S, option);
	Sa = trig_darea(a, r, r, der_Sa, option);
	Sb = trig_darea(r, b, r, der_Sb, option);
	Sc = trig_darea(r, r, c, der_Sc, option);

	dSa[0] = der_Sa[0] + (der_Sa[1]+der_Sa[2])*der_r[0];
	dSa[1] = (der_Sa[1]+der_Sa[2])*der_r[1];
	dSa[2] = (der_Sa[1]+der_Sa[2])*der_r[2];

	dSb[0] = (der_Sb[0]+der_Sb[2])*der_r[0];
	dSb[1] = der_Sb[1] + (der_Sb[0]+der_Sb[2])*der_r[1];
	dSb[2] = (der_Sb[0]+der_Sb[2])*der_r[2];

	dSc[0] = (der_Sc[0]+der_Sc[1])*der_r[0];
	dSc[1] = (der_Sc[0]+der_Sc[1])*der_r[1];
	dSc[2] = der_Sc[2] + (der_Sc[0]+der_Sc[1])*der_r[2];

	if(Sa==0) {
		Sa = S - sign_a*Sb - sign_b*Sc;
		dSa[0] = der_S[0] - sign_a*dSb[0] - sign_b*dSc[0];
		dSa[1] = der_S[1] - sign_a*dSb[1] - sign_b*dSc[1];
		dSa[2] = der_S[2] - sign_a*dSb[2] - sign_b*dSc[2];
	}
	if(Sb==0) {
		Sb = S - sign_c*Sa - sign_b*Sc;
		dSb[0] = der_S[0] - sign_c*dSa[0] - sign_b*dSc[0];
		dSb[1] = der_S[1] - sign_c*dSa[1] - sign_b*dSc[1];
		dSb[2] = der_S[2] - sign_c*dSa[2] - sign_b*dSc[2];
	}
	if(Sc==0) {
		Sc = S - sign_c*Sa - sign_a*Sb;
		dSc[0] = der_S[0] - sign_c*dSa[0] - sign_a*dSb[0];
		dSc[1] = der_S[1] - sign_c*dSa[1] - sign_a*dSb[1];
		dSc[2] = der_S[2] - sign_c*dSa[2] - sign_a*dSb[2];
	}

	*areaA = 0.5*(sign_c*Sa + sign_b*Sc);
	*areaB = 0.5*(sign_a*Sb + sign_c*Sa);
	*areaC = 0.5*(sign_b*Sc + sign_a*Sb);

	if(option==0) return;

	da_ab = -0.5*rab/(ra*rb);
	db_bc = -0.5*rbc/(rb*rc);
	dc_ac = -0.5*rac/(ra*rc);

	dSa_dab = dSa[0]*da_ab;
	dSa_dbc = dSa[1]*db_bc;
	dSa_dac = dSa[2]*dc_ac;

	dSb_dab = dSb[0]*da_ab;
	dSb_dbc = dSb[1]*db_bc;
	dSb_dac = dSb[2]*dc_ac;

	dSc_dab = dSc[0]*da_ab;
	dSc_dbc = dSc[1]*db_bc;
	dSc_dac = dSc[2]*dc_ac;

	darea[0][0] = 0.5*(sign_c*dSa_dab + sign_b*dSc_dab);
	darea[0][1] = 0.5*(sign_c*dSa_dac + sign_b*dSc_dac);
	darea[0][2] = 0.5*(sign_c*dSa_dbc + sign_b*dSc_dbc);

	darea[1][0] = 0.5*(sign_a*dSb_dab + sign_c*dSa_dab);
	darea[1][1] = 0.5*(sign_a*dSb_dac + sign_c*dSa_dac);
	darea[1][2] = 0.5*(sign_a*dSb_dbc + sign_c*dSa_dbc);

	darea[2][0] = 0.5*(sign_b*dSc_dab + sign_a*dSb_dab);
	darea[2][1] = 0.5*(sign_b*dSc_dac + sign_a*dSb_dac);
	darea[2][2] = 0.5*(sign_b*dSc_dbc + sign_a*dSb_dbc);

  }

#endif
