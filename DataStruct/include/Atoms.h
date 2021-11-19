/*********************************************************************************
 *	The Atoms class
 *********************************************************************************/

#ifndef ATOMS_H
#define ATOMS_H

  #include <gmp.h>
  #include <vector>

/*********************************************************************************
  Class that characterizes each vertex of the Delaunay/Alpha complex
 *********************************************************************************/

  class Atoms {
	public:
		double Radius;
		double Coordinates[3];
		double Weight;
		double CoefS, CoefV, CoefM, CoefG;

		Atoms() {
		}

		Atoms(double x, double y, double z, double radius, double aspS, double aspV, double aspM, double aspG);
		~Atoms();

	private:
		double truncate_real(double x, int ndigit);
  };


/*********************************************************************************
        Constructor: set coordinates to specific components
 *********************************************************************************/

Atoms::Atoms(double x, double y, double z, double radius, double aspS, double aspV, double aspM, double aspG) {


	int ndigit = 8;
	double x1 = truncate_real(x, ndigit);
	double y1 = truncate_real(y, ndigit);
	double z1 = truncate_real(z, ndigit);
	double r1 = truncate_real(radius, ndigit);
	this->Coordinates[0] = x1;
	this->Coordinates[1] = y1;
	this->Coordinates[2] = z1;
	this->CoefS = aspS;
	this->CoefV = aspV;
	this->CoefM = aspM;
	this->CoefG = aspG;
	this->Radius = r1;

	this->Weight = -r1*r1 + x1*x1 + y1*y1 + z1*z1;

  }

/*********************************************************************************
        Destructor
 *********************************************************************************/

  Atoms::~Atoms() {
  }

/*********************************************************************************
        Truncate real
 *********************************************************************************/

 double Atoms::truncate_real(double x, int ndigit)
 {
	double x_out, y;
	double fact;

	int mantissa;
	int digit;

	mantissa = (int) x;
	y = x - mantissa;

	x_out = mantissa;
	fact = 1;
	for(int i = 0; i < ndigit; i++) {
		fact *= 10;
		digit = (int) std::round(y*10);
		y = 10*(y-digit/10.0);
		x_out += digit/fact;
	}

	return x_out;
  }

#endif // ATOMS_H
