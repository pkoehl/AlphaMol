/*********************************************************************************
 *	The Vertex class
 *********************************************************************************/

#ifndef VERTEX_H
#define VERTEX_H

  #include <gmp.h>
  #include <vector>

/*********************************************************************************
  Class that characterizes each vertex of the Delaunay/Alpha complex
 *********************************************************************************/

  class Vertex {
	public:
		double Radius;
		double Coordinates[3];
		double Weight;
		double CoefS, CoefV, CoefM, CoefG;
		double gamma;

		std::bitset<8> info;
		bool status;

		Vertex() {
		}

		Vertex(double x, double y, double z, double radius, double aspS,
			double aspV, double aspM, double aspG);
		~Vertex();

	private:
		double truncate_real(double x, int ndigit);
  };


/*********************************************************************************
        Constructor: set coordinates to specific components
 *********************************************************************************/

Vertex::Vertex(double x, double y, double z, double radius, double aspS, double aspV,
	double aspM, double aspG) {


	int ndigit = 8;
	double x1 = truncate_real(x, ndigit);
	double y1 = truncate_real(y, ndigit);
	double z1 = truncate_real(z, ndigit);
	double r1 = truncate_real(radius, ndigit);
	this->Coordinates[0] = x1;
	this->Coordinates[1] = y1;
	this->Coordinates[2] = z1;
	this->Radius = r1;
	this->CoefS = aspS;
	this->CoefV = aspV;
	this->CoefM = aspM;
	this->CoefG = aspG;

	std::bitset<8> b(std::string("00000000"));
	this->info = b;
	this->info[1] = 1;

	long long ival1= std::round(r1*10000.0);
	long long ival2 = -ival1*ival1;
	ival1 = std::round(x1*10000.0); 
	ival2 += ival1*ival1;
	ival1 = std::round(y1*10000.0); 
	ival2 += ival1*ival1;
	ival1 = std::round(z1*10000.0); 
	ival2 += ival1*ival1;
	this->Weight = (double) ival2/100000000.0;
	this->gamma = 0;

  }

/*********************************************************************************
        Destructor
 *********************************************************************************/

  Vertex::~Vertex() {
  }

/*********************************************************************************
        Truncate real
 *********************************************************************************/

 double Vertex::truncate_real(double x, int ndigit)
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

#endif // VERTEX_H
