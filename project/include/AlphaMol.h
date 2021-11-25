/* ===============================================================================================
   AlphaMol: a program for computing geometric measures of a union of balls

   Author:  Patrice Koehl  (collaboration with Herbert Edelsbrunner
   Date:    9/22/2019
   Version: 1
   =============================================================================================== */

#ifndef _ALPHAMOL_H_
#define _ALPHAMOL_H_

/* ===============================================================================================
   System includes
   =============================================================================================== */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <cmath>
#include <ctime>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <cstdlib>
#include <limits>
#include <bitset>
#include <stack>
#include "jlcxx/jlcxx.hpp"

/* ===============================================================================================
   Local includes
   =============================================================================================== */

#include "Atoms.h"
#include "Vertex.h"
#include "Tetrahedron.h"
#include "Edge.h"
#include "Face.h"

#include "ReadInput.h"

#include "delcx.h"
#include "alfcx.h"
#include "volumes.h"

DELCX delcx;
ALFCX alfcx;
VOLUMES volumes;

/* ===============================================================================================
   Prototypes
   =============================================================================================== */

static void usage(char** argv);
bool parse_args(int argc, char **argv, std::string *INfile, int *flag_CA, 
	double *r_h2o, int *flag_deriv, std::string *OUTfile);

void calculate_measures(jlcxx::ArrayRef<double> outs, const jlcxx::ArrayRef<double> in_coordinates, const jlcxx::ArrayRef<double> in_radii, const double in_coefS, const double in_coefV, const double in_coefM, const double in_coefG, int8_t flag_deriv, int8_t info_out_flag);


#endif
