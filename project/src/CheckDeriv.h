
/* ===============================================================================================
   Compares analytical and numerical derivatives of the volumes of the union of balls
   =============================================================================================== */

   void CheckDeriv(int natoms, double *coord, double *radii, double *coefS, double *coefV, 
	double *coefM, double *coefG, double *dsurf, double *dvol, double *dmean, double *dgauss)
   {

	double dx = 1.e-4;
	double u1, u2, v1, v2, w1, w2, x1, x2;
	double alpha = 0;
	double err1 = 0;
	double err2 = 0;
	double err3 = 0;
	double err4 = 0;
	int flag_deriv = 0;

	std::vector<Vertex> vertices;
	std::vector<Tetrahedron> tetra;
	std::vector<Edge> edges;
	std::vector<Face> faces;

	double Surf, WSurf, Vol, WVol, Mean, WMean, Gauss, WGauss;

	int nfudge = 8;
	double *ballwsurf = new double[natoms+nfudge];
	double *dsurf_num = new double[3*(natoms+nfudge)];
	memset(dsurf_num, 0, 3*(natoms+nfudge)*sizeof(double));

	double *ballwvol;
	double *dvol_num;
	ballwvol = new double[natoms+nfudge];
	dvol_num = new double[3*(natoms+nfudge)];
	memset(dvol_num, 0, 3*(natoms+nfudge)*sizeof(double));

	double *ballwmean;
	double *dmean_num;
	ballwmean = new double[natoms+nfudge];
	dmean_num = new double[3*(natoms+nfudge)];
	memset(dmean_num, 0, 3*(natoms+nfudge)*sizeof(double));

	double *ballwgauss;
	double *dgauss_num;
	ballwgauss = new double[natoms+nfudge];
	dgauss_num = new double[3*(natoms+nfudge)];
	memset(dgauss_num, 0, 3*(natoms+nfudge)*sizeof(double));

	for(int i = 0; i < 3*natoms; i++) {

		coord[i] += dx;
		delcx.setup(natoms, coord, radii, coefS, coefV, coefM, coefG, vertices, tetra);
		delcx.regular3D(vertices, tetra);
		alfcx.alfcx(alpha, vertices, tetra);
		alfcx.alphacxEdges(tetra, edges);
		alfcx.alphacxFaces(tetra, faces);
		volumes.ball_dvolumes(vertices, tetra, edges, faces, &WSurf, &WVol,
		&WMean, &WGauss, &Surf, &Vol, &Mean, &Gauss, ballwsurf, ballwvol,
		ballwmean, ballwgauss, dsurf, dvol, dmean, dgauss, flag_deriv);
		u1 = WSurf;
		v1 = WVol;
		w1 = WMean;
		x1 = WGauss;
		coord[i] -= 2*dx;
		delcx.setup(natoms, coord, radii, coefS, coefV, coefM, coefG, vertices, tetra);
		delcx.regular3D(vertices, tetra);
		alfcx.alfcx(alpha, vertices, tetra);
		alfcx.alphacxEdges(tetra, edges);
		alfcx.alphacxFaces(tetra, faces);
		volumes.ball_dvolumes(vertices, tetra, edges, faces, &WSurf, &WVol,
		&WMean, &WGauss, &Surf, &Vol, &Mean, &Gauss, ballwsurf, ballwvol,
		ballwmean, ballwgauss, dsurf, dvol, dmean, dgauss, flag_deriv);
		u2 = WSurf;
		v2 = WVol;
		w2 = WMean;
		x2 = WGauss;
		dsurf_num[i] = (u1-u2)/(2*dx);
		err1 += (dsurf[i] - dsurf_num[i])*(dsurf[i] - dsurf_num[i]);
		dvol_num[i] = (v1-v2)/(2*dx);
		err2 += (dvol[i] - dvol_num[i])*(dvol[i] - dvol_num[i]);
		dmean_num[i] = (w1-w2)/(2*dx);
		err3 += (dmean[i] - dmean_num[i])*(dmean[i] - dmean_num[i]);
		dgauss_num[i] = (x1-x2)/(2*dx);
		err4 += (dgauss[i] - dgauss_num[i])*(dgauss[i] - dgauss_num[i]);

		coord[i] += dx;

	}

	err1 = std::sqrt(err1/(3*natoms)); 
	err2 = std::sqrt(err2/(3*natoms)); 
	err3 = std::sqrt(err3/(3*natoms)); 
	err4 = std::sqrt(err4/(3*natoms)); 

	std::cout<< "RMS error between analytical and numerical surface derivs    : " << err1 << std::endl;
	std::cout<< "RMS error between analytical and numerical volume derivs     : " << err2 << std::endl;
	std::cout<< "RMS error between analytical and numerical mean curv derivs  : " << err3 << std::endl;
	std::cout<< "RMS error between analytical and numerical Gauss curv derivs : " << err4 << std::endl;
	std::cout << " " << std::endl;


	std::cout << "Individual errors on surface derivatives: " << std::endl;
	for(int i = 0; i < natoms; i++) {
		std::cout << "Atom: " << i << " Coordinate: x dS anal: " << dsurf[3*i] << " dS num = " << dsurf_num[3*i] << std::endl;
		std::cout << "Atom: " << i << " Coordinate: y dS anal: " << dsurf[3*i+1] << " dS num = " << dsurf_num[3*i+1] << std::endl;
		std::cout << "Atom: " << i << " Coordinate: z dS anal: " << dsurf[3*i+2] << " dS num = " << dsurf_num[3*i+2] << std::endl;
	}
	std::cout << " " << std::endl;

	std::cout << "Individual errors on volume derivatives: " << std::endl;
	for(int i = 0; i < natoms; i++) {
		std::cout << "Atom: " << i << " Coordinate: x dV anal: " << dvol[3*i] << " dV num = " << dvol_num[3*i] << std::endl;
		std::cout << "Atom: " << i << " Coordinate: y dV anal: " << dvol[3*i+1] << " dV num = " << dvol_num[3*i+1] << std::endl;
		std::cout << "Atom: " << i << " Coordinate: z dV anal: " << dvol[3*i+2] << " dV num = " << dvol_num[3*i+2] << std::endl;
	}
	std::cout << " " << std::endl;

	std::cout << "Individual errors on mean curvature derivatives: " << std::endl;
	for(int i = 0; i < natoms; i++) {
		std::cout << "Atom: " << i << " Coordinate: x dM anal: " << dmean[3*i] << " dM num = " << dmean_num[3*i] << std::endl;
		std::cout << "Atom: " << i << " Coordinate: y dM anal: " << dmean[3*i+1] << " dM num = " << dmean_num[3*i+1] << std::endl;
		std::cout << "Atom: " << i << " Coordinate: z dM anal: " << dmean[3*i+2] << " dM num = " << dmean_num[3*i+2] << std::endl;
	}
	std::cout << " " << std::endl;

	std::cout << "Individual errors on Gauss curvature derivatives: " << std::endl;
	for(int i = 0; i < natoms; i++) {
		std::cout << "Atom: " << i << " Coordinate: x dG anal: " << dgauss[3*i] << " dG num = " << dgauss_num[3*i] << std::endl;
		std::cout << "Atom: " << i << " Coordinate: y dG anal: " << dgauss[3*i+1] << " dG num = " << dgauss_num[3*i+1] << std::endl;
		std::cout << "Atom: " << i << " Coordinate: z dG anal: " << dgauss[3*i+2] << " dG num = " << dgauss_num[3*i+2] << std::endl;
	}
	std::cout << " " << std::endl;


}
