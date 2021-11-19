/*********************************************************************************
 *	ReadInput:
 *	A set of procedures for reading different types of input files
 *********************************************************************************/

#ifndef READINPUT_H
#define READINPUT_H

  #include <vector>
  #include <cassert>
  #include <map>
  #include <stdio.h>
  #include "Atoms.h"

  std::map<std::string, double> opls_rad;
  std::map<std::string, double> opls_asp;
 
/* ===============================================================================================
   Class for input
   =============================================================================================== */

  class ReadInput {
	public:
		ReadInput();
		~ReadInput();

		void readFromCRD(std::string filename, std::vector<Atoms> & atoms, double rad_H2O);
		void readFromPQR(std::string filename, int flag, std::vector<Atoms> & atoms, double rad_H2O);
		void readFromPDB(std::string filename, int flag, std::vector<Atoms> & atoms, double rad_H2O);
		void centerMol(std::vector<Atoms>& vertexList);

	private:
		void setRadiusASP(std::string resname, std::string atmname, double *r, double *asp);
		void setOPLS();
  };

/* ===============================================================================================
   (Empty) constructor
   =============================================================================================== */

  ReadInput::ReadInput() {
	setOPLS();
  }

/* ===============================================================================================
   (Empty) Destructor
   =============================================================================================== */

  ReadInput::~ReadInput() {
  }

/* ===============================================================================================
   Read from PDB file for each atom selected:
	- coordinates x,y, and z
	- Assign radius, ASP based on standard geometry, using OPLS
   =============================================================================================== */

  void ReadInput::readFromPDB(std::string fileName, int flag, std::vector<Atoms>& atoms, double rad_H2O)
  {
	std::string line;

  	std::ifstream inFile;
	inFile.open(fileName);

	double x,y,z;
	double r,asp;

	std::string coord;
	std::string resname, atmname;
	while (getline(inFile, line)) // until reach the end of file 
	{
		if (line.substr(0,6) == "ATOM  ") // only read lines with ATOM at position 0
		{
			if(flag == 1) {
				atmname = line.substr(13,4);
				resname = line.substr(17,4);
				coord = line.substr(30,24);
				sscanf(coord.c_str(),"%8lf%8lf%8lf",&x,&y,&z);
				setRadiusASP(resname, atmname, &r, &asp);
				r += rad_H2O;
				Atoms atm(x, y, z, r, asp, asp, asp, asp);
				atoms.push_back(atm);
			}
			else {
				if(line.find("CA")!= std::string::npos) {
					atmname = line.substr(13,4);
					resname = line.substr(17,4);
					coord = line.substr(30,24);
					sscanf(coord.c_str(),"%8lf%8lf%8lf",&x,&y,&z);
					setRadiusASP(resname, atmname, &r, &asp);
					r += rad_H2O;
					Atoms atm(x, y, z, r, asp, asp, asp, asp);
					atoms.push_back(atm);
				}
			}
		}
	}

	inFile.close();

  }

/* ===============================================================================================
   Read from PQR file for each atom selected:
	- coordinates x,y, and z and radius r
	- Assign ASP based on standard geometry, using OPLS
   =============================================================================================== */

  void ReadInput::readFromPQR(std::string fileName, int flag, std::vector<Atoms>& atoms, double rad_H2O)
  {
	std::string line;

  	std::ifstream inFile;
	inFile.open(fileName);

	double x,y,z, r2;
	double r,asp;

	std::string coord;
	std::string resname, atmname;
	char temp[10], temp1[10], temp2[10];
	int p;
	double charge;
	while (getline(inFile, line)) // until reach the end of file 
	{
		if (line.substr(0,6) == "ATOM  ") // only read lines with ATOM at position 0
		{
			if(flag == 1) {
				sscanf(line.c_str(), "%s %d %s %s %d %lf %lf %lf %lf %lf", 
				temp, &p, temp1, temp2, &p, &x, &y, &z, &charge, &r);
				r += rad_H2O;
				atmname=temp1;
				resname=temp2;
				setRadiusASP(resname, atmname, &r2, &asp);
				Atoms atm(x, y, z, r, asp, asp, asp, asp);
				atoms.push_back(atm);
			}
			else {
				if(line.find("CA")!= std::string::npos) {
					sscanf(line.c_str(), "%s %d %s %s %d %lf %lf %lf %lf %lf", 
					temp, &p, temp1, temp2, &p, &x, &y, &z, &charge, &r);
					r += rad_H2O;
					atmname=temp1;
					resname=temp2;
					setRadiusASP(resname, atmname, &r2, &asp);
					Atoms atm(x, y, z, r, asp, asp, asp, asp);
					atoms.push_back(atm);
				}
			}
		}
	}

	inFile.close();

  }

/* ===============================================================================================
   Read from CRD file for each atom selected:
	- coordinates x,y, and z, and radius
	- Assign ASP as 1
   =============================================================================================== */

  void ReadInput::readFromCRD(std::string fileName, std::vector<Atoms>& atoms, double rad_H2O)
  {
	std::string line;

  	std::ifstream inFile;
	inFile.open(fileName);

	double x,y,z;
	double r;
	double asp = 1.0;

	while (getline(inFile, line)) // until reach the end of file 
	{
		if(line.find("#") == std::string::npos) {
			sscanf(line.c_str(), "%lf %lf %lf %lf", &x, &y, &z, &r);
			r += rad_H2O;
			Atoms atm(x, y, z, r, asp, asp, asp, asp);
			atoms.push_back(atm);
		}
	}

	inFile.close();

  }

/* ===============================================================================================
   Center molecule on 0, 0, 0
   =============================================================================================== */

   void ReadInput::centerMol(std::vector<Atoms>& atoms)
   {

	int natoms = atoms.size();

	double c[3];
	c[0] = 0.0; c[1] = 0.0; c[2] = 0.0;
	for (int i = 0; i < natoms; i++)
	{
		c[0] += atoms[i].Coordinates[0];
		c[1] += atoms[i].Coordinates[1];
		c[2] += atoms[i].Coordinates[2];
	}
	for(int i = 0; i < 3; i++) c[i] /= natoms;

	for (int i = 0; i < natoms; i++)
	{
		atoms[i].Coordinates[0] -= c[0];
		atoms[i].Coordinates[1] -= c[1];
		atoms[i].Coordinates[2] -= c[2];
	}
  }
	
/* ===============================================================================================
   Get radius and ASP for a given atom
   =============================================================================================== */

  void ReadInput::setRadiusASP(std::string resname, std::string atmname, double *r, double *asp)
  {

	std::string info = resname.substr(0,4)+atmname.substr(0,4);
	*r = opls_rad[info];
	*asp = opls_asp[info];
  }

/* ===============================================================================================
   initialize OPLS
   =============================================================================================== */

  void ReadInput::setOPLS()
  {
	opls_rad["ALA N   "] =    1.62;
	opls_rad["ALA H   "] =    1.20;
	opls_rad["ALA CA  "] =    1.90;
	opls_rad["ALA CB  "] =    1.96;
	opls_rad["ALA C   "] =    1.88;
	opls_rad["ALA O   "] =    1.48;
	opls_rad["ALA OXT "] =    1.48;
	opls_rad["ARG N   "] =    1.62;
	opls_rad["ARG H   "] =    1.20;
	opls_rad["ARG CA  "] =    1.90;
	opls_rad["ARG CB  "] =    1.95;
	opls_rad["ARG CG  "] =    1.95;
	opls_rad["ARG CD  "] =    1.95;
	opls_rad["ARG NE  "] =    1.62;
	opls_rad["ARG HE  "] =    1.20;
	opls_rad["ARG CZ  "] =    1.12;
	opls_rad["ARG NH1 "] =    1.62;
	opls_rad["ARG HH11"] =    1.20;
	opls_rad["ARG HH12"] =    1.20;
	opls_rad["ARG NH2 "] =    1.62;
	opls_rad["ARG HH21"] =    1.20;
	opls_rad["ARG HH22"] =    1.20;
	opls_rad["ARG C   "] =    1.88;
	opls_rad["ARG O   "] =    1.48;
	opls_rad["ARG OXT "] =    1.48;
	opls_rad["ASN N   "] =    1.62;
	opls_rad["ASN H   "] =    1.20;
	opls_rad["ASN CA  "] =    1.90;
	opls_rad["ASN CB  "] =    1.95;
	opls_rad["ASN CG  "] =    1.88;
	opls_rad["ASN OD1 "] =    1.48;
	opls_rad["ASN ND2 "] =    1.62;
	opls_rad["ASN HD21"] =    1.20;
	opls_rad["ASN HD22"] =    1.20;
	opls_rad["ASN C   "] =    1.88;
	opls_rad["ASN O   "] =    1.48;
	opls_rad["ASP N   "] =    1.62;
	opls_rad["ASP H   "] =    1.20;
	opls_rad["ASP CA  "] =    1.90;
	opls_rad["ASP CB  "] =    1.95;
	opls_rad["ASP CG  "] =    1.68;
	opls_rad["ASP OD1 "] =    1.48;
	opls_rad["ASP OD2 "] =    1.48;
	opls_rad["ASP HD  "] =    0.00;
	opls_rad["ASP C   "] =    1.88;
	opls_rad["ASP O   "] =    1.48;
	opls_rad["ASP OXT "] =    1.48;
	opls_rad["CYS N   "] =    1.62;
	opls_rad["CYS H   "] =    1.20;
	opls_rad["CYS CA  "] =    1.90;
	opls_rad["CYS CB  "] =    1.95;
	opls_rad["CYS SG  "] =    1.77;
	opls_rad["CYS HG  "] =    1.20;
	opls_rad["CYS C   "] =    1.88;
	opls_rad["CYS O   "] =    1.48;
	opls_rad["CYS OXT "] =    1.48;
	opls_rad["CSS N   "] =    1.62;
	opls_rad["CSS H   "] =    1.20;
	opls_rad["CSS CA  "] =    1.90;
	opls_rad["CSS CB  "] =    1.95;
	opls_rad["CSS SG  "] =    1.77;
	opls_rad["CSS C   "] =    1.88;
	opls_rad["CSS O   "] =    1.48;
	opls_rad["GLN N   "] =    1.62;
	opls_rad["GLN H   "] =    1.20;
	opls_rad["GLN CA  "] =    1.90;
	opls_rad["GLN CB  "] =    1.95;
	opls_rad["GLN CG  "] =    1.95;
	opls_rad["GLN CD  "] =    1.88;
	opls_rad["GLN OE1 "] =    1.48;
	opls_rad["GLN NE2 "] =    1.62;
	opls_rad["GLN HE21"] =    1.20;
	opls_rad["GLN HE22"] =    1.20;
	opls_rad["GLN C   "] =    1.88;
	opls_rad["GLN O   "] =    1.48;
	opls_rad["GLN OXT "] =    1.48;
	opls_rad["GLU N   "] =    1.62;
	opls_rad["GLU H   "] =    1.20;
	opls_rad["GLU CA  "] =    1.90;
	opls_rad["GLU CB  "] =    1.95;
	opls_rad["GLU CG  "] =    1.95;
	opls_rad["GLU CD  "] =    1.68;
	opls_rad["GLU OE1 "] =    1.48;
	opls_rad["GLU OE2 "] =    1.48;
	opls_rad["GLU HE  "] =    0.00;
	opls_rad["GLU C   "] =    1.88;
	opls_rad["GLU O   "] =    1.48;
	opls_rad["GLU OXT "] =    1.48;
	opls_rad["GLY N   "] =    1.62;
	opls_rad["GLY H   "] =    1.20;
	opls_rad["GLY CA  "] =    1.90;
	opls_rad["GLY C   "] =    1.88;
	opls_rad["GLY O   "] =    1.48;
	opls_rad["GLY OXT "] =    1.48;
	opls_rad["HIS N   "] =    1.62;
	opls_rad["HIS H   "] =    1.20;
	opls_rad["HIS CA  "] =    1.90;
	opls_rad["HIS CB  "] =    1.95;
	opls_rad["HIS CG  "] =    1.88;
	opls_rad["HIS ND1 "] =    1.62;
	opls_rad["HIS HD1 "] =    1.20;
	opls_rad["HIS CE1 "] =    1.88;
	opls_rad["HIS NE2 "] =    1.62;
	opls_rad["HIS HE2 "] =    0.00;
	opls_rad["HIS CD2 "] =    1.88;
	opls_rad["HIS C   "] =    1.88;
	opls_rad["HIS O   "] =    1.48;
	opls_rad["HIS OXT "] =    1.48;
	opls_rad["ILE N   "] =    1.62;
	opls_rad["ILE H   "] =    1.20;
	opls_rad["ILE CA  "] =    1.90;
	opls_rad["ILE CB  "] =    1.93;
	opls_rad["ILE CG1 "] =    1.95;
	opls_rad["ILE CG2 "] =    1.96;
	opls_rad["ILE CD1 "] =    1.95;
	opls_rad["ILE C   "] =    1.88;
	opls_rad["ILE O   "] =    1.48;
	opls_rad["ILE OXT "] =    1.48;
	opls_rad["LEU N   "] =    1.62;
	opls_rad["LEU H   "] =    1.20;
	opls_rad["LEU CA  "] =    1.90;
	opls_rad["LEU CB  "] =    1.95;
	opls_rad["LEU CG  "] =    1.93;
	opls_rad["LEU CD1 "] =    1.96;
	opls_rad["LEU CD2 "] =    1.96;
	opls_rad["LEU C   "] =    1.88;
	opls_rad["LEU O   "] =    1.48;
	opls_rad["LEU OXT "] =    1.48;
	opls_rad["LYS N   "] =    1.62;
	opls_rad["LYS H   "] =    1.20;
	opls_rad["LYS CA  "] =    1.90;
	opls_rad["LYS CB  "] =    1.95;
	opls_rad["LYS CG  "] =    1.95;
	opls_rad["LYS CD  "] =    1.95;
	opls_rad["LYS CE  "] =    1.95;
	opls_rad["LYS NZ  "] =    1.62;
	opls_rad["LYS HZ1 "] =    1.20;
	opls_rad["LYS HZ2 "] =    1.20;
	opls_rad["LYS HZ3 "] =    1.20;
	opls_rad["LYS 1HZ "] =    1.20;
	opls_rad["LYS 2HZ "] =    1.20;
	opls_rad["LYS 3HZ "] =    1.20;
	opls_rad["LYS C   "] =    1.88;
	opls_rad["LYS O   "] =    1.48;
	opls_rad["LYS OXT "] =    1.48;
	opls_rad["MET N   "] =    1.62;
	opls_rad["MET H   "] =    1.20;
	opls_rad["MET CA  "] =    1.90;
	opls_rad["MET CB  "] =    1.95;
	opls_rad["MET CG  "] =    1.90;
	opls_rad["MET SD  "] =    1.77;
	opls_rad["MET CE  "] =    1.90;
	opls_rad["MET C   "] =    1.88;
	opls_rad["MET O   "] =    1.48;
	opls_rad["MET OXT "] =    1.48;
	opls_rad["PHE N   "] =    1.62;
	opls_rad["PHE H   "] =    1.20;
	opls_rad["PHE CA  "] =    1.90;
	opls_rad["PHE CB  "] =    1.95;
	opls_rad["PHE CG  "] =    1.88;
	opls_rad["PHE CD1 "] =    1.88;
	opls_rad["PHE CD2 "] =    1.88;
	opls_rad["PHE CE1 "] =    1.88;
	opls_rad["PHE CE2 "] =    1.88;
	opls_rad["PHE CZ  "] =    1.88;
	opls_rad["PHE C   "] =    1.88;
	opls_rad["PHE O   "] =    1.48;
	opls_rad["PHE OXT "] =    1.48;
	opls_rad["PRO N   "] =    1.62;
	opls_rad["PRO CD  "] =    1.90;
	opls_rad["PRO CA  "] =    1.90;
	opls_rad["PRO CB  "] =    1.95;
	opls_rad["PRO CG  "] =    1.95;
	opls_rad["PRO C   "] =    1.88;
	opls_rad["PRO O   "] =    1.48;
	opls_rad["PRO OXT "] =    1.48;
	opls_rad["SER N   "] =    1.62;
	opls_rad["SER H   "] =    1.20;
	opls_rad["SER CA  "] =    1.90;
	opls_rad["SER CB  "] =    1.95;
	opls_rad["SER OG  "] =    1.53;
	opls_rad["SER HG  "] =    1.20;
	opls_rad["SER C   "] =    1.88;
	opls_rad["SER O   "] =    1.48;
	opls_rad["SER OXT "] =    1.48;
	opls_rad["THR N   "] =    1.62;
	opls_rad["THR H   "] =    1.20;
	opls_rad["THR CA  "] =    1.90;
	opls_rad["THR CB  "] =    1.93;
	opls_rad["THR CG2 "] =    1.96;
	opls_rad["THR OG1 "] =    1.53;
	opls_rad["THR HG1 "] =    1.20;
	opls_rad["THR C   "] =    1.88;
	opls_rad["THR O   "] =    1.48;
	opls_rad["THR OXT "] =    1.48;
	opls_rad["TRP N   "] =    1.62;
	opls_rad["TRP H   "] =    1.20;
	opls_rad["TRP CA  "] =    1.90;
	opls_rad["TRP CB  "] =    1.95;
	opls_rad["TRP CG  "] =    1.88;
	opls_rad["TRP CD1 "] =    1.88;
	opls_rad["TRP NE1 "] =    1.62;
	opls_rad["TRP HE1 "] =    1.20;
	opls_rad["TRP CE2 "] =    1.88;
	opls_rad["TRP CZ2 "] =    1.88;
	opls_rad["TRP CH2 "] =    1.88;
	opls_rad["TRP CZ3 "] =    1.88;
	opls_rad["TRP CE3 "] =    1.88;
	opls_rad["TRP CD2 "] =    1.88;
	opls_rad["TRP C   "] =    1.88;
	opls_rad["TRP O   "] =    1.48;
	opls_rad["TRP CT2 "] =    1.88;
	opls_rad["TRP CNT "] =    1.90;
	opls_rad["TRP ONT "] =    1.48;
	opls_rad["TRP OXT "] =    1.48;
	opls_rad["TYR N   "] =    1.62;
	opls_rad["TYR H   "] =    1.20;
	opls_rad["TYR CA  "] =    1.90;
	opls_rad["TYR CB  "] =    1.95;
	opls_rad["TYR CG  "] =    1.88;
	opls_rad["TYR CD1 "] =    1.88;
	opls_rad["TYR CE1 "] =    1.88;
	opls_rad["TYR CZ  "] =    1.88;
	opls_rad["TYR OH  "] =    1.53;
	opls_rad["TYR HH  "] =    1.20;
	opls_rad["TYR CE2 "] =    1.88;
	opls_rad["TYR CD2 "] =    1.88;
	opls_rad["TYR C   "] =    1.88;
	opls_rad["TYR O   "] =    1.48;
	opls_rad["TYR OXT "] =    1.48;
	opls_rad["VAL N   "] =    1.62;
	opls_rad["VAL H   "] =    1.20;
	opls_rad["VAL CA  "] =    1.90;
	opls_rad["VAL CB  "] =    1.93;
	opls_rad["VAL CG1 "] =    1.96;
	opls_rad["VAL CG2 "] =    1.96;
	opls_rad["VAL C   "] =    1.88;
	opls_rad["VAL O   "] =    1.48;
	opls_rad["VAL OXT "] =    1.48;
	opls_rad["HOH O   "] =    1.57;
	opls_rad["HOH H1  "] =    1.20;
	opls_rad["HOH H2  "] =    1.20;
	opls_rad["HOH 1HO "] =    1.20;
	opls_rad["HOH 2HO "] =    1.20;
	opls_rad["CO  C   "] =    2.49;
	opls_rad["CO  O   "] =    1.54;
	opls_rad["O2  O1  "] =    1.54;
	opls_rad["O2  O2  "] =    1.54;
	opls_rad["TERNN   "] =    1.62;
	opls_rad["TERNHT1 "] =    1.20;
	opls_rad["TERNHT2 "] =    1.20;
	opls_rad["TERNHT3 "] =    1.20;
	opls_rad["ASP 1H  "] =    1.20;
	opls_rad["ASP 2H  "] =    1.20;
	opls_rad["ASP 3H  "] =    1.20;
	opls_rad["GLN 1H  "] =    1.20;
	opls_rad["GLN 2H  "] =    1.20;
	opls_rad["GLN 3H  "] =    1.20;
	opls_rad["LYS 1H  "] =    1.20;
	opls_rad["LYS 2H  "] =    1.20;
	opls_rad["LYS 3H  "] =    1.20;
	opls_rad["TERNCA  "] =    1.90;
	opls_rad["TERCC   "] =    1.88;
	opls_rad["TERCCA  "] =    1.90;
	opls_rad["TERCO   "] =    1.48;
	opls_rad["TERCOXT "] =    1.48;
	opls_rad["TERCOCT1"] =    1.48;
	opls_rad["TERCOCT2"] =    1.48;
	opls_rad["TERCHXT "] =    0.00;

	opls_asp["ALA N   "] =    0.01;
	opls_asp["ALA H   "] =    0.00;
	opls_asp["ALA CA  "] =    0.04;
	opls_asp["ALA CB  "] =    0.04;
	opls_asp["ALA C   "] =    0.04;
	opls_asp["ALA O   "] =    0.01;
	opls_asp["ALA OXT "] =   -0.01;
	opls_asp["ARG N   "] =    0.01;
	opls_asp["ARG H   "] =    0.00;
	opls_asp["ARG CA  "] =    0.04;
	opls_asp["ARG CB  "] =    0.04;
	opls_asp["ARG CG  "] =    0.04;
	opls_asp["ARG CD  "] =    0.04;
	opls_asp["ARG NE  "] =    0.01;
	opls_asp["ARG HE  "] =    0.00;
	opls_asp["ARG CZ  "] =    0.04;
	opls_asp["ARG NH1 "] =   -0.05;
	opls_asp["ARG HH11"] =    0.00;
	opls_asp["ARG HH12"] =    0.00;
	opls_asp["ARG NH2 "] =   -0.05;
	opls_asp["ARG HH21"] =    0.00;
	opls_asp["ARG HH22"] =    0.00;
	opls_asp["ARG C   "] =    0.04;
	opls_asp["ARG O   "] =    0.01;
	opls_asp["ARG OXT "] =   -0.01;
	opls_asp["ASN N   "] =    0.01;
	opls_asp["ASN H   "] =    0.00;
	opls_asp["ASN CA  "] =    0.04;
	opls_asp["ASN CB  "] =    0.04;
	opls_asp["ASN CG  "] =    0.04;
	opls_asp["ASN OD1 "] =   -0.01;
	opls_asp["ASN ND2 "] =    0.01;
	opls_asp["ASN HD21"] =    0.00;
	opls_asp["ASN HD22"] =    0.00;
	opls_asp["ASN C   "] =    0.04;
	opls_asp["ASN O   "] =    0.01;
	opls_asp["ASP N   "] =    0.01;
	opls_asp["ASP H   "] =    0.00;
	opls_asp["ASP CA  "] =    0.04;
	opls_asp["ASP CB  "] =    0.04;
	opls_asp["ASP CG  "] =    0.04;
	opls_asp["ASP OD1 "] =   -0.01;
	opls_asp["ASP OD2 "] =   -0.01;
	opls_asp["ASP HD  "] =    0.00;
	opls_asp["ASP C   "] =    0.04;
	opls_asp["ASP O   "] =    0.01;
	opls_asp["ASP OXT "] =   -0.01;
	opls_asp["CYS N   "] =    0.01;
	opls_asp["CYS H   "] =    0.00;
	opls_asp["CYS CA  "] =    0.04;
	opls_asp["CYS CB  "] =    0.04;
	opls_asp["CYS SG  "] =    0.04;
	opls_asp["CYS HG  "] =    0.00;
	opls_asp["CYS C   "] =    0.04;
	opls_asp["CYS O   "] =    0.01;
	opls_asp["CYS OXT "] =   -0.01;
	opls_asp["CSS N   "] =    0.01;
	opls_asp["CSS H   "] =    0.00;
	opls_asp["CSS CA  "] =    0.04;
	opls_asp["CSS CB  "] =    0.04;
	opls_asp["CSS SG  "] =    0.04;
	opls_asp["CSS C   "] =    0.04;
	opls_asp["CSS O   "] =    0.01;
	opls_asp["GLN N   "] =    0.01;
	opls_asp["GLN H   "] =    0.00;
	opls_asp["GLN CA  "] =    0.04;
	opls_asp["GLN CB  "] =    0.04;
	opls_asp["GLN CG  "] =    0.04;
	opls_asp["GLN CD  "] =    0.04;
	opls_asp["GLN OE1 "] =   -0.01;
	opls_asp["GLN NE2 "] =    0.01;
	opls_asp["GLN HE21"] =    0.00;
	opls_asp["GLN HE22"] =    0.00;
	opls_asp["GLN C   "] =    0.04;
	opls_asp["GLN O   "] =    0.01;
	opls_asp["GLN OXT "] =   -0.01;
	opls_asp["GLU N   "] =    0.01;
	opls_asp["GLU H   "] =    0.00;
	opls_asp["GLU CA  "] =    0.04;
	opls_asp["GLU CB  "] =    0.04;
	opls_asp["GLU CG  "] =    0.04;
	opls_asp["GLU CD  "] =    0.04;
	opls_asp["GLU OE1 "] =   -0.01;
	opls_asp["GLU OE2 "] =   -0.01;
	opls_asp["GLU HE  "] =    0.00;
	opls_asp["GLU C   "] =    0.04;
	opls_asp["GLU O   "] =    0.01;
	opls_asp["GLU OXT "] =   -0.01;
	opls_asp["GLY N   "] =    0.01;
	opls_asp["GLY H   "] =    0.00;
	opls_asp["GLY CA  "] =    0.04;
	opls_asp["GLY C   "] =    0.04;
	opls_asp["GLY O   "] =    0.01;
	opls_asp["GLY OXT "] =   -0.01;
	opls_asp["HIS N   "] =    0.01;
	opls_asp["HIS H   "] =    0.00;
	opls_asp["HIS CA  "] =    0.04;
	opls_asp["HIS CB  "] =    0.04;
	opls_asp["HIS CG  "] =    0.04;
	opls_asp["HIS ND1 "] =    0.01;
	opls_asp["HIS HD1 "] =    0.00;
	opls_asp["HIS CE1 "] =    0.04;
	opls_asp["HIS NE2 "] =    0.01;
	opls_asp["HIS HE2 "] =    0.00;
	opls_asp["HIS CD2 "] =    0.04;
	opls_asp["HIS C   "] =    0.04;
	opls_asp["HIS O   "] =    0.01;
	opls_asp["HIS OXT "] =   -0.01;
	opls_asp["ILE N   "] =    0.01;
	opls_asp["ILE H   "] =    0.00;
	opls_asp["ILE CA  "] =    0.04;
	opls_asp["ILE CB  "] =    0.04;
	opls_asp["ILE CG1 "] =    0.04;
	opls_asp["ILE CG2 "] =    0.04;
	opls_asp["ILE CD1 "] =    0.04;
	opls_asp["ILE C   "] =    0.04;
	opls_asp["ILE O   "] =    0.01;
	opls_asp["ILE OXT "] =   -0.01;
	opls_asp["LEU N   "] =    0.01;
	opls_asp["LEU H   "] =    0.00;
	opls_asp["LEU CA  "] =    0.04;
	opls_asp["LEU CB  "] =    0.04;
	opls_asp["LEU CG  "] =    0.04;
	opls_asp["LEU CD1 "] =    0.04;
	opls_asp["LEU CD2 "] =    0.04;
	opls_asp["LEU C   "] =    0.04;
	opls_asp["LEU O   "] =    0.01;
	opls_asp["LEU OXT "] =   -0.01;
	opls_asp["LYS N   "] =    0.01;
	opls_asp["LYS H   "] =    0.00;
	opls_asp["LYS CA  "] =    0.04;
	opls_asp["LYS CB  "] =    0.04;
	opls_asp["LYS CG  "] =    0.04;
	opls_asp["LYS CD  "] =    0.04;
	opls_asp["LYS CE  "] =    0.04;
	opls_asp["LYS NZ  "] =   -0.05;
	opls_asp["LYS HZ1 "] =    0.00;
	opls_asp["LYS HZ2 "] =    0.00;
	opls_asp["LYS HZ3 "] =    0.00;
	opls_asp["LYS 1HZ "] =    0.00;
	opls_asp["LYS 2HZ "] =    0.00;
	opls_asp["LYS 3HZ "] =    0.00;
	opls_asp["LYS C   "] =    0.04;
	opls_asp["LYS O   "] =    0.01;
	opls_asp["LYS OXT "] =   -0.01;
	opls_asp["MET N   "] =    0.01;
	opls_asp["MET H   "] =    0.00;
	opls_asp["MET CA  "] =    0.04;
	opls_asp["MET CB  "] =    0.04;
	opls_asp["MET CG  "] =    0.04;
	opls_asp["MET SD  "] =    0.04;
	opls_asp["MET CE  "] =    0.04;
	opls_asp["MET C   "] =    0.04;
	opls_asp["MET O   "] =    0.01;
	opls_asp["MET OXT "] =   -0.01;
	opls_asp["PHE N   "] =    0.01;
	opls_asp["PHE H   "] =    0.00;
	opls_asp["PHE CA  "] =    0.04;
	opls_asp["PHE CB  "] =    0.04;
	opls_asp["PHE CG  "] =    0.04;
	opls_asp["PHE CD1 "] =    0.04;
	opls_asp["PHE CD2 "] =    0.04;
	opls_asp["PHE CE1 "] =    0.04;
	opls_asp["PHE CE2 "] =    0.04;
	opls_asp["PHE CZ  "] =    0.04;
	opls_asp["PHE C   "] =    0.04;
	opls_asp["PHE O   "] =    0.01;
	opls_asp["PHE OXT "] =   -0.01;
	opls_asp["PRO N   "] =    0.01;
	opls_asp["PRO CD  "] =    0.04;
	opls_asp["PRO CA  "] =    0.04;
	opls_asp["PRO CB  "] =    0.04;
	opls_asp["PRO CG  "] =    0.04;
	opls_asp["PRO C   "] =    0.04;
	opls_asp["PRO O   "] =    0.01;
	opls_asp["PRO OXT "] =   -0.01;
	opls_asp["SER N   "] =    0.01;
	opls_asp["SER H   "] =    0.00;
	opls_asp["SER CA  "] =    0.04;
	opls_asp["SER CB  "] =    0.04;
	opls_asp["SER OG  "] =    0.01;
	opls_asp["SER HG  "] =    0.00;
	opls_asp["SER C   "] =    0.04;
	opls_asp["SER O   "] =    0.01;
	opls_asp["SER OXT "] =   -0.01;
	opls_asp["THR N   "] =    0.01;
	opls_asp["THR H   "] =    0.00;
	opls_asp["THR CA  "] =    0.04;
	opls_asp["THR CB  "] =    0.04;
	opls_asp["THR CG2 "] =    0.04;
	opls_asp["THR OG1 "] =    0.01;
	opls_asp["THR HG1 "] =    0.00;
	opls_asp["THR C   "] =    0.04;
	opls_asp["THR O   "] =    0.01;
	opls_asp["THR OXT "] =   -0.01;
	opls_asp["TRP N   "] =    0.01;
	opls_asp["TRP H   "] =    0.00;
	opls_asp["TRP CA  "] =    0.04;
	opls_asp["TRP CB  "] =    0.04;
	opls_asp["TRP CG  "] =    0.04;
	opls_asp["TRP CD1 "] =    0.04;
	opls_asp["TRP NE1 "] =    0.01;
	opls_asp["TRP HE1 "] =    0.00;
	opls_asp["TRP CE2 "] =    0.04;
	opls_asp["TRP CZ2 "] =    0.04;
	opls_asp["TRP CH2 "] =    0.04;
	opls_asp["TRP CZ3 "] =    0.04;
	opls_asp["TRP CE3 "] =    0.04;
	opls_asp["TRP CD2 "] =    0.04;
	opls_asp["TRP C   "] =    0.04;
	opls_asp["TRP O   "] =    0.01;
	opls_asp["TRP CT2 "] =    0.04;
	opls_asp["TRP CNT "] =    0.04;
	opls_asp["TRP ONT "] =    0.01;
	opls_asp["TRP OXT "] =   -0.01;
	opls_asp["TYR N   "] =    0.01;
	opls_asp["TYR H   "] =    0.00;
	opls_asp["TYR CA  "] =    0.04;
	opls_asp["TYR CB  "] =    0.04;
	opls_asp["TYR CG  "] =    0.04;
	opls_asp["TYR CD1 "] =    0.04;
	opls_asp["TYR CE1 "] =    0.04;
	opls_asp["TYR CZ  "] =    0.04;
	opls_asp["TYR OH  "] =    0.01;
	opls_asp["TYR HH  "] =    0.00;
	opls_asp["TYR CE2 "] =    0.04;
	opls_asp["TYR CD2 "] =    0.04;
	opls_asp["TYR C   "] =    0.04;
	opls_asp["TYR O   "] =    0.01;
	opls_asp["TYR OXT "] =   -0.01;
	opls_asp["VAL N   "] =    0.01;
	opls_asp["VAL H   "] =    0.00;
	opls_asp["VAL CA  "] =    0.04;
	opls_asp["VAL CB  "] =    0.04;
	opls_asp["VAL CG1 "] =    0.04;
	opls_asp["VAL CG2 "] =    0.04;
	opls_asp["VAL C   "] =    0.04;
	opls_asp["VAL O   "] =    0.01;
	opls_asp["VAL OXT "] =   -0.01;
	opls_asp["HOH O   "] =    0.00;
	opls_asp["HOH H1  "] =    0.00;
	opls_asp["HOH H2  "] =    0.00;
	opls_asp["HOH 1HO "] =    0.00;
	opls_asp["HOH 2HO "] =    0.00;
	opls_asp["CO  C   "] =    0.04;
	opls_asp["CO  O   "] =    0.01;
	opls_asp["O2  O1  "] =    0.01;
	opls_asp["O2  O2  "] =    0.01;
	opls_asp["TERNN   "] =   -0.05;
	opls_asp["TERNHT1 "] =    0.00;
	opls_asp["TERNHT2 "] =    0.00;
	opls_asp["TERNHT3 "] =    0.00;
	opls_asp["ASP 1H  "] =    0.00;
	opls_asp["ASP 2H  "] =    0.00;
	opls_asp["ASP 3H  "] =    0.00;
	opls_asp["GLN 1H  "] =    0.00;
	opls_asp["GLN 2H  "] =    0.00;
	opls_asp["GLN 3H  "] =    0.00;
	opls_asp["LYS 1H  "] =    0.00;
	opls_asp["LYS 2H  "] =    0.00;
	opls_asp["LYS 3H  "] =    0.00;
	opls_asp["TERNCA  "] =    0.04;
	opls_asp["TERCC   "] =    0.04;
	opls_asp["TERCCA  "] =    0.04;
	opls_asp["TERCO   "] =    0.01;
	opls_asp["TERCOXT "] =   -0.01;
	opls_asp["TERCOCT1"] =   -0.01;
	opls_asp["TERCOCT2"] =   -0.01;
	opls_asp["TERCHXT "] =    0.00;

  }

#endif //ReadInput.h
