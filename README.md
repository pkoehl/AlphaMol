# AlphaMolWrap

# **1. Background**
Cancel changes
This repository is a fork of {https://github.com/pkoehl/AlphaMol}[pkoehl/AlphaMol] which is

# **2. The Program: License**

I distribute the source code for the whole program, including all
subroutines, under the LGPL licensing code. Please read the text of
the license (LGPL-LICENSE.txt provided with the distribution). 
Basically, you are allowed to use and modify this code freely. 
You are also allowed to distribute it, provided that you provide 
the source code for all functions involded in Alphamol.

# **3. The Program: Content**

In the archive, you will find:

1. **Makefile**: the makefile for compiling AlphaMol

2. **bin**: the directory that will contain the executable AlphaMol
            (currently empty)

3. **project**: the directory that contains the main for the program; 
it includes a subdirectory "include" and a directory "src"

4. **InOut**: directory for reading procedures

5. **DataStruct**: directory for the different data structures
  (mostly for Delaunay triangulation and Alpha Complex)

6. **Delcx**: directory for computing the regular (i.e. weighted)
   Delaunay triangulation of a set of weighted balls in R^3

7. **Alphacx**: directory for computing the dual complex
(alpha complex with alpha=0) once the Delaunay is known

8. **Volume**: directory for computing the geometric
measures of the union of balls, once the dual complex is known

# **4. The Program: Installation**

The Makefile provided with the distribution is intentionally kept
simple. 

To install, just type:

*make clean*
*make*

Note that this will use the standard g++ compiler
that should be installed by default on your Linux / Unix / MacOS X

# **5. Running the Program**

The name of the executable is AlphaMol; it is located in the bin
subdirectory. Here I define:
- the inputs needed by the program
- typical commands to run the program

*5.a. Input file*

AlphaMol can read in three types of data:

i) A simple file containing information about the union of
balls under study, with one line per ball, and on each
line the x, y, z coordinates of the center of the ball
and its radius (all four values separated by space)

ii) A PDB file (needs to have the extension .pdb)
Radii are then assigned based on the OPLS force field

iii) A PQR file(needs to have the extension .pqr)
Radii are read directly from the file

*5.b. Typical runs*

In the directory "examples" just type:

**../bin/AlphaMol**

to see the list of options. A typical run on a CRD file

**../bin/AlphaMol -i Vol100.crd -r 0 -d 1**

with:
	-r 0: indicates that we do not inflate the balls
        -d 1: indicates that we compute the derivatives

A typical run on a PDB file:

**../bin/AlphaMol -i 1tim.pdb -r 1.4 -d 1**

with:
	-r 1.4: indicates that we inflate the balls by 1.4
		(a typical procedure for computing accessible
		surface areas, where 1.4 is the probe radius
		for a water molecule)
        -d 1  : indicates that we compute the derivatives

A typical run on a PQR file:

**../bin/AlphaMol -i mol1.pqr -r 1.4 -d 1**

with:
	-r 1.4: indicates that we inflate the balls by 1.4
		(a typical procedure for computing accessible

In the directory examples, there are several .crd, one .pdb
and one .pqr file. The files .info shows the output of the
program AlphaMol when run on those files, with r=0 for the
crd files, and r=1.4 for the pdb and pqr files.

*5.c Output files*

This is a misnomer! At this stage, AlphaMol does not output any files
and all results are provided on the screen. AlphaMol is more
designed to be used as a library in other programs than as 
a standalone program. As it is provided under LGPL format,
you can modify / use it as you want!

**4. Disclaimer**

The program is provided "as is". Contact koehl@cs.ucdavis.edu
if you have any problems.
