This repository contains Phantom of Ramses (PoR) code with the Milky Way galaxy model in Modified Newtonian dynamics (MOND), implemented as a background analytic density.

The Milky Way model is adapted from Banik et al (2023): https://arxiv.org/pdf/2204.09687. The template is based on the initial conditions generated from mondified version of the DICE code, available here: https://bitbucket.org/SrikanthTN/bonnpor/src/master/

This patch can be used as follows:

1) Compile the Makefile with
   _PATCH = ../patch/phantom_MW_staticparts_
   
if the compilation fails due to pointer and array mismatch errors in clump_finder.o and other files, do not modify the files, just use an additional flag
_-fallow-argument-mismatch_ for the compiler.

3) A new namelist file is added in PoR_main/POR_namelist/MW_analytic.nml.
   This namelist file contains a new params block:

  **&MOND_DISK_PARAMS**
  
  _mass_disk=9.15d10_ ! --- Total mass of the galaxy --- ! 
  
  _f_inner=0.8236_    ! -- Fraction of mass in the inner disk ---!
  
  _f_outer=0.1764_    ! --- Fraction of mass in the outer disk --!
  
  _rd_inner= 1.29_    ! -- Scale length of the inner disk --- !
  
  _rd_outer = 4.20_   ! -- Scale length of the outer disk --- !
  
  _rmax = 300_          ! -- maximum radius upto which density is added --!!
  
  /

  Make sure that _gravity_type_=-1 and _mond=.true_. flags are on, else this routine will not be called.

  The parameters in this block are read and a series of subroutines are called in mond_disk_density.f90 file. 
  All the subroutines first initialise a MOND disk and then compute the density profile upto rmax and store it as a table. 
  The analytic density routine then reads the table and uses cubic interpolation function to interpolate densities to desired radii. 
  This is a less expensive way to compute densitites than to call all the subroutines to compute densities for each cell, which becomes too expensive for bigger boxes.

  Soon after launch, the code generates two output files 'mond_density_table.dat' and 'mond_rotation_table.dat' files are generated. 
  The former contains radius and density values while the latter contains quantities like radius, gN/a0, Cumulative mass, MOND Vc, Newtonian Vc etc,
  which can be used to verify the analytic and RAMSES quantities. 

  This analytic density template is based on my previous work on Coma cluster analytic density in POR, available here: https://github.com/SrikanthNagesh/Coma_analytic_density_PoR

  CITATION POLICY: Users can cite Nagesh et al (2024) (https://ui.adsabs.harvard.edu/abs/2024A%26A...690A.149N/abstract), for now.

  CONTACT: srikanth.nagesh@epfl.ch; tnsrikanth1998@gmail.com
  

  
