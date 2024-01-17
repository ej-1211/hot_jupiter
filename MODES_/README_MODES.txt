

          /-----------------\
           MODES README FILE 
          \-----------------/
          
     *******************************************************************************
     *******************************************************************************
     **    December 2015                                                          **
     **                                                                           **
     **    Faculty of Mathematics and Physics                                     **
     **    University of Ljubljana, Slovenia                                      **
     **    http://www.fmf.uni-lj.si/en/                                           **
     **                                                                           ** 
     **    MODES group                                                            **
     **    http://meteo.fmf.uni-lj.si/MODES                                       **
     **    Funded by the European Research Council, Grant Agreement no. 28015     **
     **                                                                           **
     *******************************************************************************
     *******************************************************************************


==========================================================================================
MODES is a software for the analysis of global dynamical fields in (re)analyses, global 
weather forecasts and climate models. 
It is a new diagnostic tool, developed within the MODES project,
which allows one to diagnose properties of the balanced and inertio-gravity circulation
across many scales. In particular, the inertio-gravity spectrum, which has only recently 
become observable, can be studied simultaneously in the mass and wind fields while 
considering the whole atmosphere depth as represented by the models.  

The MODES software is described in the open-access paper published in Geosci. Model. Dev.:
Nedjeljka Zagar et al., 2015: Normal-mode function representation of global 3-D data sets: 
open-access software for the atmospheric research community.  
http://www.geosci-model-dev.net/8/1169/2015/gmd-8-1169-2015.pdf

MODES software and technical support for its application are available from 
http://meteo.fmf.uni-lj.si/MODES/

The software consists of a new code written mainly in Fortran 90 which needs to be 
combined with several basic libraries available in the public domain. The external
libraries needed for the software implementation are the libraries for handling the 
input data in GRIB and NetCDF formats and the LAPACK and ALFPACK libraries 
for solving the eigenvalue problem. While a version of the LAPACK (version 3.4) and 
a somewhat modified ALFPACK source code are provided with the package, the NetCDF
libraries need to be installed by the user. The ALFPACK package, which is used for
the computation of the  associated Legendre functions of the first kind, originates 
from NCAR.

NOTE: in this release only NetCDF input is supported by default, the software version 
      for reading the GRIB format is available on request. The reason is that a majority of 
      users seem to have troubles with the installation of the grib api software. 


==========================================================================================

PRE-INSTALLATION
 After unpacking the modes tar file lapack and alfpack will need to be installed. 
 Their tar files are provided including lapack version 3.4. 
 Place them where suitable and in folder lapack-3.4.0 type make and lapack and blas 
 libraries will be installed generating liblapack.a, librefblas.a, libtmglib.a files. 
 Unpack alfpack and repeat command. The default compiler is gfortran.
 For more details check original README files contained in according folders.
 NetCDF library's are not suply'd by this software so user needs to install them and add
 correct paths upon the installation.

NOTE: since software expects both 'libnetcdff.a' ('fortran' support library) and 
      'libnetcdf.a' ('C' support library) make sure that both library's are in same
      './lib' directory or modify ./NMF_MODES/Mkinclude file accordingly. 

==========================================================================================
MODES INSTALLATION 
- MODES can be installed with make. Configuration have to be set in the 
  ./NMF_MODES/Mkinclude file. A Mkinclude.gfortran examples for gfortran compiler is 
  provided. Ifort compiler is also possible but user is encouraged to use gfortran compiler
  as it is well tested and in use by MODES group. There is also Mkinclude.gfortran.grib
  which is used for compilation if the user wants to analyze input data in GRIB format. 
  BUT, in this case a few changes are needed in the code to activate the use of grib_api
  and the user should contact us. 
- It is expected from user to understand the meaning and purpose of Mkincloude file, so
  in short all paths to library's in Mkinclude file need to fit owners computer.
  However, ti is not imperative to understand Mkinclude file and "how to" support can be
  provided on request.
- After setting correct paths in same folder just type make and installation will start.
  One can check if everything went well. No errors should be given on screen and five 
  executable files should be generated in the ./NMF_MODES/bin folder with the following 
  names: 3Dexpansion, 3Dinversion, gauss, hort_struct, vsf. If not, check all paths
  inside Mkinclude file again and make sure all pre-installation requirements are met. Two
  of the most common errors are absence of the required libraries or wrong path to them.


==========================================================================================
PREPARATION OF INPUT DATA 
This section will explain how to prepare input fields for use in the MODES software.
In order to have a successful projection of input data to spectral space as well as an 
inversion user will need to have:
 - stability profile (BINARY file) oriented bottom to top
 - sigma levels (BINARY file) oriented bottom to top
 - model levels (TXT file) oriented top to bottom 
   (example: http://www.ecmwf.int/en/forecasts/documentation-and-support/91-model-levels)
The input files from any model need to contain the following fields: velocity components 
u and v, temperature T, specific humidity q, surface pressure (sp or lnsp), and surface
geopotential z. These fields can be all stored in one file per time step or one field 
per file per time step, on the user-specified Gaussian grid and with a specific file name
syntax. MODES software expects that each file ending has the form: yyyymmddhhmmss0

NOTE : If user have surface geopotential as an constant files (not date depended) file 
       can be given independently as a single file without yyyymmddhhmmss0 extension

########################################################################################
WARNING: MODES can recognize orientation of input files in NetCDF and rearrange them to 
         fit software requirements. However, if levels are indexed with numbers 1 -> N 
         MODES presumes 1 == top and N == surface. To verify if MODES loaded your data
         correctly pleas check temperature profile displayed while computing expansion.run
         MODES will print from top to bottom, consequently last line needs to be ~ 285°
########################################################################################

==========================================================================================
PREPARATIONS OF NAME LISTS
To be able to make 3D-NMF expansion one needs to prepare the input files and compute 
Gaussian grid, Vertical Structure Functions (VSF) and Hough functions. All computations
should be done from the NMF run directory. Each of sub directories contains the name-list 
(filename.cnf file) for one of the executable files in NMF_MODES/bin directory. One 
needs to adjust the name list for its own requirements, make the necessary links and 
execute the executable file. Detailed description of ”how to” is given further on.

RECOMMENDED: do not change the original ”.cnf” file but we advice to use suitable 
             extension to separate various experiments and use ”ln -sf” 
             (as shown below) to link with original ”*.cnf” file
IMPORTANT: in following text, all parameters in the namelists are described but the user 
           modifies only a number of them as described. 


==========================================================================================
Step-by-step execution of various MODES programs in order to project the 3D data based on
test input files provided with MODES so user can have solid how to reference. After test
is complete successfully user can create new .cnf files to fit his input data.
==========================================================================================
Contact and support: Damjan Jelic,  damjan.jelic@fmf.uni-lj.si

GAUSSIAN GRID

Location : NMF_RUN/Grid
File     : gauss.cnf

Source
*************************************************************************
1...  &Gaussian
2...  N           = 48
3...  gauss fname = '../../NMF_DATABASE/gauss/gauss48.data'
4...  /
*************************************************************************

Description:
2.  defines the number of meridional grid points
3.  defines the path of BINARY output file

Execution:
TEST: - link test file with original file: ln -sf gauss.cnf_test gauss.cnf
      - execute the executable file: ./gauss.run 
      - check if new file appear at specify location

PREPARE and link the new file with original file: ln -sf gauss.cnf_new gauss.cnf

==========================================================================================

VERTICAL STRUCTURE FUNCTIONS

Location : NMF_RUN/Vert
File     : vsfcalc.cnf

Source
*************************************************************************
1...  &vsfcalc cnf
2...  stab fname       = '../../NMF_DATABASE/vsf/erainterim_sigma_stability.data'
3...  vgrid fname      = '../../NMF_DATABASE/vsf/erainterim_sigma_levels.data'
4...  vsf fname        = '../../NMF_DATABASE/vsf/vsf.data'
5...  equiheight fname = '../../NMF_DATABASE/vsf/equivalent_height.data'
6...  num vmode        = 40
7...  zgrid            = 'Sigma'
8...  mp               = 60
9...  top              = 1.0
10... hstd             = 8000.0d0
11... given stability  = .true.
12... given vsf        = .false.
13... ocheck           = .true.
14... /
*************************************************************************

Description:
2. & 3. are the input BINARY files for computation of VSF. One needs to calculate both 
    stability and sigma levels.
4. & 5. are the output BINARY files of ”VSF” and ”equivalent height” respectively
6.  is number of vertical modes one wishes to use (eq. or less then model levels)
7.  type of vertical grid for vertical expansion. Sigma is currently the only option
8.  number of model levels
9.  model top in hPa (don't change)
10. scaling constant for equivalent heights - please leave it at 8000m  (don't change)
11. switch if you want to use precomputed stability file .true. (default) or .false.
12. switch if you want to use precomputed vertical structure functions .true. or
    .false. (default)
13. is a logical variable. .true. or .false. should be chosen. If .true. is chosen, 
    the information about the orthogonality will be displayed.

Execution:
TEST: - link test file with original file: ln -sf vsfcalc.cnf_test vsfcalc.cnf
      - execute the executable file: ./vsf.run 
      - check if new files appear at specify location

PREPARE and link the new file with original file: ln -sf vsfcalc.cnf_new vsfcalc.cnf

IMPORTANT : if using high vertical resolution: check the equivalent heights. Very small 
            values  (e.g. smaller than 1 meter) will be severily trapped to the equator 
            and will hardly contribute to the quality of the projection. Depending on 
            purpose, we suggest not to use equivalent depths under 5 meters. 
IMPORTANT : The provided binary file with sigma levels should be given from bottom to
            top, i.e. from σ = 1 towards σ = 0.
            
In this version, also netCDF format file with vertical structure functions is generated 
automatically.

==========================================================================================
HOUGH STRUCTURE FUNCTIONS

Location : NMF_RUN/Hort
File     : houghcalc.cnf

Source
*************************************************************************
1...  &houghcalc cnf
2...  szw              =  0
3...  ezw              = 40
4...  maxl             = 30
5...  my               = 48
6...  freq fname       = '../../NMF_DATABASE/hough/freq.data'
7...  ks mode          = 'K'
8...  ocheck           = .false.
9...  &meridional grid
10... ygrid fname      = '../../NMF_DATABASE/gauss/gauss48.data'
11... &vsf cnf
12... equiheight fname = '../../NMF_DATABASE/vsf/equivalent_height.data'
13... num vmode        = 20
14... &output
15... output gmt       = .false.
16... ofname gmt       = '../../NMF_DATABASE/hough/hough_gtm'
17... ofname bin       = '../../NMF_DATABASE/hough/hough'
18... bin combine      = 'zonal'
19... /
*************************************************************************

Description:
2. & 3. first and last zonal wavenumber included in calculation
4.  is a total number of each type of hough modes used (EIG, WIG, ROT)
5.  Number of meridional grid points
6.  Laplace tidal frequency BINARY output
7.  is a selection for obtaining the eigenvector of rotational mode in zonal 
    wavenumber 0. 'K' or 'S' should be chosen. The default value is 'K'.
8.  switch for orthogonality check of hough functions (.true. - with check 
    or .false. - without check).
10. BINARY input of meridional grid
12. BINARY input of equivalent height
13. is number of vertical modes VSF-s one wishes to use
15. whether to save output Hough functions as a text file
16. If 'output gmt ' is .true. , set the filename for it
17. Filename for Hough functions BINARY output
18. 'zonal' means that Hough function output is sorted by zonal wavenumber

Execution:
TEST: - link test file with original file: ln -sf houghcalc.cnf_test houghcalc.cnf
      - execute the executable file: ./hort_struct.run 
      - check if new files appear at specify location

PREPARE and link the new file with original file: ln -sf houghcalc.cnf_new houghcalc.cnf

HINT : the computation of Hough functions can take a while depending on the power of 
       the computer and 1.-5. parameters so be prepared to go for coffee or lunch ;)

In this version, also netCDF format files with horizontal structure functions should be 
generated automatically.

==========================================================================================
3D NMF EXPANSION

Location : NMF_RUN/Proj
File     : normal.cnf

Source
*************************************************************************
1...  &normal cnf
2...  nx              = 96
3...  ny              = 48
4...  nz              = 60
5...  nstep           =  1
6...  coef3DNMF fname = '../../NMF_DATABASE/coef/Hough_coeff_'
7...  output 3DNMF    = .true.
8...  saveps          = .false.
9...  savemeant       = .false.
10... ps fname        = '../../NMF_DATABASE/coef/Ps_'
11... meant fname     = '../../NMF_DATABASE/coef/Tmean_'
12... saveasci        = .false.,
13... afname          = '../../NMF_DATABASE/coef/InputData_'
14... aformat         = '(96E20.4,1x)',
      savehoughwind   = .false.,
15... savenc          = .true.
16... ncname          = '../../NMF_DATABASE/coef/Hough_coeff_'
17... &time
18... datetype = 'yyyymmddhh'
19... syear    =  2000
20... smon     =    01
21... sday     =    01
22... shour    =    00
23... smins    =    00
24... ssec     =    00
25... slen     =    00
26... eyear    =  2000
27... emon     =    01
28... eday     =    01
29... ehour    =    00
30... emins    =    00
31... esec     =    00
32... elen     =    00
33... dt       = 21600
34... &input data
35... dataformat input  = 'netcdf'
36... orig              = 'ECMWF'
37... zgrid type        = 'hybrid'
38a.. ifile_orog        = ''
38b.. numoffile         = 1,
39... ifile head(1)     = '../../NMF_DATABASE/input_data/erainterim_moda_L60_N24_'
40... ifile head(2)     = ''
41... &hough cnf
42... hough fname       = '../../NMF_DATABASE/hough/hough'
43... num zw = 40,
44... maxl   =  30,
45... my     = 48,
46... &meridional grid
47... ygrid fname       = '../../NMF_DATABASE/gauss/gauss48.data'
48... &vsf cnf
49... vsf fname         = '../../NMF_DATABASE/vsf/vsf.data'
50... vgrid fname       = '../../NMF_DATABASE/vsf/erainterim_sigma_levels.data'
51... equiheight fname  = '../../NMF_DATABASE/vsf/equivalent_height.data'
52... model level fname = '../../NMF_DATABASE/vsf/erainterim_model_levels.dat'
53... num vmode         = 20
54... zgrid             = 'Sigma'
55... mp                = 60
56... top               = 1.0
57... given vsf         = .true.
58... /
*************************************************************************

Description:
2.  number of zonal grid points
3.  number of meridional grid points
4.  number of vertical grid points
5.  number of time steps you want to compute (use 1).
6.  3D-NMF expansion BINARY output file (spectral space). Date is added automatically.
7.  whether to save the output file (true) or not (falce)
8.– 11. whether to save the surface pressure and mean vertical temperature files or not.
    If ”true” set the filenames
12. & 13. allows one to save input files (u,v,P,t,q) in ASCII format. By default 
    u, v & P are saved. To change one needs to go into the code.
14. is format of ASCII file from 13. Structure is: (ny*nz*3,nx) in double precision
    (3 = u + v + P)
15. & 16. switch and path for saving NetCDF output.
17. in these part date will be generated and placed as an extension on all output and
    input files within name list. Thus one needs to prepare input files in expected
    sequence yyyymmddhhmmss0.
18. – 32. is date definition. ”s” mean starting date(year, month, day, hour, minute,
    second) ”e” means end. ”len” is forecast length.
33. is time step in seconds. For 6 h time step it is ”21600”, for month it is ”2678400”
35. format of data files for expansion is ”netcdf”, grib can be on request
36. source of data, can be 'ECMWF', 'MERRA', 'CCSM', 'CMIP5'. See apendix OROGIN
37. defines the vertical grid type: ”hybrid” or ”sigma”
38a. if surface geopotential ia prowided as contant field just add a path and file name
38b. how many input files one has. If all variables are in one file use value ”1”.
39. & 40. definition of location of input files. The software will use all files 
    connected to dates defined in 17.
42. location of Hough functions
43. is a total number of zonal wavenumbers (including zonal wavenumber 0)
44. is a number of each hough modes (EIG, WIG, ROT)
45. number of meridional grid points
47. BINARY input of meridional grid
49. BINARY input of VSF
50. BINARY input of sigma levels
51. BINARY input of equivalent height
52. TXT file containing number of model half levels, ”a” and ”b” coefficients for hybrid
    coordinate, half level pressures and full level pressures
53. number of vertical modes one wishes to use
55. number of model levels

Execution:
TEST: - link test file with original file: ln -sf normal.cnf_test normal.cnf
      - execute the executable file: ./expansion.run
      - check if new files appear at specify location

PREPARE and link the new file with original file: ln -sf normal.cnf_new normal.cnf

IMPORTANT : during the visualization of data be sure to read them correctly because the
            dimensions of output data are: ”double” for expansion and ”real*4” for
            inversion
In this version, also netCDF format files with the Hough expansion coefficients can be 
generated if savenc is set .true.

==========================================================================================
3D NMF INVERSION

Location : NMF_RUN/Inv
File     : normal inverse.cnf

Source
*************************************************************************
1...  &normal cnf inverse
2...  nx              = 96
3...  ny              = 48
4...  nz              = 60
5...  nstep           =  1
6...  coef3DNMF fname = '../../NMF_DATABASE/coef/Hough_coeff_'
7...  inverse fname   = '../../NMF_DATABASE/inverse/Erai_allk_'
8...  inv2hybrid      = .false.
9...  ps fname        = '../../NMF_DATABASE/coef/Ps_'
10... meant fname     = '../../NMF_DATABASE/coef/Tmean_'
11... saveasci        = .true.,
12... afname          = '../../NMF_DATABASE/inverse/Erai_ascii_allk_'
13... aformat         = '(96E16.4,1x)'
      savenc          = .true.
      ncname          = '../../NMF_DATABASE/inverse/Erai_allk_'
14... &time
15... datetype = 'yyyymmddhh'
16... syear    =  2000
17... smon     =    01
18... sday     =    01
19... shour    =    00
20... smins    =    00
21... ssec     =    00
22... slen     =    00
23... eyear    =  2000
24... emon     =    01
25... eday     =    01
26... ehour    =    00
27... emins    =    00
28... esec     =    00
29... elen     =    00
30... dt       = 21600
31... &hough cnf
32... hough fname    = '../../NMF_DATABASE/hough/hough'
33... num zw   =   40
34... maxl     =   30
35... my       =   48
36... shough   =    1
37... ehough   =   90
38... &meridional grid
39... ygrid fname   = '../../NMF DATABASE/gauss/gauss128.data'
40... &filter cnf
41... eig_n_s  =   610
42... eig_n_e  =   100
43... wig_n_s  =   610
44... wig_n_e  =   100
45... rot_n_s  =   610
46... rot_n_e  =   100
47... kmode_s  =   220
48... kmode_e  =    85
49... vmode_s  =   410
50... vmode_e  =    70
51... &vsf cnf
52... vsf fname         = '../../NMF_DATABASE/vsf/vsf.data'
53... vgrid fname       = '../../NMF_DATABASE/vsf/erainterim_sigma_levels.data'
54... equiheight fname  = '../../NMF_DATABASE/vsf/equivalent_height.data'
55... model level fname = '../../NMF_DATABASE/vsf/erainterim_model_levels.dat'
56... num vmode         = 20
57... zgrid             = 'Sigma'
58... mp                = 60
59... top               = 1.0
60... given vsf         = .true.
61... /
*************************************************************************

Description:
2.  number of zonal grid points
3.  number of meridional grid points
4.  number of vertical grid points
5.  number of time steps you want to compute (use 1).
6.  input BINARY file of 3D-NMF expansion (spectral space)
7.  output BINARY file of 3D-NMF inversion (physical space)
8.  to interpolate output filesfrom sigma levels to hybrid levels (unavaulable)
9. & 10.  to save mean surface pressure and mean temperature
11. – 13. to save ASCII file of inverse (u,v,P) fields in format: (ny*nz*3,nx) 
    single precision (3 = u + v + P).
15.– 30. is using inputs from normal.cnf!
32. location of Hough functions
33. is index of maximum zonal wavenumber, which includes zonal wavenumber 0
34. is a number of each hough modes (EIG, WIG, ROT)
35. number of meridional grid points
36. & 37. begining and end of all hough modes.
39. BINARY input of meridional grid
40.– 50. filtering of 3D expansion. For more informations see Appendix FILTERING
51.– 60. is using inputs from normal.cnf!

Execution:
TEST: - link test file with original file: ln -sf normal_inverse.cnf_test normal_inverse.cnf
      - link test file with original file: ln -sf ../Proj/normal.cnf_test normal.cnf
      - execute the executable file: ./inversion.run
      - check if new files appear at specify location

PREPARE and link the new file with original file: ln -sf normal_inverse.cnf_new normal_inverse.cnf
PREPARE and link the new file with original file: ln -sf ../Proj/normal.cnf_new normal.cnf

HINT : 3D-NMF inversion is using time from normal.cnf name list so be sure to adjust 
       the date and time correctly

IMPORTANT : during the visualization of data be sure to read them correctly because the
            dimensions of output data are: ”double” for expansion and ”real*4” for
            inversion

In this version, also netCDF format files with the results of filtering be 
generated if savenc is set .true.


==========================================================================================
GENERAL NOTES
Don't hesitate to ask for help, modes team is always at your disposal.
Gave us feedback and help us improve software to fit new user even better.

==========================================================================================
----------------------------------APENDIX-------------------------------------------------
==========================================================================================
ORIGIN
Each data set comes with different variable names. In this version we covered some frequent 
data sets like ECMWF's data, MERRA, CCSM and CMIP5 project data. If MODES complains about 
finding some variable its most likely that issue is in 'origin'. If your variable names 
don't mach with proposed origin, pleas contact us and we will male new origin for your needs.

ECMWF : "lon", "lat", "levels", "time", "var129", "var130", "var131", "var132", "var133", "var152" 
MERRA : "lon", "lat", "levels", "time", "u", "v", "t", "qv", "ps", "PHIS"
CCSM  : "lon", "lat", "levels", "time", "U", "V", "T", "Q", "PS", "PHIS"
CMIP5 : "lon", "lat", "levels", "time", "ua", "va", "ta", "hus", "ps", "orog"

==========================================================================================
FILTERING
This MODES version contains filter which filters out data one dose not need. To obtain balanced 
circulation in our example on needs to pull out all unbalanced circulation using:

   eig_n_s =   1 -> first mode to be filtered
   eig_n_e =  30 ->  last mode to be filtered
   wig_n_s =   1
   wig_n_e =  30
   rot_n_s = 610
   rot_n_e = 100
   kmode_s = 220
   kmode_e =  85
   vmode_s = 410
   vmode_e =  70

As long as starting mode number is greater then end mode number no filtering will occur for
that group of modes.

  













