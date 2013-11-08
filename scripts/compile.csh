#!/bin/csh -f

#--------------------------------------------------------------------------
# Script compile.csh
#
# Purpose: Top level script for compiling driver and subroutines within 
#          current working directory.
#
# Method: 1) Link driver and subroutines to current directory
#         2) Compile routines in current directory
#
# History: 08/01/08  Original version. Henry R. Winterbottom
#          05/29/09  Updated to include dynamic linking and updating of 
#                    internal routines for external library paths. 
#                    Henry R. Winterbottom
#
#-------------------------------------------------------------------------

# Define shell environment

unlimit

# Define current working directory

set pwd = `pwd`

# Obtain configure file name for respective executable

set configurefile = `ls ${pwd}/configure/configure.*`

# Link configure file to current working directory

ln -sf ${configurefile} .

# Link makefile for executable to current working directory

ln -sf ${pwd}/scripts/Makefile ${pwd}/Makefile

# Define configure file for respective executable

set configurefile = `ls ${pwd} | grep "configure\."`

# Define netcdf path using environment variable

set netcdf = `cat ${configurefile} | grep "NETCDFPATH" | awk '{print $3}'`

# Define FFTW path using environment variable

set fftw = `cat ${configurefile} | grep "FFTWPATH" | awk '{print $3}'`

# Link HYCOM Fortran 90 routines to current working directory

ln -sf ${pwd}/src/*.F .

# Link all Fortran 90 routines to current working directory

ln -sf ${pwd}/src/*.f90 .

# Link all Fortran 77 routines to current working directory

ln -sf ${pwd}/src/*.f .

# Link all include files to current working directory

ln -sf ${pwd}/src/*.inc .

# Compile routines

make all

# Remove all links to fortran routines in current working directory

rm ${pwd}/*.F   >& /dev/null
rm ${pwd}/*.f   >& /dev/null
rm ${pwd}/*.f90 >& /dev/null
rm ${pwd}/*.h   >& /dev/null
rm ${pwd}/*.inc >& /dev/null

# Remove all object files created during compilation

rm ${pwd}/*.o >& /dev/null

# Remove all module files created during compilation

rm ${pwd}/*.mod >& /dev/null

# Remove links to netcdf files if required

rm ${pwd}/include >& /dev/null
rm ${pwd}/lib     >& /dev/null

# Remove links to NEMSIO files if required

rm ${pwd}/nemsio_include >& /dev/null
rm ${pwd}/nemsio_lib     >& /dev/null

# Remove links to w3 files if required

rm ${pwd}/w3lib >& /dev/null

# Remove configure file for executable

rm ${configurefile} >& /dev/null

# Remove makefile for executable

rm ${pwd}/Makefile >& /dev/null
