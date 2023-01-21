##
##  CMakeFile for fukui.x
##  Vigo, 2022-05-10
##
##  Copyright (C) 2022 Nicolás Otero Martínez     -    nom05 (at) uvigo.es
##                     Marcos Mandado Alonso      -  mandado (at) uvigo.es
##                     Ricardo A. Mosquera Castro - mosquera (at) uvigo.es
## 
##  fukui.x is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##  
##  fukui.x is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public License
##  along with fukui.x.  If not, see <http://www.gnu.org/licenses/>.
##
##  To install this software on your system, see the INSTALL file
##  or the manual in the root directory of this package.
##
# CMake project file for FUKUI

cmake_minimum_required(VERSION 2.8.12)

## >> DEFAULT COMPILER << ## Force compiler you want
# cmake -DCMAKE_Fortran_COMPILER=gfortran src_directory

message("CMAKE_Fortran_COMPILER is ${CMAKE_Fortran_COMPILER}")
message("CMAKE_GENERATOR_FC is ${CMAKE_GENERATOR_FC}")

## >> SETTINGS         << ##
project(fukui)
enable_language (Fortran)

set(CMAKE_BUILD_TYPE Release)
message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
# Of course, change “Release” (O3) with “Debug” (no optimizations) for debug builds.

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

find_package(OpenMP)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -Wall -Wextra")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -warn all -static")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

message("CMAKE_Fortran_COMPILER_VERSION: " "${CMAKE_Fortran_COMPILER_VERSION}")

if (OPENMP_FOUND)
   message("OPENMP FOUND")
   message("OPENMP FLAGS: " ${OpenMP_C_FLAGS})
   set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
   set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
   if (CMAKE_BUILD_TYPE MATCHES "RELEASE")
       message("${CMAKE_Fortran_FLAGS_RELEASE}")
   elseif (CMAKE_BUILD_TYPE MATCHES "DEBUG")
       message("${CMAKE_Fortran_FLAGS_DEBUG}")
   endif (CMAKE_BUILD_TYPE MATCHES "RELEASE")
   set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

##################################################################################
### >> MODULE  subrout. << ##
add_subdirectory(modules)
set(mdl modules)
set (MODULES
   ${mdl}/computation.f90
   ${mdl}/grid.f90
   ${mdl}/main.f90
   ${mdl}/pi.f90
   ${mdl}/wfn.f90
   ${mdl}/output.f90
)

add_subdirectory(openmp)
set(omp openmp)
set (OMP
      ${omp}/init_par.f90
)
message("OPENMP FOUND")
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
if (CMAKE_BUILD_TYPE MATCHES "RELEASE")
       message("${CMAKE_Fortran_FLAGS_RELEASE}")
elseif (CMAKE_BUILD_TYPE MATCHES "DEBUG")
       message("${CMAKE_Fortran_FLAGS_DEBUG}")
endif (CMAKE_BUILD_TYPE MATCHES "RELEASE")
   set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

## >> MAIN COMPILATION << ##
add_executable(fukui.x 
      ${MODULES}
      fukui.f90
      help.f90
      goon.f90
      ${OMP}
)
