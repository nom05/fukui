! ** ./onlygrid.f90 >> Subroutine to treat with grid file only
!
!  Copyright (c) 2023  Nicolás Otero Martínez - Marcos Mandado Alonso - Ricardo A. Mosquera Castro
!  This file is part of the fukui program available in:
!      https://github.com/nom05/fukui
!
!  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
!  by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!  
!  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License along with this code.  If not, see 
!  <http://www.gnu.org/licenses/>.

subroutine onlygrid(ntime)
  use,intrinsic :: iso_fortran_env, only: output_unit,real64,int32
  use           :: grid           , only: readgrid
  use           :: main           , only: print_main,partial_time
  use           :: output         , only: shortprint_and_deallocate
  use           :: wfn            , only: read_nuclei_coord_from_wfn

  implicit none

  integer           , parameter        :: i4 = int32
  integer(kind = i4)                   :: ntime
!
! Read wfn  file (nuclei Cartesian coordinates)
!
  call print_main("Reading wfn file")
  call read_nuclei_coord_from_wfn
  call partial_time(ntime,'WFN file read in')
!
! Read grid file
!
  call print_main("Reading grid file with a set of points:")
  call readgrid
  call partial_time(ntime,'Grid file read in')
!
! Print results and deallocate
!
  call shortprint_and_deallocate

end subroutine onlygrid
