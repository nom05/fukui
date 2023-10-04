! ** ./goon.f90 >> Main subroutine of the fukui program
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

subroutine goon(ntime)

  use,intrinsic :: iso_fortran_env,   only: output_unit,real64,int32
  use           :: main           ,   only: print_main,sigpi,nproc,partial_time,itype,txttype,ipartdef,iato
  use           :: wfn            ,   only: readwfn,nato
  use           :: grid           ,   only: readgrid
  use           :: pi             ,   only: readpi,set_pi_mo_only
  use           :: computation    ,   only: orde_ic,skip_gaussians &
                                          , Dens                   & !! Density
                                          , DensAom                & !! Density + AOMs
                                          , DensDipMom             & !! Density + Dipole Moment
                                          , DensDipMomAom            !! Density + AOMs + Dipole Moment
  use           :: output         ,   only: print_and_deallocate

  implicit  none

  integer             ,parameter         :: i4        = int32

  integer  (kind = i4)                   :: ntime,ipart
!
! Read files
!
!  * WFN
!
  call print_main("Reading wave function file:")
  call readwfn
  call partial_time(ntime,'Wave function read in')
  if (iato.GT.nato) stop 'ERROR: atom label > # atoms!!!'

!  * Grid file
!
  call print_main("Reading grid file with a set of points:")
  call readgrid
  call partial_time(ntime,'Grid file read in')

!  * Sigma/pi separation
!
  if (sigpi) then
     call print_main("Opening sigma/pi separation file:")
     call readpi
  endif !! (sigpi) then
!
! Previous preparation
!
!   * Obtain the first and last gaussians centered on each atom
!
  call orde_ic
!
!   * Skipping gaussian functions
!
  call print_main("Preparing calculation:")
  call skip_gaussians
!
!   * Only pi if specified
!
  call set_pi_mo_only
!
! Calculation starts
!
!     allocate(poom(nom))                                               !! R1
!     call cpu_time(ptime)
  if (nproc.EQ.0)      then
     ipart = ipartdef
  else if (nproc.LT.8) then
     ipart = 2*nproc
  else
     ipart =   nproc
  endif !! (nproc.EQ.0) then

  call partial_time(ntime,'Intermediate steps in')
  call print_main("Calculation type: "//trim(txttype(itype)))
  call print_main("Computing grid points:")

  select case(itype)
  !! Types of calculations:
  !!   xx
  !!   ||-> AOM
  !!   |--> Dip
  !!
  !! bin-dec
  !!   *  00-  0: e density
    case(0)
      call             Dens(ipart)

  !!   *  01-  1: e density + AOM
    case(1)
      call          DensAom(ipart)

  !!   *  10-  2: e density + dip moment
    case(2)
      call       DensDipMom(ipart)

  !!   *  11-  3: e density + dip moment + AOM
    case(3)
      call    DensDipMomAom(ipart)

    case default
      stop 'ERROR: The requested calculation is not implemented'
  end select !! case(itype)

  call partial_time(ntime,'Grid points computed in')

  call print_and_deallocate

end subroutine goon
