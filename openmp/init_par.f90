! ** openmp/init_par.f90 >> Subroutine to check OpenMP parallelization
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

subroutine init_par

  use,intrinsic :: iso_fortran_env,only: output_unit,real64,int32
  use omp_lib !! Enabling OpenMP functions !!
  use           :: main           ,only: nproc
  implicit none
  integer             ,parameter      :: dp        = real64
  integer             ,parameter      :: i4        = int32  &
                                       , ou        = output_unit
  integer(kind = i4)                  :: nprocmax,n_threads,id
  real   (kind = dp)                  :: wtime

  wtime = omp_get_wtime()
  write(ou,'(/,X, &
    & "**************** THREAD INFORMATION ******************" &
    &                                                            )')
  write(ou,'(5X,">> Job running using OpenMP.")')
  nprocmax = omp_get_num_procs()
  write(ou,'(5X,">> The number of processors is",X,5("."),I4)') nprocmax
  if (nprocmax.LT.nproc) then
     write(ou,'(/,15X, &
    &     "** WARNING **",I4,"=spec.proc >",I4,"=av.proc.")') nproc,nprocmax
     write(ou,'(25X,"av.proc. will be set",/)')
     nproc = nprocmax
  else if (nproc.EQ.0) then
     nproc = nprocmax
  endif !! (nprocmax.LT.nproc) then
  call    OMP_SET_NUM_THREADS(nproc)
  n_threads = omp_get_max_threads()
  write(ou,'(5X,">> The number of threads is"   ,X,8("."),I4)') n_threads
!
!  INSIDE THE PARALLEL REGION, have each thread say hello.
!
!$omp parallel private ( id )
!
!  Have each thread say hello.
!
  id = omp_get_thread_num()
  write(ou,'(25X,"Hello from process",I8)') id
!$omp end parallel
!
!  Finish up by measuring the elapsed time.
!
  wtime = omp_get_wtime()-wtime
  write(ou,'(5X,">> Elapsed wall clock time"   ,X,8("."),G11.4)' ) wtime
  write(ou,'(/,X,54("*"))')

end subroutine init_par
