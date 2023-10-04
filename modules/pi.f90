! ** modules/pi.f90 >> Main module to process pi MOs in fukui program
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

module pi

  use,intrinsic :: iso_fortran_env,        only: output_unit,real64,int32
  use           :: omp_lib

  implicit none

  integer             ,parameter,private      :: dp = real64     &
                                               , i4 = int32      &
                                               , ou = output_unit 

  integer(kind = i4),allocatable,dimension(:) :: nompi
  integer(kind = i4)                          :: npi,nbac,nbal

  contains
    
    subroutine readpi
      use :: main,                   only: isom,print_main,lsg,sudgfchk,prepnumb,int2str,ifuk &
                                         , print_rdarr
      use :: wfn ,                   only: nom

      implicit none

      character(len  =   17)            :: text
      character(len  = 1000)            :: numb
      integer  (kind =   i4)            :: itmp,i,j,k,icde

      allocate(nompi(nom))
      if (lsg) then
         read (isom,*,iostat=i) itmp,npi
         read (isom,*,iostat=j) nompi(:npi)
         if (i.NE.0.OR.j.NE.0) stop 'ERROR: Reading .sg file'
      else
         text = 'PI   respecto al '
         npi  = 0
         k    = 0
         icde = 1
         do while (k.EQ.0.AND.icde.EQ.1)
            call sudgfchk(text,isom,icde,1  ,k)
            if (k.EQ.0.AND.icde.EQ.1) then
               npi = npi+1
               backspace(isom)
               read (isom,'(2X,I4)') nompi(npi)
            endif !! (k.EQ.0.AND.icde.EQ.1) then
         enddo !! while (k.EQ.0.AND.icde.EQ.1)
      endif !! (lsg) then
      call prepnumb(npi ,nompi  ,numb)
      call print_rdarr("PI orbitals read")
      write(ou  ,'(A)',advance='NO') trim(numb)
      write(ou  ,'(X,A,A,A)') 'PI(',int2str(npi),')'
      write(ifuk,'(A)',advance='NO') trim(numb)
      write(ifuk,'(X,A,A,A)') 'PI(',int2str(npi),')'
    end subroutine readpi

    subroutine set_pi_mo_only
      use :: wfn ,      only: nom,matom,ialfa,uhf,nbac
      use :: main,      only: lpi,sigpi,print_rdarr,ifuk,int2str,prepnumb,ifuk

      implicit none

      integer  (kind =  i4) :: i
      character(len  = 100) :: numb

      if (allocated(matom)) stop 'ERROR: matom allocated uncorrectly'
      if (lpi.AND.sigpi) then
         nbac = nom
         nbal = ialfa
         nom  = npi
         allocate(matom(nom))
         ialfa= 0
         call print_rdarr("# MO sets to ")
         write(ou  ,'(A)') int2str(nom)
         write(ifuk,'(A)') int2str(nom)
         matom(:nom) = nompi(:nom)
         if (uhf) then
!$omp parallel default(none)             &
!$omp&              if(uhf)              &
!$omp&          shared(ialfa,matom,nbal)
!$omp   workshare
            ialfa = count(matom.LE.nbal)
!$omp   end workshare
!$omp end parallel
            call prepnumb(ialfa,nompi,numb)
            call print_rdarr("Alpha PI orbitals")
            write(ou  ,'(A,X,"(",A,")")') trim(numb),int2str(ialfa)
            write(ifuk,'(A,X,"(",A,")")') trim(numb),int2str(ialfa)
         endif !! (uhf) then
      else
         allocate(matom(nom))
         nbac = nom
         matom = (/ (i, i = 1,nom) /)
      endif !! (lpi.AND.sigpi) then
    end subroutine set_pi_mo_only

end module pi
