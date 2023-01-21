! ** modules/output.f90 >> Main module to print the output file in fukui program
!
!  Copyright (c) 2022  Nicolás Otero Martínez - Marcos Mandado Alonso - Ricardo A. Mosquera Castro
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

module output

  use :: main                                 ,only: i4,dp,ou,ifuk,lskpr,igrid,lskp,sigpi,isom,lpi,laom
  use :: computation                          ,only: popul,poom,s,int2str
  use :: grid                                 ,only: rhpt,rhptsk
  use :: pi                                   ,only: nompi,npi
  use :: wfn                                  ,only: nom,matom,pel,nbac

  implicit none

  integer  (kind = i4),dimension(2),parameter     :: outmat=(/ ou, ifuk /)
  integer  (kind = i4)                            :: iom,jom
  real     (kind = dp)                            :: tmp
  
  contains

    subroutine print_and_deallocate

      implicit none

      integer(kind = i4)                          :: i
      character( len = 50)                        :: label

      do i = 1,2
         write (outmat(i),'(/,X,"RESULTS OF THE INTEGRATION")' )
         write (outmat(i),'(4X,"POP = ",7X,1PD21.14," au")'    ) popul
         write (outmat(i),'(4X,"REF POP = ",3X,1PD21.14," au")') rhpt
         write (outmat(i),'(20X,"FUKUI = ",F6.3," au")'        ) popul-rhpt
         if (lskpr.OR.(igrid.EQ.0.AND.lskp)) then
            write (outmat(i),'(/,4X,"Considering SKIPPING:")'         )
            write (outmat(i),'(4X,"REF POP = ",3X,1PD21.14," au")'    ) rhptsk
            write (outmat(i),'(4X,"DIF REF POP = ",   1PD9.3  ," au")') rhpt-rhptsk
            write (outmat(i),'(20X,"FUKUI = ",F6.3," au")'            ) popul-rhptsk
         endif !! (lskp) then
      enddo !! i = 1,2
      write (ifuk,'(/,4X,"ORBITAL CONTRIBUTIONS:")')
      label(:) = ''
      label    = '(4X,"N(",A'//int2str(ceiling(log10(dble(nom))))//',") = ",F12.6," au")'
      do i=1,nom
         write (ifuk,label) int2str(matom(i)),poom(i)*pel(i)
      enddo !! j=1,nom
      if (sigpi.AND..NOT.lpi) then
!$omp parallel default(none) shared(tmp,poom,nompi,npi,pel)
!$omp   workshare
         tmp = sum(poom(nompi(:npi))*pel(nompi(:npi)))
!$omp   end workshare
!$omp end parallel
         do i = 1,2
            write (outmat(i),'(/,10X,"SIGMA/PI CONTRIBUTIONS:",2X,I4," molecular orbitals with PI symmetry", &
                              &/,10X,"SIGMA",F10.5,5X,"PI",F10.5)') npi,popul-tmp,tmp
         enddo !! i = 1,2
         if (allocated(nompi)) deallocate(nompi)
      else if (sigpi) then
         write (ifuk,'(/,X,">> SIGMA/PI CONTRIBUTIONS skipped, only working with pi density as requested",/)') 
         if (allocated(nompi)) deallocate(nompi)
      endif !! (sigpi) then
      if (allocated(matom)) deallocate(matom)
      if (allocated(poom))  deallocate(poom)
      if (allocated(pel))   deallocate(pel)

      if (laom) then
         write (ifuk,'(/,10X,"The Atomic Overlap Matrix",/,/,/)')
         do iom=1,nbac
            write (ifuk,'(8E20.12)') (s(jom+iom*(iom-1)/2),jom=1,iom)
         enddo !! iom=1,nom
         deallocate(s)
      endif !! (laom) then
    end subroutine print_and_deallocate
  
end module output
