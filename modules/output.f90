! ** modules/output.f90 >> Main module to print the output file in fukui program
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

module output

  use :: main                                 ,only: i4,dp,ou,ifuk,lskpr,igrid,lskp,sigpi,isom &
                                                   , lpi,laom,ldip,iato,coord,dipgrid,dipgridsk
  use :: computation                          ,only: popul,poom,s,int2str,dmi
  use :: grid                                 ,only: rhpt,rhptsk
  use :: pi                                   ,only: nompi,npi
  use :: wfn                                  ,only: nom,matom,pel,nbac,uhf,ialfa,nuclei

  implicit none

  integer  (kind = i4),dimension(2),parameter     :: outmat=(/ ou, ifuk /)
  integer  (kind = i4)                            :: iom,jom
  real     (kind = dp)                            :: tmp,tmpa,tmpb,tmpaa,tmpbb
  real     (kind = dp),dimension(3)               :: ctgrid,ctgridsk
  
  contains

    subroutine print_and_deallocate

      implicit none

      integer  (kind = i4)                        :: i
      character( len = 50)                        :: label

      if (ldip) then ! Charge transfer contribution to the dip. mom.
        ctgrid   = -rhpt   * coord ! Original grid
        ctgridsk = -rhptsk * coord ! Original grid skipping points
        coord    = -popul  * coord ! Integration
      endif !! (ldip) then

      do i = 1,2
         write (outmat(i),'(/,X,"RESULTS OF THE INTEGRATION")',advance='no')
         if (iato.GT.0) write (outmat(i),'(X,"(ATOM",X,A,")")',advance='no') int2str(iato)

! >> Total(+ alpha + beta) 
         write (outmat(i),'(/,4X,"POP",15X,"=",3X,1PD21.12," au")'         ) popul
         if (uhf) then
!$omp parallel default(none) shared(tmpa,tmpb,ialfa,poom,nompi,npi,pel)
!$omp   workshare
            tmpa = sum(poom(       :ialfa) * pel(       :ialfa))
            tmpb = sum(poom(ialfa+1:     ) * pel(ialfa+1:     ))
!$omp   end workshare
!$omp end parallel
            write(outmat(i),'(&
              & 6X,"* ALPHA POP",5X,"=",3X,1PD21.12," au",/,           &
              & 6X,"* BETA  POP",5X,"=",3X,1PD21.12," au"     )'     ) &
                  tmpa,tmpb
         endif !! (uhf) then
         write (outmat(i),'(4X,"REF POP",11X,"=",3X,1PD21.12," au")') rhpt
         write (outmat(i),'(/,20X,"FUKUI = ",F6.3," au")'           ) popul-rhpt

! >> Dipole moment
         if (ldip) then
            if (iato.GT.0) then
               write (outmat(i),'(2(/),4X,"DIP MOM")')
               write (outmat(i),'(6X,"* INTRINSIC",7X ,"= ",3(1PD17.8,X)," au")'   ) dmi(:)
               write (outmat(i),'(6X,"* CHARGE TRANSFER = ",3(1PD17.8,X)," au")'   ) coord(:)
               write (outmat(i),'(6X,"* TOTAL",11X    ,"= ",3(1PD17.8,X)," au")'   ) coord(:) + dmi(:)
               write (outmat(i),'(     4X,"REF DIP MOM")')
               write (outmat(i),'(6X,"* INTRINSIC",7X ,"= ",3(1PD17.8,X)," au")'   ) dipgrid(:)
               write (outmat(i),'(6X,"* CHARGE TRANSFER = ",3(1PD17.8,X)," au")'   ) ctgrid(:)
               write (outmat(i),'(6X,"* TOTAL",11X    ,"= ",3(1PD17.8,X)," au")'   ) ctgrid(:) + dipgrid(:)
            else
               write (outmat(i),'(/,4X,    "DIP MOM = ",3(1PD17.8,X)," au")'       ) dmi(:)
               write (outmat(i),'(  4X,"REF DIP MOM = ",3(1PD17.8,X)," au")'       ) dipgrid(:)
               write (outmat(i),'(10X,"** NOTE: Atom not specified, intrinsic dip. mom. unavailable")')
            endif !! (iato.GT.0) then
         endif !! (ldip) then

! >> Comparison with skipping
         if (lskpr.OR.(igrid.EQ.0.AND.lskp)) then
            write (outmat(i),'(/,4X,"Considering SKIPPING:")'         )
            write (outmat(i),'(/,4X,"REF POP   = ",3X,1PD21.14," au")') rhptsk
            write (outmat(i),'(4X,"DIF REF POP = ",   1PD9.3  ," au")') rhpt-rhptsk

! >> Fukui index if applicable
            write (outmat(i),'(20X,"FUKUI = ",F6.3," au")'            ) popul-rhptsk

! >> Dipole moment with skipping
            if (ldip) then
               if (iato.GT.0) then
                  write (outmat(i),'(/,4X,"REF DIP MOM")')
                  write (outmat(i),'(6X,"* INTRINSIC",7X ,"= ",3(1PD17.8,X)," au")'    ) dipgridsk(:)
                  write (outmat(i),'(6X,"* DIF INTRINSIC",3X ,"= ",3(1PD17.8,X)," au")') dipgrid(:)-dipgridsk(:)
                  write (outmat(i),'(6X,"* CHARGE TRANSFER     = ",3(1PD17.8,X)," au")') ctgridsk(:)
                  write (outmat(i),'(6X,"* DIF CHARGE TRANSFER = ",3(1PD17.8,X)," au")') ctgrid(:)-ctgridsk(:)
                  write (outmat(i),'(6X,"* TOTAL",11X    ,"= ",3(1PD17.8,X)," au")'    ) dipgridsk(:)+ctgridsk(:)
                  write (outmat(i),'(6X,"* DIF TOTAL",7X ,"= ",3(1PD17.8,X)," au")'    ) dipgrid(:)+ctgrid(:) &
                                                                                        -ctgridsk(:)-dipgridsk(:)
               else
                  write (outmat(i),'(  4X,"REF DIP MOM     = ",3(1PD17.8,X)," au")'    ) dipgridsk(:)
                  write (outmat(i),'(  4X,"DIF REF DIP MOM = ",3(1PD17.8,X)," au")'    ) dipgrid(:)-dipgridsk(:)
               endif !! (iato.GT.0) then
            endif !! (ldip) then
         endif !! (lskp) then

      enddo !! i = 1,2

! >> (alpha and beta) Molecular orbital contributions
      label(:) = ''
      label    = '(6X,"N(",A'//int2str(ceiling(log10(dble(nom))))//',") = ",F12.6," au")'
      if (uhf) then
         write (ifuk,'(/,4X,"ALPHA ORBITAL CONTRIBUTIONS:")')
         do i=1,ialfa
            write (ifuk,label) int2str(matom(i)),poom(i)*pel(i)
         enddo !! j=1,nom
         write (ifuk,'(/,4X,"BETA  ORBITAL CONTRIBUTIONS:")')
         do i=ialfa+1,nom
            write (ifuk,label) int2str(matom(i)),poom(i)*pel(i)
         enddo !! j=1,nom
      else
         write (ifuk,'(/,4X,"ORBITAL CONTRIBUTIONS:")')
         do i=1,nom
            write (ifuk,label) int2str(matom(i)),poom(i)*pel(i)
         enddo !! j=1,nom
      endif !! (uhf) then

! >> Sigma/Pi separation
      if (sigpi.AND..NOT.lpi) then
!$omp parallel default(none) shared(tmp,poom,nompi,npi,pel)
!$omp   workshare
         tmp = sum(poom(nompi(:npi))*pel(nompi(:npi)))
!$omp   end workshare
!$omp end parallel
         if (uhf) then
!$omp parallel default(none) shared(tmpaa,tmpbb,ialfa,poom,nompi,npi,pel)
!$omp   workshare
            tmpaa = sum(poom(nompi(:npi))*pel(nompi(:npi)),MASK=nompi(:npi).LE.ialfa)
            tmpbb = sum(poom(nompi(:npi))*pel(nompi(:npi)),MASK=nompi(:npi).GT.ialfa)
!$omp   end workshare
!$omp end parallel
         endif !! (uhf) then
         do i = 1,2
            write (outmat(i),'(/,4X,"SIGMA/PI CONTRIBUTIONS:",2X,I4," molecular orbitals with PI symmetry", &
                              &/,4X,"SIGMA",F10.5,5X,"PI",F10.5)') npi,popul-tmp,tmp
         enddo !! i = 1,2
         if (uhf) then
            do i = 1,2
               write (outmat(i),'(6X,"*",1X,A," alpha spin orbitals with PI symmetry", &
                                 &/,8X,"ALPHA SIGMA",F10.5,5X,"ALPHA PI",F10.5)')      &
                                                  int2str(count(nompi(:npi).LE.ialfa)),&
                                                  tmpa-tmpaa                          ,&
                                                  tmpaa
               write (outmat(i),'(6X,"*",1X,A," beta  spin orbitals with PI symmetry", &
                                 &/,8X,"BETA  SIGMA",F10.5,5X,"BETA  PI",F10.5)')      &
                                                  int2str(count(nompi(:npi).GT.ialfa)),&
                                                  tmpb-tmpbb                          ,&
                                                  tmpbb
            enddo !! i = 1,2
         endif !! (uhf) then
         if (allocated(nompi)) deallocate(nompi)
      else if (sigpi) then
         write (ifuk,'(/,X,">> SIGMA/PI CONTRIBUTIONS skipped, only working with pi density as requested",/)') 
         if (allocated(nompi)) deallocate(nompi)
      endif !! (sigpi) then

      if (allocated(matom)) deallocate(matom)
      if (allocated(poom))  deallocate(poom)
      if (allocated(pel))   deallocate(pel)
      if (allocated(coord)) deallocate(coord)

! >> Atomic overlap matrix
      if (laom) then
         write (ifuk,'(/,10X,"The Atomic Overlap Matrix",/,/,/)')
         do iom=1,nbac
            write (ifuk,'(8E20.12)') (s(jom+iom*(iom-1)/2),jom=1,iom)
         enddo !! iom=1,nom
         deallocate(s)
      endif !! (laom) then
    end subroutine print_and_deallocate

    subroutine shortprint_and_deallocate

      implicit none

      integer(kind = i4) :: i


      ctgrid   = -rhpt   * coord ! Original grid
      if (lskpr.OR.(igrid.EQ.0.AND.lskp)) ctgridsk = -rhptsk * coord
      do i = 1,2
         write (outmat(i),'(/,X,"RESULTS OF THE INTEGRATION")',advance='no')
         if (iato.GT.0) write (outmat(i),'(X,"(ATOM",X,A,")",/)'           ) int2str(iato)
         write (outmat(i),'(4X,"REF POP",11X,"=",3X,1PD21.12," au")') rhpt
         if (iato.GT.0) then
            write (outmat(i),'(     4X,"REF DIP MOM")')
            write (outmat(i),'(6X,"* INTRINSIC",7X ,"= ",3(1PD17.8,X)," au")'   ) dipgrid(:)
            write (outmat(i),'(6X,"* CHARGE TRANSFER = ",3(1PD17.8,X)," au")'   ) ctgrid(:)
            write (outmat(i),'(6X,"* TOTAL",11X    ,"= ",3(1PD17.8,X)," au")'   ) ctgrid(:) + dipgrid(:)
         else
            write (outmat(i),'(  4X,"REF DIP MOM = ",3(1PD17.8,X)," au")'       ) dipgrid(:)
            write (outmat(i),'(10X,"** NOTE: Atom not specified, intrinsic dip. mom. unavailable")')
         endif !! (iato.GT.0) then

! >> Comparison with skipping
         if (lskpr.OR.(igrid.EQ.0.AND.lskp)) then
            write (outmat(i),'(/,4X,"Considering SKIPPING:")'         )
            write (outmat(i),'(/,4X,"REF POP   = ",3X,1PD21.14," au")') rhptsk
            write (outmat(i),'(4X,"DIF REF POP = ",   1PD9.3  ," au")') rhpt-rhptsk

! >> Dipole moment with skipping
            if (ldip) then
               if (iato.GT.0) then
                  write (outmat(i),'(/,4X,"REF DIP MOM")')
                  write (outmat(i),'(6X,"* INTRINSIC",7X ,"= ",3(1PD17.8,X)," au")'    ) dipgridsk(:)
                  write (outmat(i),'(6X,"* DIF INTRINSIC",3X ,"= ",3(1PD17.8,X)," au")') dipgrid(:)-dipgridsk(:)
                  write (outmat(i),'(6X,"* CHARGE TRANSFER     = ",3(1PD17.8,X)," au")') ctgridsk(:)
                  write (outmat(i),'(6X,"* DIF CHARGE TRANSFER = ",3(1PD17.8,X)," au")') ctgrid(:)-ctgridsk(:)
                  write (outmat(i),'(6X,"* TOTAL",11X    ,"= ",3(1PD17.8,X)," au")'    ) dipgridsk(:)+ctgridsk(:)
                  write (outmat(i),'(6X,"* DIF TOTAL",7X ,"= ",3(1PD17.8,X)," au")'    ) dipgrid(:)+ctgrid(:) &
                                                                                        -ctgridsk(:)-dipgridsk(:)
               else
                  write (outmat(i),'(  4X,"REF DIP MOM     = ",3(1PD17.8,X)," au")'    ) dipgridsk(:)
                  write (outmat(i),'(  4X,"DIF REF DIP MOM = ",3(1PD17.8,X)," au")'    ) dipgrid(:)-dipgridsk(:)
               endif !! (iato.GT.0) then
            endif !! (ldip) then
         endif !! (lskp) then
      enddo !! i = 1,2

    end subroutine shortprint_and_deallocate
  
end module output
