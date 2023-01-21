! ** utils/bstk_astk.f90 >> Utility to transform stk files from ASCII to binary and vice versa.
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


program stk
!
!   Programa para convertir stk's binarios en stk's ascii e viceversa
!
      implicit real*8 (a-h,o-z) 
      character fileg*1000,at*5,filestu*100,filestd*100
      logical kin
!
! Se kin é true, le binary e escribe ASCII
! Se kin é false, le ASCII e escribe binary
!
! Input de datos:
!
      at(:) = ''
      call getarg(2,at)
      if (at .EQ. '') then
         write (*,*) 'binary DATA or ASCII data? data/ascii?: '
         read (*,*) at
      endif
      if (at(1:4).EQ.'DATA'.OR.at(1:4).EQ.'data'.OR.at(1:5).EQ.'ASCII'.OR.at(1:5).EQ.'ascii') goto 1
      stop 'Format not correctly selected'
    1 if (at(1:4).EQ.'DATA'.or.at(1:4).EQ.'data') then
         kin=.TRUE.
      else
         kin=.FALSE.
      endif

      call getarg(1,fileg)
      if (fileg .EQ. '') then
         write (*,*) 'File (.stk) with stk data (without extension)): '
         read (*,'(a)') fileg
      endif
      filestu = fileg(1:index(fileg,' ')-1)//'.stk'
      filestd = fileg(1:index(fileg,' ')-1)//'.stkd'
      write (*,*) at,' .stk file: ', filestu
      istu=20
      istd=30
      if(kin) then
        open (istu, FILE=filestu, STATUS='UNKNOWN',form='unformatted')
        open (istd, FILE=filestd, STATUS='UNKNOWN')
      else
        open (istu, FILE=filestu, STATUS='UNKNOWN')
        open (istd, FILE=filestd, STATUS='UNKNOWN',form='unformatted')
      endif
!
! Obtendo o número de puntos,lendo e escribindo:
!
      write (*,*) kin
      call long(istu,istd,kin)

      stop
end program stk
!
! --------------------------------------------------------------
!
subroutine long(istu,istd,kin)
!
! Surrutina que obtén o número de liñas para así definir os arrays
! dinámicos
!
      implicit real*8 (a-h,o-z)
      logical kin

      if (kin) then
         i=1
    1    continue
         read (istu,end=2) xp,yp,zp,wp,rhp
         write(istd,'(5x,3(f15.8,1x),2(2x,e15.8))') xp,yp,zp,wp,rhp
         i=i+1
         goto 1
    2    continue
      else
         i=1
    3    continue
         read (istu,*,end=4) xp,yp,zp,wp,rhp
         write(istd) xp,yp,zp,wp,rhp
         i=i+1
         goto 3
    4    continue
      endif
      nptos=i-1
      write (*,'(I10," points")') nptos

end subroutine long
