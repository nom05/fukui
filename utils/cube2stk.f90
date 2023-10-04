! ** utils/cube2stk.f90 >> Utility to transform cube files into stk files
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


program cube2stk

  use,intrinsic :: iso_fortran_env, only: output_unit,real64,int32
  implicit none

  integer,parameter                          :: dp = real64
  integer,parameter                          :: i4 = int32

  character(len=99)                          :: filename,fileroot,output
  character(len=81)                          :: name
  logical                                    :: lexist
  integer(i4)                                :: i,i1,i2,i3,j,iarg &
                                              , icube = 14        & !! Ref.   cube file
                                              , iout  = 15        & !! Output cube file
                                              , iato  =  0
  integer(i4)                                :: natoms,nslow,ninter,nfast,ipar(0:3)
  integer(i4),dimension(:      ),allocatable :: iatnumb
  real(dp)                                   :: vectors(0:3,3),scalefactor
  real(dp)   ,dimension(:,:    ),allocatable :: cooratom
  real(dp)   ,dimension(:,:,:  ),allocatable :: densgrid
  real(dp)   ,dimension(:,:,:,:),allocatable :: coord
  real(dp),parameter                         :: bohr2a=.5291772086_dp

  equivalence(ipar(0),natoms)
  equivalence(ipar(1),nslow )
  equivalence(ipar(2),ninter)
  equivalence(ipar(3),nfast )

  iarg    = iargc() !! Number of arguments.
  name(:) = ''

  if (iarg.NE.2) then
     call getarg(0,name)
     write(output_unit,'(" ** HELP:",1X,A,1X,"cubefile.cube atom_number",4X,"(w/ or w/o extension, number as integer)",/)') &
          trim(name)
     stop
  endif !! (iarg.NE.2) then
  call getarg(1,name)
  if (len_trim(name).EQ.81) stop 'File name too long'
  inquire(file=trim(name),exist=lexist)
  filename(:) = ''
  if (lexist) then ! filename*
     fileroot(:) = ''
     filename    = trim(name)
     if (index(name,'.',.TRUE.).GT.0) then
        fileroot = name(:index(name,'.',.TRUE.)-1)
     else
        fileroot = trim(name)
     endif !! (index(name,'.',.TRUE.).GT.0) then
  else
     fileroot    = trim(name)
     filename    = trim(name)
     inquire(file=trim(filename)//'.cub',exist=lexist)
     if (lexist) then ! filename.cub
        filename = trim(filename)//'.cub'
     else
        inquire(file=trim(filename)//'.cube',exist=lexist)
        if (lexist) then ! filename.cube
           filename = trim(filename)//'.cube'
        endif !! (lexist) then ! filename.cube
     endif !! (lexist) then ! filename.cub
  endif !! (lexist) then ! filename*
  if (.NOT.lexist) then
     write(output_unit,'(2x,"Requested file does not exist: ",a)') trim(name)
     stop 'Fatal error!!'
  endif !! (.NOT.lexist) then
  name(:) = ''
  if (iarg.EQ.2) then
     call getarg(2,name)
     read (name(verify(name,' '):verify(name,' ',.TRUE.)),'(I5)') iato
  endif !! (iarg.GT.1) then
  open(unit=icube,file=trim(filename),status='OLD')
  output(:) = ''
  output    = trim(fileroot)//'.stk'
  inquire(file=trim(output),exist=lexist)
  if (lexist) stop '.stk file exists and will not be overwritten. Remove it by hand.'
  open(unit=iout,file=trim(output),status='NEW',form='UNFORMATTED')

! Collecting data:
  rewind(icube)
  read (icube,'(1(/))')
  do i=0,3
     read (icube,'(I5,3(F12.6))') ipar(i),vectors(i,:)
  enddo !! i=0,3
  scalefactor = vectors(1,1)*vectors(2,2)*vectors(3,3)
  write(output_unit,'(1X,"Requested atom ...",1X,I5   )') iato
  write(output_unit,'(1X,"# atoms ..........",1X,I5   ,/,&
                     &1X,"# slow points ....",1X,I5   ,/,&
                     &1X,"# inter.points ...",1X,I5   ,/,&
                     &1X,"# fast points ....",1X,I5   ,/,&
                     &1X,"Scale factor .....",1X,F14.8,/,&
                     &1X,"Total # points ...",1X,I9   )') natoms,nslow,ninter,nfast,scalefactor,nslow*ninter*nfast

  allocate(densgrid(nfast,ninter,nslow),coord(nfast,ninter,nslow,3))
  allocate(iatnumb(natoms),cooratom(natoms,3))

  do i=1,natoms
     read (icube, '(I5,12X,3(F12.6))') iatnumb(i),(cooratom(i,j),j=1,3)
  enddo !! i=1,natoms
  print *,'---------------------------'
  print *,'First 5 atomic coordinates:'
  i1 = min(natoms,5)
  do i = 1,i1
     write(output_unit,'(I2,I3,3(F12.6))') i,iatnumb(i),(cooratom(i,j),j=1,3)
  enddo !! i = 1,i1
  print *,'---------------------------'

  vectors = vectors * bohr2a
  write(output_unit,'(   "Reading d.val. ...")',advance='no')
  do I1 = 1,nslow
     Do I2 = 1,ninter
        read (icube,'(6E13.5)') (densgrid(I3,I2,I1),I3=1,nfast)
        do I3 = 1,nfast
           do i=1,3
              coord(I3,I2,I1,i)=vectors(0,i)+(I1-1)*vectors(1,i)+(I2-1)*vectors(2,i)+(I3-1)*vectors(3,i)
           enddo !! i=1,3
        enddo !! I3 = 1,nfast
     enddo !! I2 = 1,ninter
  enddo !! I1 = 1,nslow
  print *,'Finished correctly'
  close(icube)

  write(output_unit,'(   "Writing stk file .")',advance='no')
  do I1 = 1,nslow
     Do I2 = 1,ninter
        do I3 = 1,nfast
!          write(iout ,'(5(1PE13.5))') coord(I3,I2,I1,:3),1._dp,densgrid(I3,I2,I1) !! For testing
           write(iout) coord(I3,I2,I1,:3),scalefactor,densgrid(I3,I2,I1)
        enddo !! I3 = 1,nfast
     enddo !! I2 = 1,ninter
  enddo !! I1 = 1,nslow
  print *,'Done'
  close(iout)

end program
