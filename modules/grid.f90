! ** modules/grid.f90 >> Main module to process (stk or cube) grid files
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

module grid

  use,intrinsic :: iso_fortran_env,   only: output_unit,real64,int32
  use           :: main           ,   only: istk,lskp,lskpr,igrid,ifuk

  implicit none

  integer                   ,parameter,private:: dp = real64
  integer                   ,parameter,private:: i4 = int32       &
                                               , ou = output_unit 

  type,public :: xyzw
    real (kind = dp)                                 :: x,y,z,w
  end type       xyzw

  integer(kind = i4)                         ,public :: nptos
  logical                                    ,public :: is_stk_bin
  type   (xyzw     ),allocatable,dimension(:),public :: points
  real   (kind = dp)                                 :: rhpt,rhptsk

  contains

  subroutine openstk(inp,filename,is_grid_bin,iter)
    implicit none
    integer  (kind = i4),intent(in ) :: inp
    integer  (kind = i4)             :: i
    character(len  =  *),intent(in ) :: filename
    character(len  =  1)             :: txttmp
    logical             ,intent(out) :: is_grid_bin
    integer  (kind = i4),intent(out) :: iter
!   See:
!   https://stackoverflow.com/questions/38651487/how-to-detect-whether-a-file-is-formatted-or-unformatted
!   to understand the following part of the code.
    i          = 0
    iter       = 0
    is_stk_bin = .TRUE.
    open (unit=inp,file=trim(filename),status='old',form='unformatted',recl=1)
    do while ((i.EQ.0).AND.is_stk_bin)
       read (inp,iostat=i) txttmp
       is_stk_bin = is_stk_bin.AND.(iachar(txttmp)<=127)
       iter       = iter+1
    enddo !! while ((iii.EQ.0).AND.is_stk_bin)
    if (is_stk_bin) then
       close(inp)
       open (unit=inp,file=trim(filename),status='old',form=  'formatted')
    else
       close(inp)
       open (unit=inp,file=trim(filename),status='old',form='unformatted')
    endif !! (is_stk_bin) then
    is_grid_bin = is_stk_bin
  end subroutine openstk

  integer function numbnptos(inp)
    implicit none
    integer(kind = i4),intent(in) :: inp
    integer(kind = i4)            :: i
    real   (kind = dp)            :: a
    character                     :: txt

    rewind(inp)
    numbnptos = 0
    if (is_stk_bin) then
       read (inp,*) txt
       do 
          read (inp,*,iostat=i) a
          if (i.NE.0) exit
          numbnptos = numbnptos+1
       enddo 
    else
       do 
          read (inp,iostat=i) a
          if (i.NE.0) exit
          numbnptos = numbnptos+1
       enddo 
    endif !! (is_stk_bin) then

  end function numbnptos

  integer function nptcube(inp)
    implicit none
    integer    ,intent(in) :: inp
    integer(i4)            :: ipar(4)
    rewind(inp)
    read (inp,'(1(/))')
    read (inp,'(4(I5,/))') ipar(:)
    nptcube = product(ipar(2:))
  end function

  subroutine readgrid
    use :: main                , only: int2str,npmx,print_second
    implicit none
    integer  (kind = i4)            :: nptst

    select case (igrid)
       case(0)
         nptos  = numbnptos(istk)
       case(1)
         if (.NOT.is_stk_bin) stop 'ERROR: the cube file is unformatted!'
         nptos  = nptcube(istk)
       case default
         stop 'Unknown grid format'
    end select !! case (igrid)
    call print_second(int2str(nptos)//" points (total)")
    if (npmx  .LT.nptos) then
       write(ou,'(" ** Increase    npmx **",A,"=old value < nptos=",A)') int2str(npmx),int2str(nptos)
       stop ' ** Edit source code in main module or create a file with extension .mxal following the manual **'
    endif !! (npmx  .LT.nptos) then
    allocate(points(nptos))                 !! R1
    select case(igrid)
      case(0)           !! STK  file
         if (lskp.OR.lskpr) then
            call    le1_stksk(rhpt,rhptsk,nptst)
            write(ou  ,'(10X," * ",A," of ",A," points were considered -> ",F6.2,"% skipped")') &
                           int2str(nptos),int2str(nptst),100*dble(nptst-nptos)/dble(nptst)
            write(ifuk,'(10X," * ",A," of ",A," points were considered -> ",F6.2,"% skipped")') &
                           int2str(nptos),int2str(nptst),100*dble(nptst-nptos)/dble(nptst)
         else
            call le1_stk(rhpt)
            call print_second(int2str(nptos)//" points read")
         endif !! (lskp.OR.lskpr) then
      case(1)           !! CUBE file
         call readcube(rhpt,rhptsk,nptst)
! NOTE: ctffg is incompatible with cube files. The weight of each point
! corresponds to the scalar triple product of the parallelepiped vectors
         write(ou  ,'(10X," * ",A," of ",A," points were considered -> ",F6.2,"% skipped")') &
                        int2str(nptos),int2str(nptst),100*dble(nptst-nptos)/dble(nptst)
         write(ifuk,'(10X," * ",A," of ",A," points were considered -> ",F6.2,"% skipped")') &
                        int2str(nptos),int2str(nptst),100*dble(nptst-nptos)/dble(nptst)
      case default
         stop 'Unknown grid file format'
    end select !! case(igrid)

  end subroutine readgrid

  subroutine le1_stk(totpop)
    use :: main                 ,only : int2str,minusculas,print_second,print_rdarr
    implicit none
    integer  (kind = i4)             :: i,j
    real     (kind = dp),intent(out) :: totpop
    real     (kind = dp)             :: rhp
    character(len  = 51)             :: customfmt

    rewind(istk)
    totpop       = 0._dp
    customfmt(:) = ''
    if (is_stk_bin) then
       read (istk,'(A)',iostat=i) customfmt
       if (i.NE.0) stop 'ERROR while reading stk file. Abnormal format line'
       if (len_trim(customfmt).GT.len(customfmt)-1) then
          write(ou,*) 'ERROR: Format line exceeded the maximum length: ',int2str(len(customfmt)-1)," characters"
          stop
       endif !! (len_trim(customfmt).GT.len(customfmt)-1)
       if (index(minusculas(customfmt),'free').GT.0) then
          call print_second("Free-format ordinary text stk file")
          do j = 1,nptos
             read (istk,*,iostat=i,err=2) points(j)%x,points(j)%y,points(j)%z,points(j)%w,rhp
             totpop = totpop + points(j)%w*rhp
          enddo 
       else
          call print_rdarr("Custom-format ordinary text stk file: ")
          write(ou,'(A,A,A)') &
                    '"',customfmt(                                             &
                                    index(trim(customfmt),' ',.TRUE.)+1        &
                                              :                                &
                                    len_trim(customfmt)                        &
                                 ),'"'
          if (index(                                                           &
                        customfmt(                                             &
                                    index(trim(customfmt),' ',.TRUE.)+1        &
                                              :                                &
                                    len_trim(customfmt)                        &
                                 )                                             &
                        ,'('                                                   &
                   ).LE.0)                                                     &
             stop 'ERROR while reading stk file. Add initial and final parentheses "( )" to the custom format'
          if (index(                                                           &
                        customfmt(                                             &
                                    index(trim(customfmt),' ',.TRUE.)+1        &
                                              :                                &
                                    len_trim(customfmt)                        &
                                 )                                             &
                        ,')'                                                   &
                   ).LE.0)                                                     &
             stop 'ERROR while reading stk file. Add final parenthesis "( )" to the custom format'
          do j = 1,nptos
             read (istk,trim(customfmt),end=2,iostat=i) points(j)%x,points(j)%y,points(j)%z,points(j)%w,rhp
             totpop = totpop + points(j)%w*rhp
          enddo 
       endif !! (index(minusculas(customfmt),'free').GT.0) then
    else
       call print_second("Binary stk file")
       do j = 1,nptos
          read (istk,end=2,iostat=i)  points(j)%x,points(j)%y,points(j)%z,points(j)%w,rhp
          totpop  = totpop + points(j)%w*rhp
       enddo 
    endif !! (lfmtstk) then
2   if (i.NE.0) stop ' ERROR: Unknown problem with stk file while reading'
  end subroutine le1_stk

  subroutine le1_stksk(totpop,totpopsk,nptst)
    use :: main,                  only: int2str,minusculas,cutoff=>ctffg,cutoffrho=>ctffr
    implicit none
    integer  (kind = i4)             :: nptst,i,j
    real     (kind = dp),intent(out) :: totpop,totpopsk
    real     (kind = dp)             :: rhp
    type     (xyzw     ),allocatable :: tmp(:)
    character(len  = 51)             :: customfmt

    rewind(istk)
    nptst    = nptos
    nptos    = 0
    totpop   = 0._dp
    totpopsk = 0._dp
    i        = 0
    if (is_stk_bin) then
       read (istk,'(A)',iostat=i) customfmt !! With free format, the code jumps next line
       if (i.NE.0) stop 'ERROR while reading stk file. Abnormal format line'
       if (len_trim(customfmt).GT.len(customfmt)-1) then
          write(ou,*) 'ERROR: Format line exceeded the maximum length: ',int2str(len(customfmt)-1)," characters"
          stop
       endif !! (len_trim(customfmt).GT.len(customfmt)-1)
       i = 0
       if (index(minusculas(customfmt),'free').GT.0) then
          write(ou,'(10X,"* Free-format ordinary text stk file")')
          do j = 1,nptst
             nptos = nptos+1
             read (istk,*,end=2,iostat=i) points(nptos)%x,points(nptos)%y,points(nptos)%z,points(nptos)%w,rhp
             totpop = totpop + points(nptos)%w*rhp
             if (points(nptos)%w.LT.cutoff.OR.rhp.LT.cutoffrho) then
                nptos = nptos-1
             else
                totpopsk  = totpopsk + points(nptos)%w*rhp
             endif !! (points(nptos)w.LT.cutoff) then
          enddo 
       else
          write(ou,'(10X," * Custom-format ordinary text stk file:",X,A,A,A)') &
                    '"',customfmt(                                             &
                                    index(trim(customfmt),' ',.TRUE.)+1        &
                                              :                                &
                                    len_trim(customfmt)                        &
                                 ),'"'
          if (index(                                                           &
                        customfmt(                                             &
                                    index(trim(customfmt),' ',.TRUE.)+1        &
                                              :                                &
                                    len_trim(customfmt)                        &
                                 )                                             &
                        ,'('                                                   &
                   ).LE.0)                                                     &
             stop 'ERROR while reading stk file. Add initial and final parentheses "( )" to the custom format'
          if (index(                                                           &
                        customfmt(                                             &
                                    index(trim(customfmt),' ',.TRUE.)+1        &
                                              :                                &
                                    len_trim(customfmt)                        &
                                 )                                             &
                        ,')'                                                   &
                   ).LE.0)                                                     &
             stop 'ERROR while reading stk file. Add final parenthesis "( )" to the custom format'
          do j = 1,nptst
             nptos = nptos+1
             read (istk,trim(customfmt),end=2,iostat=i) &
                    points(nptos)%x,points(nptos)%y,points(nptos)%z,points(nptos)%w,rhp
             totpop = totpop + points(nptos)%w*rhp
             if (points(nptos)%w.LT.cutoff.OR.rhp.LT.cutoffrho) then
                nptos = nptos-1
             else
                totpopsk = totpopsk + points(nptos)%w*rhp
             endif !! (points(nptos)%w.LT.cutoff) then
          enddo 
       endif !! (index(minusculas(customfmt),'free').GT.0) then
    else
       write(ou,'(10X," * Binary stk file")')
       do j = 1,nptst
          nptos = nptos+1
          read (istk,end=2,iostat=i) points(nptos)%x,points(nptos)%y,points(nptos)%z,points(nptos)%w,rhp
          totpop  = totpop + points(nptos)%w*rhp
          if (points(nptos)%w.LT.cutoff.OR.rhp.LT.cutoffrho) then
             nptos = nptos-1
          else
             totpopsk = totpopsk + points(nptos)%w*rhp
          endif !! (points(nptos)%w.LT.cutoff) then
       enddo 
    endif !! (is_stk_bin) then
 2  if (i.GT.0) stop ' ** Some problem with stk file **'
 ! Reducing arrays size. increasing memory fragmention?
    call move_alloc(points,tmp)
    allocate(points,source=tmp(:nptos))
    deallocate(tmp)

  end subroutine le1_stksk

  subroutine readcube(totpop,totpopsk,nptst)
    use :: main,                              only: int2str,dbl2str,minusculas,cutoffrho=>ctffr,print_rdarr
    implicit none
    integer  (kind= i4)                          :: nptst
    integer  (kind= i4)                          :: i1,i2,i3,n
    real     (kind= dp),intent(out)              :: totpop,totpopsk
    integer  (kind= i4)                          :: ipar(0:3),natoms,nslow,ninter,nfast
    real     (kind= dp)                          :: vectors(0:3,3),scalefactor
    real     (kind= dp),allocatable,dimension(:) :: tmparray
    character( len=100)                          :: tmptxt

    equivalence(ipar(0),natoms)
    equivalence(ipar(1),nslow )
    equivalence(ipar(2),ninter)
    equivalence(ipar(3),nfast )

    rewind(istk)
    read (istk,'(1(/))')
    do i1=0,3
       read (istk,'(I5,3(F12.6))') ipar(i1),vectors(i1,:)
    enddo !! i=0,3
    scalefactor = vectors(1,1)*vectors(2,2)*vectors(3,3) &
                + vectors(2,1)*vectors(3,2)*vectors(1,3) &
                + vectors(3,1)*vectors(1,2)*vectors(2,3) &
                - vectors(3,1)*vectors(2,2)*vectors(1,3) &
                - vectors(3,2)*vectors(2,3)*vectors(1,1) &
                - vectors(3,3)*vectors(2,1)*vectors(1,2)
    nptst       = nptos
    call print_rdarr("# atoms")
    write(ou  ,'(A)') int2str(natoms)
    write(ifuk,'(A)') int2str(natoms)
    call print_rdarr("# slow points")
    write(ou  ,'(A)') int2str(nslow)
    write(ifuk,'(A)') int2str(nslow)
    call print_rdarr("# inter. points")
    write(ou  ,'(A)') int2str(ninter)
    write(ifuk,'(A)') int2str(ninter)
    call print_rdarr("# fast points")
    write(ou  ,'(A)') int2str(nfast)
    write(ifuk,'(A)') int2str(nfast)
    call print_rdarr("Scale factor")
    call dbl2str(scalefactor,6,ifuk)
    call print_rdarr("Total # points")
    write(ou  ,'(A)') int2str(nptst)
    write(ifuk,'(A)') int2str(nptst)

    do i1=1,natoms
       read (istk,*)
    enddo !! i=1,natoms
    allocate(tmparray(nfast))
    n          = 0
    nptos      = 0
    totpop     = 0._dp
    totpopsk   = 0._dp
!! Check where cube file comes from
    tmptxt(:)  = ''
    call print_rdarr("Format type")
    read (istk,'(A)') tmptxt
    backspace(istk)
    select case(len_trim(tmptxt))
      case(78)
        write(ou  ,'(A)') 'G09 Cubegen or compatible'
        write(ifuk,'(A)') 'G09 Cubegen or compatible'
      case(84)
        write(ou  ,'(A)') 'Multiwfn 3.7 or compatible'
        write(ifuk,'(A)') 'Multiwfn 3.7 or compatible'
      case default
        write(ou  ,'(A)') 'Unknown'
        write(ifuk,'(A)') 'Unknown'
    end select
    if (lskpr) then
       do i1 = 1,nslow
          do i2 = 1,ninter
             read (istk,*) tmparray(:)
!$omp parallel default (none) shared ( totpop,totpopsk,tmparray,cutoffrho,nptos )
  !$omp workshare
             totpop   = totpop   +   sum(tmparray)
             totpopsk = totpopsk +   sum(tmparray,mask=tmparray.GT.cutoffrho)
             nptos    = nptos    + count(tmparray.GT.cutoffrho)
  !$omp end workshare
!$omp end parallel
             do i3 = 1,nfast
                if (tmparray(i3).GT.cutoffrho) then
                   n     = n+1
                   points(n)%x = vectors(0,1)+(i1-1)*vectors(1,1)+(i2-1)*vectors(2,1)+(i3-1)*vectors(3,1)
                   points(n)%y = vectors(0,2)+(i1-1)*vectors(1,2)+(i2-1)*vectors(2,2)+(i3-1)*vectors(3,2)
                   points(n)%z = vectors(0,3)+(i1-1)*vectors(1,3)+(i2-1)*vectors(2,3)+(i3-1)*vectors(3,3)
                endif !! (tmparray(i3).GT.cutoffrho) then
             enddo !! I3 = 1,nfast
          enddo !! I2 = 1,ninter
       enddo !! I1 = 1,nslow
    else
       nptos   = 0
       totpop  = 0._dp
       if (nptst.NE.nslow*ninter*nfast) stop 'ERROR: initial # points != # points read.'
       do i1 = 1,nslow
          do i2 = 1,ninter
             read (istk,'(6(E13.5))') tmparray(:)
             totpop = totpop + sum(tmparray)
             do i3 = 1,nfast
                nptos     = nptos + 1
                points(nptos)%x = vectors(0,1)+(i1-1)*vectors(1,1)+(i2-1)*vectors(2,1)+(i3-1)*vectors(3,1)
                points(nptos)%y = vectors(0,2)+(i1-1)*vectors(1,2)+(i2-1)*vectors(2,2)+(i3-1)*vectors(3,2)
                points(nptos)%z = vectors(0,3)+(i1-1)*vectors(1,3)+(i2-1)*vectors(2,3)+(i3-1)*vectors(3,3)
             enddo !! I3 = 1,nfast
          enddo !! I2 = 1,ninter
       enddo !! I1 = 1,nslow
    endif !! (lskpr) then
    points(:nptos)%w = scalefactor
    totpop     = totpop * scalefactor
    if (lskpr) totpopsk  = totpopsk*scalefactor
    deallocate(tmparray)
  end subroutine

end module grid
