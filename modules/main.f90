! ** modules/main.f90 >> Main module where default values and arrays are set + Subroutines of general use in fukui program
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

module main

  use,intrinsic :: iso_fortran_env,   only: output_unit,real64,int32

  implicit none

  integer             ,parameter,public  :: dp        = real64      &
                                          , i4        = int32       &
                                          , ou        = output_unit

!! Default parameters in program
!! * DEFAULT CUT-OFF VALUES
  real     (kind = dp),parameter,public  :: ctffddef  =      18._dp & !! * Skipping gaussian functions (max distance in au for further primitives).
                                          , ctffgdef  = 1.0E-10     & !! * Skipping points (neglibible gaussian quadrature weights).
                                          , ctffrdef  = 5.0E-09       !! * Skipping points (negligible reference electron dens.).
!! *         CUT-OFF VARS
  real     (kind = dp)          ,public  :: ctffg = 0._dp,ctffr = 0._dp,ctffd
!! * MAX DEFAULT DIMENSIONS                             7654321
  integer  (kind = i4),parameter,public  :: natomxdef =      50     & !! * max # atoms.
                                          , nommxdef  =     250     & !! * max # molecular orbitals.
                                          , npmxdef   = 1000000     & !! * max # points.
                                          , nprimxdef =    1000       !! * max # primitive functions.
!!                                                      7654321
!! * MAX         DIMENSIONS                                    
  integer  (kind = i4)          ,public  :: natomx    = natomxdef   &
                                          , nommx     = nommxdef    &
                                          , npmx      = npmxdef     &
                                          , nprimx    = nprimxdef
!! * PROGRAM VERSION
  character(len  =  6),parameter,public  :: version   = '231003'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer  (kind = i4)          ,public  :: itype     =   0          !! Types of calculations:
                                                                     !!      xx
                                                                     !!      ||-> AOM
                                                                     !!      |--> Dip
                                                                     !!
                                                                     !!     bin-dec
                                                                     !!   *  00-  0: e density
                                                                     !!   *  01-  1: e density + AOM
                                                                     !!   *  10-  2: e density + dip moment
                                                                     !!   *  11-  3: e density + dip moment + AOM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!0123456789012345678901234567890
  character(len  = 22),parameter,dimension(0:3),public &
                                         :: txttype   = (/ 'e dens                ' &
                                                         , 'e dens + AOM          ' &
                                                         , 'e dens + dip mom      ' &
                                                         , 'e dens + dip mom + AOM' /)

  integer  (kind = i4)          ,public  :: igrid     =   0          !! Grid format (stk=0 default, cube=1).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... File units
  integer  (kind = i4),parameter,public  :: iwfn      =  15         &
                                          , istk      =  16         &
                                          , ifuk      =  17         &
                                          , isom      =  18         &
                                          , iskp      =  19         &
                                          , iproc     =  20         &
                                          , imxdm     =  21         &
                                          , iaom      =  22         &
                                          , idip      =  23

! ... Options in program if file units exist + paralelization
  logical                       ,public  :: loptions(11)            &
                                          , runparal  = .TRUE.      & !!  1
                                          , lproc     = .FALSE.     & !!  2
                                          , sigpi     = .FALSE.     & !!  3
                                          , lsg       = .FALSE.     & !!  4
                                          , lskpd     = .FALSE.     & !!  5
                                          , lskp      = .FALSE.     & !!  6
                                          , lskpr     = .FALSE.     & !!  7
                                          , lpi       = .FALSE.     & !!  8
                                          , laom      = .FALSE.     & !!  9
                                          , ldip      = .FALSE.     & !! 10
                                          , lonlygrid = .FALSE.       !! 11
  equivalence                              (loptions( 1),runparal ) &
                                          ,(loptions( 2),lproc    ) &
                                          ,(loptions( 3),sigpi    ) &
                                          ,(loptions( 4),lsg      ) &
                                          ,(loptions( 5),lskpd    ) &
                                          ,(loptions( 6),lskp     ) &
                                          ,(loptions( 7),lskpr    ) &
                                          ,(loptions( 8),lpi      ) &
                                          ,(loptions( 9),laom     ) &
                                          ,(loptions(10),ldip     ) &
                                          ,(loptions(11),lonlygrid)

! ... # args., # procs., specified atom, default # part to divide grid (w/o parallelization)
  integer  (kind = i4)          ,public  :: nproc, iato = 0, ipartdef = 8

! ... Partial times
  real     (kind = dp),dimension(100)    :: ptime = 0._dp
  real     (kind = dp),public            :: wtim0,cptim0,totim0

! ... Dipole moment: Cartesian coordinates backup if atom is included, dipole moments from grid
  real     (kind = dp),allocatable,dimension(:),public &
                                         :: coord,dipgrid,dipgridsk

  contains

    integer (kind = i4) function nlines(input)
      implicit none
      integer (kind = i4),intent(in) :: input
      integer (kind = i4)            :: ios

      rewind(input)
      nlines = 0
      do 
         read (input,*,iostat=ios)
         if (ios.NE.0) exit
         nlines    = nlines   +1
      enddo 
      if (ios.NE.0) then
         if (is_iostat_end(ios)) then
            rewind(input)
            return !! correct way to finish
         else
            stop 'ERROR: End-Of-File state not acquired and something wrong with file!'
         endif !! (is_iostat_end(ios)) then
      else
         stop 'ERROR: End-Of-Loop state acquired with wrong iostat exit code?!'
      endif !! (ios.NE.0) then
    end function nlines

    function minusculas(texto)
!
! Based on function "lowercase" in libSUFR: libsufr.sourceforge.net
!
      implicit none
      character,intent(in) :: texto*(*)
      character            :: minusculas*(len(texto))
      integer              :: i,ich
      minusculas = texto
      do i = 1,len_trim(minusculas)
         ich = ichar(minusculas(i:i))
         if (ich.GE.65.AND.ich.LE.91) ich = ich+32
         minusculas(i:i) = char(ich)
      enddo !! i = 1,len_trim(minusculas)
    end function minusculas

    function int2str(n)
      implicit none
      integer(kind = i4), intent(in) :: n
      character :: int2str*(max(ceiling(log10(dble(abs(n)+1))),1) - (sign(1,n)-1)/2)  ! 0-9 -> 1; 10-99 -> 2; +1 if <0
      
      write(int2str,'(I0)') n
    end function int2str

    subroutine dbl2str(n,dec,iprint,avanza,post)
      implicit none
      real     (kind = dp),intent(in)          :: n
      integer  (kind = i4),intent(in)          :: dec
      integer  (kind = i4),intent(in),optional :: iprint
      logical             ,intent(in),optional :: avanza
      character :: text2print*(max(ceiling(                                    & !! Considering rounding
                                              log10(                           & !! e.g. 99.95, 99.94 with dec=1
                                                      abs(n)                   & !! ceiling(log10(n)) = 2
                                                      + 6._dp*10._dp**(-dec-1) & !! +.06 = 100.01,100.00 => 3, 2
                                                   )                           &
                                          ),1) - (sign(1,int(n,i4))-1)/2 + dec + 1)
      character                                :: formato*(5+max(ceiling(log10(dble(abs(dec)+1))),1))
      character(len  =  *),intent(in),optional :: post
      character(len  =  3)                     :: avanzatxt

      avanzatxt = 'YES'
      if (present(avanza)) then
         avanzatxt(:) = ''
         if (.NOT.avanza) then
            avanzatxt = 'NO'
         else
            avanzatxt = 'YES'
         endif !! (.NOT.avanza) then
      endif !! (present(avanza)) then
      write(formato,'(A,I0,A)') '(F0.',max(dec,0),')'
      write(text2print,trim(formato)) n
      if (text2print(1:1).EQ.'.') then
         text2print = '0'//trim(text2print)
      else if (text2print(1:2).EQ.'-.') then
         text2print = '-0'//trim(text2print(2:))
      endif !! (pepito(1:1).EQ.'.') then
      if (.NOT.present(post)) then
         write(ou,'(A)',advance=trim(avanzatxt)) text2print
         if (present(iprint)) write(iprint,'(A)',advance=trim(avanzatxt)) text2print
      else
         write(ou,'(A)',advance='NO') text2print
         if (present(iprint)) write(iprint,'(A)',advance='NO') text2print
         write(ou,'(A)',advance=trim(avanzatxt)) trim(post)
         if (present(iprint)) write(iprint,'(A)',advance=trim(avanzatxt)) trim(post)
      endif !! (.NOT.present(post)) then
    end subroutine dbl2str

    subroutine sudgfchk(text,iout,icde,irw,jjj)
      implicit none
      integer(kind = i4),intent(in ) :: iout,irw
      integer(kind = i4),intent(out) :: icde,jjj
      integer(kind = i4)             :: iii
      character                      :: line*80,text*17

      icde = 0
      jjj  = 0
      if (irw.EQ.0) rewind(iout)
      do
        read (iout,'(A)',iostat=iii) line
        if (iii.NE.0) then
          jjj = 1
          return
        endif
        if (index(line,text).GT.0) then
           exit
        endif
      enddo
      icde = 1

    end subroutine sudgfchk

    subroutine prepnumb(imos,moc    ,nome)
      implicit none
      character(len  = 7)             :: charint,charintd
      character(len  = *),intent(out) :: nome
      integer  (kind = 4),intent(in ) :: imos
      integer  (kind = 4)             :: i,j,n
      integer  (kind = 4),intent(in ) :: moc(imos)

      nome(:) = ''
      j = 0
      n = 0
      do i = 1,imos
         if (moc(i)-j.NE.n) then
            n = moc(i)
            write(charint ,'(I7)') n
            if (i.NE.1) then
               if (j.EQ.1) then
                  nome = nome(1:len(trim(nome)))//','//trim(charint (verify(charint ,' '):7))
               else
                  write(charintd,'(I7)') moc(i-1)
                  nome = nome(1:len(trim(nome)))//'-'//trim(charintd(verify(charintd,' '):7))//','// &
                                                       trim(charint (verify(charint ,' '):7))
               endif !! (j.EQ.1) then
            else
                  nome = nome(1:len(trim(nome)))//trim(charint(verify(charint,' '):7))
            endif !! (i.NE.1) then
            j = 1
         else
            j = j+1
            if (i.EQ.imos) then
                  write(charint ,'(I7)') moc(i)
                  nome = nome(1:len(trim(nome)))//'-'//trim(charint (verify(charint ,' '):7))  
            endif !! (i.EQ.imos) then
         endif !! (moc(i)-j.NE.n) then
      enddo !! i = 1,imos

    end subroutine prepnumb

    subroutine print_main(texto)
      implicit none
      character(len  =  *),intent(in)   :: texto
      write(ou  ,'(X,">>",X,A)') trim(texto)
      write(ifuk,'(X,">>",X,A)') trim(texto)
    end subroutine print_main

    subroutine print_second(texto)
      implicit none
      character(len  =  *),intent(in)   :: texto
      write(ou  ,'(11X,"*",X,A)') trim(texto)
      write(ifuk,'(11X,"*",X,A)') trim(texto)
    end subroutine print_second

    subroutine print_rdarr(texto)
      implicit none
      integer  (kind = i4)              :: i
      character(len  =  *),intent(in)   :: texto
      character(len  = 30)              :: innerlabel
      i = 40-len_trim(texto)
      if (i.GT.0) then
         innerlabel(:) = ''
         innerlabel    = '(11X,"*",X,A,X,'//int2str(i)//'("."),X)'
         write(ou  ,innerlabel,advance='no') trim(texto)
         write(ifuk,innerlabel,advance='no') trim(texto)
      else
         write(ou  ,'(11X,"* ",A,":")',advance='no') trim(texto)
         write(ifuk,'(11X,"* ",A,":")',advance='no') trim(texto)
      endif !! (i.GT.0) then
    end subroutine print_rdarr

    subroutine partial_time(n,text,lfinal)

      use omp_lib
  
      implicit none

      integer  ( kind=i4  )            :: n,nlength,nproccalc
      integer  ( kind=i4  ),parameter  :: ntotchar=30
      real     ( kind=dp  )            :: tmp
      character(  len=*   ),intent(in) :: text
      character(  len=100 )            :: label
      logical     ,optional,intent(in) :: lfinal
      logical                          :: lfinaldef=.FALSE.

      if (present(lfinal)) lfinaldef = lfinal
      if (.NOT.lfinaldef) then
         if (n.GT.size(ptime)) then
            write(ou,'(2(/),10X,"ERROR: ** Increase dimension of ptime **")')
            stop
         endif !! (n.GT.size(ptime)) thenb
         if (n.EQ.0) then
            wtim0 = omp_get_wtime()
         else
            ptime(n) = omp_get_wtime()-wtim0-sum(ptime(:n))
            nlength  = len_trim(text)
            if (nlength.GT.ntotchar) then
               write(ou,'(2(/),10X,"ERROR: ** Text too long to be printed with time **")')
               stop
            else              
                label(:) = '' 
                label    = '(/,5X,"[TT]",X,A,1X,'//int2str(43-nlength)//'("."),1X)'
                write(ou  ,label,advance='NO') trim(text)
                write(ifuk,label,advance='NO') trim(text)
                call dbl2str(ptime(n),2,iprint=ifuk,post=' seconds',avanza=.FALSE.)
                write(ou  ,'(/)')
                write(ifuk,'(/)')
            endif !! (nlength.GT.ntotchar) then
         endif !! (n.EQ.0) then
         n = n+1
      else 
                label(:) = '' 
                label    = '(/,5X,"[TT]",X,"TOTAL ELAPSED TIME",1X,25("."),1X)'
                write(ou  ,label,advance='NO')
                write(ifuk,label,advance='NO')
                tmp = sum(ptime(:n-1))
                call dbl2str(tmp,2,iprint=ifuk,post=' seconds',avanza=.TRUE.)
                label    = '(10X,"PERFORMANCE",1X,32("."),1X)'
                write(ou  ,label,advance='NO')
                write(ifuk,label,advance='NO')
                nproccalc = nproc
                if (nproc.EQ.0) nproccalc = 1
                call dbl2str(100._dp*(totim0-cptim0)/tmp/nproccalc,2,iprint=ifuk,post='%',avanza=.TRUE.)
                write(ou  ,'(/)')
                write(ifuk,'(/)')
      endif !! (.NOT.lfinaldef) then

    end subroutine partial_time

end module main
