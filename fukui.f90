! ** ./fukui.f90 >> Program to calculate electron densities from a set of points and calculate atomic overlap matrices 
!                   (read from a .stk or cube files) by means of the use of a .wfn file.
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

program fukui

! VERSIONS and CHANGELOG
!  * Version 1 (Mars 10 2007).
!  * Version 2 (December 1 2008) prints Sij atomic component. Based on the program desglose-e7.f.
!  * Version 151203. Restricted/unrestricted, it prints more info in the output,
!    sigma/pi or only pi calculations, skips points with a very low gaussian quadrature
!    weight (optional), skip points with a low reference electron density, reads max dims 
!    of arrays from a external file in the source code.
!  * Version 151207. Each subroutine/function is a file, CMake is hereinafter
!    employed.
!  * Version 160708: - New structure of code.
!                    - New subroutine 'goon'
!                    - OpenMP parallelization.
!                    - Serial and parallel are practically independent codes sharing
!                      common parts (thanks to CMake).
!                    - Dynamic definition of arrays (the cut-offs to define the max
!                      size of arrays will not be removed to avoid huge memory use 
!                      in small computers.
!                    - Help is in an external file.
!                    - AOM calc is an option.
! NOTE (http://stackoverflow.com/questions/13870564/gfortran-openmp-segmentation-fault-occurs-on-basic-do-loop):
! Be careful with huge arrays. When default arrays are used from the beginning, FORTRAN compilers automatically 
! placed in the static storage (known as the .bss segment) while small arrays are allocated on the stack. When you 
! switch OpenMP support on and in some cases with a normal compilations, no automatic static allocation is used and
! these arrays gets allocated on the stack of the routine. The default stack limits on the operative systems are
! very restrictive, hence the segmentation fault is observed. To avoid this, I have defined these arrays by means of
! allocatable variables.
!  * Version 161103: - Adapted to work with wfn file obtained by means of AIMAll.
!  * Version 211126: - New quicksort code is added from https://github.com/mcuntz/jams_fortran.git (MIT license)
!                    - Support to read .sg files (Sigma/Pi separation or MO specification). .sg files have more 
!                      priority than .som files (they are searched first)
!  * Version 211229: - fukui detects if the stk is binary or text. Utils updated.
!  * Version 211230: - cube file support
!  * Version 220322: - First f90 modular version
!                    - Quicksort dependency is removed from the code
!                    - Cube files support
!
!!! ... Define default parameters

  use           :: main                   ,  only: ou,i4,version,runparal,nproc,igrid        &
                                                 , iwfn,ifuk,iproc,isom,iskp,istk               &
                                                 , nlines,imxdm,natomx,nommx,npmx,nprimx        &
                                                 , sigpi,lsg,lproc,lpi,lskpd,lskp,lskpr,laom    &
                                                 , ctffd,ctffddef,ctffg,ctffgdef,ctffr,ctffrdef &
                                                 , sudgfchk,prepnumb,int2str,dbl2str            &
                                                 , print_rdarr,print_second,print_main          &
                                                 , partial_time,cptim0,totim0
  use           :: grid                   ,  only: openstk

  implicit none

! ... # args., # procs., specified atom
  integer  (kind = i4)                          :: iarg
! ... Temporary variables
  integer  (kind = i4)                          :: i,iato
  character(len  =  4)                          :: ext
  character(len  =  4),parameter,dimension(0:1) :: gridext=(/'stk ','cube'/)
! ... Logical variables
  logical                                       :: lexist,lfile2,lfmtstk,ltmp
! ... Time variables
  integer  (kind = i4)                          :: ntime = 0

! ... File names and so on
  character                                     :: file1*100,fproc*105,fwfn*104,file2*105,fstk*109,ffuk*209 &
                                                 , fmxdm*105,fsom*104,fskpd*104,fskp*104,fskpr*104,fpi*104  &
                                                 , faom*104

  interface
    subroutine help
    end subroutine help
    subroutine goon(i,j)
      use :: main,only: i4
      integer(kind = i4) :: i,j
    end subroutine goon
    subroutine init_par
    end subroutine init_par
  end interface

! ... # arguments
  iarg     = iargc()
!
! The program starts
!
  write (ou,'(3(/),30X,"PROGRAM FUKUI",X,A,2(/))') trim(version)
  call partial_time(ntime,'')
!
! Print help if # args is < 2
!
  if (iarg.LT.1) then
     write(ou,'(/,10X,"** WFN file is not included as argument **",/)')
     call help
     stop
  endif
!
! Does WFN file exist?
!
  call getarg (1,file1)
  fwfn(:)  = ''
  fwfn     = trim(file1(:index(file1,' ',.TRUE.)-1))
  inquire(file=fwfn,exist=lexist)
  if (.NOT.lexist) then
     fwfn  = trim(fwfn)//'.wfn'
     inquire(file=fwfn,exist=lexist)
     if (.NOT.lexist) stop '.wfn file NOT FOUND.'
  else
     file1 = trim(file1(:index(file1,'.',.TRUE.)-1))
  endif !! (.NOT.lexist) then
!
! grid file
!
  if (iarg.LT.2) then
     write(ou,'(/,10X,"** STK file is not included as argument **",/)')
     call help
     stop
  endif
  call getarg (2,file2)
  fstk(:)  = ''
  ext(:)   = ''
  fstk     = trim(file2(:index(file2,' ',.TRUE.)-1))
  inquire(file=fstk,exist=lfile2) !! I will check for final existence later
  if (.NOT.lfile2) then           !! Extension for stk file is optional
     igrid = 0 
     fstk  = trim(fstk)//'.stk'
     ext   = '.stk'
     inquire(file=fstk,exist=lfile2)
     if (.NOT.lfile2) stop '.stk/.cube file NOT FOUND.'
  else
     file2  = trim(file2(:index(file2,'.',.TRUE.)-1))
     ext(:) = ''
     ext    = trim(fstk(index(fstk,'.',.TRUE.)+1:))
     select case(trim(ext))
       case('stk')
          igrid = 0
       case('cub')
          igrid = 1
       case('cube')
          igrid = 1
       case default
          stop 'ERROR: Unknown extension for grid file'
     end select !! case(trim(ext))
  endif !! (.NOT.lfile2) then
!
! Output file
!
  ffuk = trim(file1(:index(file1,' ',.TRUE.)-1))//'_'//trim(file2(:index(file2,' ',.TRUE.)-1))//'.fuk'
  open (unit=ifuk,file=ffuk,status='unknown',form='formatted')
  write (ifuk,'(3(/),30X,"PROGRAM FUKUI",X,A,2(/))') trim(version)
  call print_main("Authors: Nicolas Otero, Marcos Mandado, Ricardo Mosquera")
  call print_main("License: GPLv3")
  call print_main("Wave function file found: "//trim(fwfn))
  call print_main("Grid file using ."//trim(gridext(igrid))//" format: "//trim(fstk))

  call cpu_time(cptim0)
!
! Check parallelization
!
  nproc = 0
  fproc = file1(1:index(file1,' ')-1)//'.proc'
  inquire(file=fproc,exist=lproc)
  if (lproc) then
     open (unit=iproc,file=fproc,status='old',form='formatted')
     read (iproc,*,iostat=i) nproc
     if (is_iostat_end(i)) then
        nproc = 0
     else
        if (i.NE.0) stop 'ERROR: Some strange character(s) in proc file. Integer value only allowed'
     endif !! (is_iostat_end(i)) then
     close(iproc)
     call print_main(trim(fproc)//" was found. Parallelization enabled!!")
     call init_par
     call print_rdarr("# threads to be used")
     write(ou  ,'(A)') int2str(nproc)
     write(ifuk,'(A)') int2str(nproc)
  else
     runparal = .FALSE.
  endif !! (lproc) then
!
! Check if the default max array dimensions will be modified
!
  fmxdm = file1(1:index(file1,' ')-1)//'.mxal'
  inquire(file=fmxdm,exist=lexist)
  if (lexist) then
     open (unit=imxdm,file=fmxdm,status='unknown',form='formatted')
     call print_main(trim(fmxdm)//" was found. Reading new max array allocations")
     if (nlines(imxdm).NE.4) stop 'Check .mxal file. # lines != 4!!'
     read (imxdm,*) natomx,nommx,npmx,nprimx
     close(imxdm)
     write(ou  ,'(10X,5(A10))') 'natomx','nommx','npmx','nprimx'
     write(ou  ,'(10X,5(I10))')  natomx , nommx , npmx , nprimx
     write(ifuk,'(10X,5(A10))') 'natomx','nommx','npmx','nprimx'
     write(ifuk,'(10X,5(I10))')  natomx , nommx , npmx , nprimx
  else
     call print_main("Using default max array allocations:")
     write(ou  ,'(10X,5(A10))') 'natomx','nommx','npmx','nprimx'
     write(ou  ,'(10X,5(I10))')  natomx , nommx , npmx , nprimx
     write(ifuk,'(10X,5(A10))') 'natomx','nommx','npmx','nprimx'
     write(ifuk,'(10X,5(I10))')  natomx , nommx , npmx , nprimx
  endif !! (lexist) then
!
! Activate sigma/pi separation if a .SG or .SOM file exists. Detect .PI file too.
!
  sigpi   = .FALSE.
  lsg     = .FALSE.
  fsom(:) = ''
  fsom    = file1(1:index(file1,' ')-1)//'.sg'
  inquire(file=fsom,exist=sigpi)
  if (sigpi) then ! sg  files
     open(unit=isom,file=fsom,status='old',form='formatted')
     call print_main(".sg file detected. Sigma/Pi separation will be performed by using "//trim(fsom))
     lsg  = .TRUE.
  else            ! som files
     fsom(:) = ''
     fsom    = file1(1:index(file1,' ')-1)//'.som'
     inquire(file=fsom,exist=sigpi)
     if (sigpi) call print_main(".som file detected. Sigma/Pi separation will be performed by using "//trim(fsom))
  endif !! (sigpi) then
  if (sigpi) then
     fpi=file1(1:index(file1,' ')-1)//'.pi'
     inquire(file=fpi,exist=lpi)
     if (lpi) call print_main(".pi file detected ("//trim(fpi)//"). Only pi orbitals will be employed")
  endif !! (lpi) then
!
! Activate primitive skipping
!
  if (iarg.GE.3) then
     call getarg(3,file2)
     read (file2,*,iostat=i) iato
     if (i.NE.0) stop 'ERROR: optional atom label is uncorrectly specified'
     fskpd = file1(1:index(file1,' ')-1)//'.skd'
     inquire(file=fskpd,exist=lskpd)
     if (lskpd) then
        open (unit=iskp,file=fskpd,status='old',form='formatted')
        read (iskp,*,iostat=i) ctffd
        if (i.NE.0) ctffd = ctffddef
        close(iskp)
     else
        ctffd = ctffddef
     endif !! (lskpd) then
     call print_main("Atom "//int2str(iato)//" was specified as 3rd argument. Skipping gaussian functions.")
     call print_rdarr("Distance in AU from this nucleus")
     call dbl2str(ctffd,3,ifuk)
  else
     iato = 0
  endif !! (iarg.GE.3) then
 
! Activate points skipping if .SKG and/or .SKR files are detected. Read cut-offs if needed
!
  fskp  = file1(1:index(file1,' ')-1)//'.skg'
  inquire(file=fskp ,exist=lskp )
  fskpr = file1(1:index(file1,' ')-1)//'.skr'
  inquire(file=fskpr,exist=lskpr)
  if (lskp ) then
     open (unit=iskp,file=fskp,status='old',form='formatted')
     read (iskp,*,iostat=i) ctffg
     if (i.NE.0) ctffg = ctffgdef
     close(iskp)
     if (igrid.EQ.0) write(ou ,'(" >> .skg file detected (",A,"). Skipping points with low weight(" &
        &                     ,1PE10.3,"=cut-off)")') trim(fskp),ctffg
     if (lfile2)    write(ifuk,'(" >> .skg file detected (",A,"). Skipping points with low weight(" &
        &                     ,1PE10.3,"=cut-off)")') trim(fskp),ctffg
  endif !! (lskp ) then
  if (lskpr) then
     open (unit=iskp,file=fskpr,status='old',form='formatted')
     read (iskp,*,iostat=i) ctffr
     if (i.NE.0) ctffr = ctffrdef
     close(iskp)
     write(ou  ,'(" >> .skr file detected (",A,"). Skipping points with negligible ref. dens.("     &
        &        ,1PE10.3,"=cut-off)")') trim(fskpr),ctffr
     if (lfile2) write(ifuk,'(" >> .skr file detected (",A,"). Skipping points with negligible ref. dens.(", &
        &                   1PE10.3,"=cut-off)")') trim(fskpr),ctffr
  endif !! (lskpr) then
!
! Activate Atomic Overlap Matrix calculation if .AOM file is detected.
! 
  faom  = file1(1:index(file1,' ')-1)//'.aom'
  inquire(file=faom ,exist=laom )
  if (laom ) call print_main(".aom file detected ("//trim(faom)//"). Atomic overlap matrix will be calculated")
!
! Opening files
!
  open(unit=iwfn,file=trim(fwfn),status='old',form='formatted')
!
! Does grid file exist?
!
  call print_main("Opening grid file: "//trim(fstk))
  inquire(file=fstk,exist=lexist)
  if (.NOT.lexist) stop 'ERROR: grid file NOT FOUND.'
  call openstk(istk,trim(fstk),lfmtstk,i)
  if (lfmtstk) then
     call print_second("formatted grid file detected after "//int2str(i)//" iterations")
  else
     call print_second("unformatted grid file detected after "//int2str(i)//" iterations")
  endif !! (lfmtstk) then
!
! File with sigma/pi MO separation
!
  if (sigpi) open(unit=isom,file=trim(fsom),status='old',form='formatted')

  call partial_time(ntime,'Settings prepared in')
!
! Using main computing subroutine
!
  call       goon(iato,ntime)
!
! Close files:
!
  inquire(unit=iwfn,opened=ltmp)
  if (ltmp) close(unit=iwfn)
  inquire(unit=isom,opened=ltmp)
  if (ltmp) close(unit=isom)
  inquire(unit=istk,opened=ltmp)
  close(unit=istk)

!
! Print total and partial times, parallelization performance
!
  call cpu_time(totim0)
  call partial_time(ntime,'',.TRUE.)
  close(ifuk)

!!!! Ver lo de cerrar ficheros, comprobar allocates


end program
