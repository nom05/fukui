! ** modules/wfn.f90 >> Main module to read wfn file
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

module wfn

  use ,intrinsic :: iso_fortran_env,                only: output_unit,real64,int32
  use            :: main           ,                only: print_rdarr,ifuk
  use            :: omp_lib

  implicit none

  integer                            ,parameter,private:: dp = real64
  integer                            ,parameter,private:: i4 = int32
  integer                            ,parameter,private:: ou = output_unit

  integer(kind = i4)                           ,public :: nom,nprim,nato,nheavy,ialfa,nbac
  logical                                      ,public :: uhf
  integer(kind = i4),allocatable,dimension(:  ),public :: matom
  real   (kind = dp),allocatable,dimension(:  ),public :: pel,eom
  real   (kind = dp),allocatable,dimension(:,:),public :: c
  type,public :: xyzZ
    real   (kind = dp)                                 :: x,y,z
    integer(kind = i4)                                 :: ZZ
  end type       xyzZ
  type,public :: prim_info
    real   (kind = dp)                                 :: expo
    integer(kind = i4)                                 :: centers,types
  end type       prim_info
  type   (xyzZ     ),allocatable,dimension(:  ),public :: nuclei
  type   (prim_info),allocatable,dimension(:  ),public :: primitives

  interface     readblock
    module procedure rdarr_i4, rdarr_dp, rdarr_mo, rdarr_xyzq
  end interface readblock

  contains

    subroutine rdarr_i4(arr,inp,ndim,text)
      implicit none
      integer  (kind = i4)              :: i
      integer  (kind = i4),intent(in)   :: ndim,inp
      integer  (kind = i4),dimension(:) :: arr
      character(len  =  *),intent(in)   :: text
      call     print_rdarr(text)
      read (inp,'(20X,20I3)',iostat=i) arr(:ndim)
      if (i.NE.0) then
         write(ou,'("ERROR in ",A,"!!")') trim(text)
         stop
      endif !! (i.NE.0) then
!     write(ou ,'(20X,20I3)'         ) arr(:ndim)               !! DEBUG
      write(ou  ,'("OK")')
      write(ifuk,'("OK")')
    end subroutine rdarr_i4

    subroutine rdarr_dp(arr,inp,ndim,text)
      implicit none
      integer  (kind = i4)              :: i
      integer  (kind = i4),intent(in)   :: ndim,inp
      real     (kind = dp),dimension(:) :: arr
      character(len  =  *),intent(in)   :: text
      call     print_rdarr(text)
      read (inp,'(10X,5D14.7)',iostat=i) arr(:ndim)
      if (i.NE.0) then
         write(ou,'("ERROR in ",A,"!!")') trim(text)
         stop
      endif !! (i.NE.0) then
!     write(ou ,'(5D14.7)'         ) arr(:ndim)       !! DEBUG
      write(ou  ,'("OK")')
      write(ifuk,'("OK")')
    end subroutine rdarr_dp

    subroutine rdarr_mo(arr1,arr2,arr3,inp,ndim1,ndim2,text)
      implicit none
      integer  (kind = i4)                :: i,j
      integer  (kind = i4),intent(in)     :: ndim1,ndim2,inp
      real     (kind = dp),dimension(:,:) :: arr1
      real     (kind = dp),dimension(:  ) :: arr2,arr3
      character(len  =  *),intent(in)     :: text
      call     print_rdarr(text)
      do i=1,ndim1
         read (inp,'(37X,F10.7,15X,F12.6)',iostat=j,err=52) arr2(i),arr3(i)
    !    write(ou ,'(37X,F10.7,15X,F12.6)'                ) arr2(i),arr3(i) !! DEBUG
         read (inp,'(5D16.8)',iostat=j,err=52) arr1(i,:ndim2)
    !    write(ou ,'(5D16.8)'                ) arr1(i,1),arr1(i,ndim2)      !! DEBUG
      enddo !! i=1,nom
 52   if (j.NE.0) then
         write(ou,'("ERROR in ",A,"!!")') trim(text)
         stop
      endif !! (i.NE.0) then
      write(ou  ,'("OK")')
      write(ifuk,'("OK")')
    end subroutine rdarr_mo

    subroutine rdarr_xyzq(arrx,arry,arrz,arrq,inp,ndim,text)
      implicit none
      integer  (kind = i4)                :: i,j
      integer  (kind = i4),intent(in)     :: ndim,inp
      real     (kind = dp),dimension(:  ) :: arrx,arry,arrz
      integer  (kind = i4),dimension(:  ) :: arrq
      character(len  =  *),intent(in)     :: text
      call     print_rdarr(text)
      do i=1,ndim
         read (inp,'(24X,   3(F12.8),10X,I3)',err=53,iostat=j) arrx(i),arry(i),arrz(i),arrq(i)
!        write(ou ,'(24X,   3(F12.8),10X,I3)'                ) arrx(i),arry(i),arrz(i),arrq(i) !! DEBUG
      enddo !! i=1,ndim
 53   if (j.NE.0) then
         write(ou,'("ERROR in ",A,"!!")') trim(text)
         stop
      endif !! (i.NE.0) then
      write(ou  ,'("OK")')
      write(ifuk,'("OK")')
    end subroutine rdarr_xyzq

    subroutine readwfn
!
!  LE TODA A INFORMACION DO ARQUIVO .WFN INCLUINDO OS VALORES DA ENERXIA
!  ELECTRONICA TOTAL (enerx) E DA RELACION VIRIAL (virial) XUNTO COA
!  INFORMACION XEOMETRICA (nuclei%x,nuclei%y,nuclei%z,nuclei%ZZ), DE PRIMITIVAS (ic,it,ex),
!  OCUPACION E ENERXIA DOS OM (pel, eom) e MATRIZ DE COEFICIENTES DOS OM (c)
!  
!  R.A.Mosquera Castro,    Vigo  , 22 de febreiro de 1998
!  Nicolás Otero Martínez, Pau   ,  7 de decembro do 2015 ** ACTUALIZACIÓNS e BUGFIXES
!  Nicolás Otero Martínez, Pau   , 30 de    xunho do 2016 ** ADAPTACION DE CODIGO e
!                                                          DIMENSION PREESTABLECIDA DE
!  Nicolás Otero Martínez, Cangas,  9 de  xaneiro do 2022 ** MODULAR F90 VERSION
!
      use :: main           ,  only: i4,ou,iwfn,natomx,nommx,nprimx,int2str,print_main,print_second

      implicit none

      character                   :: lin1*80
      integer(kind = i4)          :: nunpair,i,k
      logical                     :: lcheck
! real   (kind = dp)          :: enerx,virial

      uhf   = .FALSE.
      ialfa = 0
      rewind(iwfn)
!
! Localiza unha palabra clave para ler informacion (compatibilidade con GAMESS)
!
      do 
         read (iwfn,'(A80)',iostat=i) lin1
!         write(*   ,'(A80)') lin1                                     !! DEBUG
         if ( index(trim(lin1),'GAUSSIAN').GT.0) exit
         if ( index(trim(lin1),'GTO'     ).GT.0) exit
         if (i.NE.0) stop 'ERROR: Wrong wfn file? Unable to find line with dimensions'
     enddo 
     backspace (iwfn)
!
! Read # of MOs, primitives and atoms
!
     read (iwfn,'(19X,I4,15X,I5,17X,I3)',iostat=i) nom,nprim,nato
     call print_second(int2str(nom)//" MOs, "//int2str(nprim)//" gaussian functions, "//int2str(nato)//" atoms.")
!
! Check dimensions.
!
     lcheck = .TRUE.
     if (natomx.LT.nato) then
        lcheck = .FALSE.
        write(ou,'(" ** Increase  natomx **",A,"=old value < nato =",A)') int2str(natomx),int2str(nato)
     endif !! (natomx.LT.nato) then
     if (nommx .LT.nom) then
        lcheck = .FALSE.
        write(ou,'(" ** Increase   nommx **",A,"=old value < nom  =",A)') int2str(nommx),int2str(nom)
     endif !! (nommx .LT.nom) then
     if (nprimx.LT.nprim) then
        lcheck = .FALSE.
        write(ou,'(" ** Increase  nprimx **",A,"=old value < nprim=",A)') int2str(nprimx),int2str(nprim)
     endif !! (nprimx.GE.nprim) then
     if (.NOT.lcheck) stop ' ** Edit source code in main module or create a file with extension .mxal following the manual **'
!
! Allocate main arrays
!
     allocate(nuclei(nato),primitives(nprim))
     allocate(pel(nom),eom(nom),c(nom,nprim))
!
! Reading data blocks
!
     call readblock(nuclei%x,nuclei%y,nuclei%z,nuclei%ZZ,iwfn,nato,"Coordinates and nuclear charges")

!$omp parallel default (none) shared ( nheavy,nuclei )
  !$omp workshare
     nheavy = count(nuclei%ZZ>1)
  !$omp end workshare
!$omp end parallel

     call readblock(primitives%centers,iwfn,nprim,"Primitive centers"  )

     call readblock(primitives%types  ,iwfn,nprim,"Primitive types"    )

     call readblock(primitives%expo   ,iwfn,nprim,"Primitive exponents")
!
!  Reading occupation numbers, MO energies and MO coefficients matrix
!
     call readblock(c,pel,eom,iwfn,nom,nprim,"MO coefficients")
!
!  130404 For G09 versions > A02. Detecting UHF. Improvable code
!
     if (nint(pel(1)).NE.2) then
        uhf   = .TRUE.
        call print_second("Unrestricted calculation detected.")
        ialfa = int(nom/2)+mod(nom,2)
        k     = ialfa
        do while (nint(eom(1)).NE.nint(eom(k)))  !! It sucks
           k = k+1
        enddo !! while (nint(eom(1).NE.nint(eom(k))
        ialfa = k-1
        nunpair = 2*ialfa-nom
        call print_second(                                                  &
             "# alpha e ="//int2str(ialfa)//" "//                           &
             int2str(nom-ialfa)//"= # beta e (mult="//int2str(2*ialfa-nom+1)&
                    //")")
        write (ou  , &
          '(10X,"ORB. E(last_alpha)=",F11.6,5X,"ORB. E(1st_beta)=",F11.6)') &
                                                    eom(ialfa),eom(ialfa+1)
        write (ifuk, &
          '(10X,"ORB. E(last_alpha)=",F11.6,5X,"ORB. E(1st_beta)=",F11.6)') &
                                                    eom(ialfa),eom(ialfa+1)
     endif !! (nint(pel(1)).NE.2) then
!
!  Read electronic energy and virial relationship
!
!    read (iwfn,*)  !! Not needed

   end subroutine readwfn

   subroutine deallocate_wfn
     implicit none
     if (allocated(nuclei))     deallocate(nuclei)
     if (allocated(primitives)) deallocate(primitives)
     if (allocated(c))          deallocate(c)
     if (allocated(eom))        deallocate(eom)
   end subroutine deallocate_wfn
   
   subroutine read_nuclei_coord_from_wfn
     use :: main,            only: i4,ou,iwfn,natomx,print_second,int2str
     implicit none
     integer  (kind = i4)       :: i
     character( len = 80)       :: lin1
     LOGICAL                    :: lcheck
     rewind(iwfn)
     do 
        read (iwfn,'(A80)',iostat=i) lin1
        if ( index(trim(lin1),'GAUSSIAN').GT.0) exit
        if ( index(trim(lin1),'GTO'     ).GT.0) exit
        if (i.NE.0) stop 'ERROR: Wrong wfn file? Unable to find line with dimensions'
     enddo 
     backspace(iwfn)
!
! Read # of MOs, primitives and atoms
!
     read (iwfn,'(60X,I3)',iostat=i) nato
     call print_second(int2str(nato)//" atoms.")
!
! Check dimensions.
!
     lcheck = .TRUE.
     if (natomx.LT.nato) then
        lcheck = .FALSE.
        write(ou,'(" ** Increase  natomx **",A,"=old value < nato =",A)') int2str(natomx),int2str(nato)
     endif !! (natomx.LT.nato) then
     if (.NOT.lcheck) stop ' ** Edit source code in main module or create a file with extension .mxal following the manual **'
!
! Allocate main arrays
!
     allocate(nuclei(nato))
!
! Reading data blocks
!
     call readblock(nuclei%x,nuclei%y,nuclei%z,nuclei%ZZ,iwfn,nato,"Coordinates and nuclear charges")

!$omp parallel default (none) shared ( nheavy,nuclei )
  !$omp workshare
     nheavy = count(nuclei%ZZ>1)
  !$omp end workshare
!$omp end parallel

   end subroutine read_nuclei_coord_from_wfn

end module wfn
