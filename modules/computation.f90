! ** modules/computation.f90 >> Main computation module of fukui program
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

module computation

  use :: main,                             only: i4,dp,ou,int2str
  use :: omp_lib

  implicit none

  integer(kind = i4),allocatable,dimension(:) :: icma,icmi,inddist
  integer(kind = i4)                          :: icalcat
  real   (kind = dp),allocatable,dimension(:) :: poom,om,s
  real   (kind = dp)                          :: popul,rhon

  contains

  subroutine orde_ic
  !
    use :: wfn        ,    only: nato,primitives
  
    implicit none
  
    integer(kind = i4) :: icc,i
  
    allocate(icma(nato),icmi(nato))                                   !! I1
    icmi(1) = 1
!$omp parallel default(none) &
!$omp&         shared(icma,primitives )
!$omp   workshare
    icma(1) = count(primitives%centers.EQ.1)
!$omp   end workshare
!$omp end parallel
    do i = 2,nato
!$omp parallel default(none) &
!$omp&         shared(primitives,icc,i )
!$omp   workshare
       icc     = count(primitives%centers.EQ.i)
!$omp   end workshare
!$omp end parallel
       icmi(i) = icma(i-1) + 1
       icma(i) = icma(i-1) + icc
    enddo !! i = 1,nato
  
  ! do i=1,nato                                                                                     !! DEBUG
  !    write (*,'(10X,"Gaussian functions ",I4," to ",I4," centred on atom",I3)') icmi(i),icma(i),i !! DEBUG
  ! enddo !! i=1,nato                                                                               !! DEBUG
  
  end subroutine orde_ic

  subroutine skip_gaussians(iato)

    use :: main           ,   only: ctffd,print_rdarr,dbl2str,ifuk
    use :: wfn            ,   only: nato,nuclei,nprim

    implicit none

    real   (kind = dp),allocatable,dimension(:) :: dist
    real   (kind = dp)                          :: vec(3)
    integer(kind = i4),intent(in)               :: iato
    integer(kind = i4)                          :: i,iskprim

    allocate(   dist(nato))                                           !! R1
    allocate(inddist(nato))                                           !! I1
    call print_rdarr("Percentage of gaussian functions used")
    if (iato.GT.0) then
!$omp parallel default(none)                  &
!$omp&          shared(nuclei,nato,iato,dist) &
!$omp&         private(vec,i)
!$omp do schedule(dynamic)
       do i = 1,nato
          vec(1)  = nuclei(i)%x-nuclei(iato)%x
          vec(2)  = nuclei(i)%y-nuclei(iato)%y
          vec(3)  = nuclei(i)%z-nuclei(iato)%z
          dist(i) = sqrt(dot_product(vec,vec))
       enddo !! i = 1,nato
!$omp end do
!$omp end parallel
       icalcat = 0
       do i = 1,nato
          if (dist(i).LE.ctffd) then
             icalcat          = icalcat + 1
             inddist(icalcat) = i
          endif !! (dist(i).LE.ctffd) then
       enddo !! i = 1,nato
!$omp parallel default(none)                    &
!$omp&          shared(dist,inddist,ctffd,icalcat,icma,icmi,iskprim)
!$omp workshare
       iskprim = sum(icma(inddist(:icalcat))-icmi(inddist(:icalcat))+1)
!$omp end workshare
!$omp end parallel
       call dbl2str(100*dble(iskprim)/dble(nprim),2,iprint=ifuk,post='%')
    else
       icalcat = nato
       inddist = (/ (i, i = 1,nato) /)
       call dbl2str(100._dp,2,iprint=ifuk,post='%')
    endif !! (iato.GT.0) then
    deallocate(dist)

  end subroutine skip_gaussians

  subroutine maincalc(ipart)

    use :: wfn ,                              only: nom,ialfa,nprim,primitives,nuclei,uhf,c,pel,matom,nbac &
                                                  , deallocate_wfn
    use :: grid,                              only: nptos,points

    implicit none

    real   (kind = dp)                          :: rrr,vec(3)
    integer(kind = i4)                          :: ipart
    integer(kind = i4)                          :: i,j,k,l,iom,jom,ifrag,id
    real   (kind = dp),allocatable,dimension(:) :: ggg

    allocate(s(nbac*(nbac+1)/2),poom(nom),ggg(nprim),om(nom))
    popul = 0._dp
!$omp parallel default(none) &
!$omp&          shared(poom,s)
!$omp workshare
    poom  = 0._dp
    s     = 0._dp
!$omp end workshare
!$omp end parallel

    if (uhf) then !!!!!!!!! UHF !!!!!!!!!!!!!!

!$omp parallel default(none) &
!$omp&         private(ifrag,i,j,vec,rrr,k,ggg,om,rhon,l,iom,jom, id ) &
!$omp&          shared(ipart,nptos,icalcat,points,nuclei,inddist,icmi  &
!$omp&                ,icma,nprim,primitives,nom,c,matom,pel,ialfa   ) &
!$omp&         reduction(+:popul,poom,s)
!$omp do schedule(dynamic)
       do ifrag = 1,ipart
          do i=(ifrag-1)*nptos/ipart+1,ifrag*nptos/ipart
             do j=1,icalcat
                vec(1) = points(i)%x-nuclei(inddist(j))%x
                vec(2) = points(i)%y-nuclei(inddist(j))%y
                vec(3) = points(i)%z-nuclei(inddist(j))%z
                rrr    = dot_product(vec,vec)
                do k=icmi(inddist(j)),icma(inddist(j))
                   call gaussi(vec(1),vec(2),vec(3),rrr,k,ggg(k))
                enddo !! k=icmi(inddist(j)),icma(inddist(j))
             enddo !! j=1,icalcat
             om   = 0._dp
             rhon = 0._dp
             do j=1,nom
                do l=1,icalcat
                   do k=icmi(inddist(l)),icma(inddist(l))
                      om(j)=om(j)+c(matom(j),k)*ggg(k)
                   enddo !! k=icmi(inddist(l)),icma(inddist(l))
                enddo !! l=1,icalcat
                rhon    = rhon   + om(j)*om(j)*pel(matom(j))
                poom(j) = poom(j)+ om(j)*om(j)*points(i)%w
             enddo !! j=1,nom
             popul = popul+rhon*points(i)%w

! ... Unrestricted Atomic Overlap Matrix
             do iom=1,ialfa
                do jom=1,iom
                   k = matom(jom) + matom(iom)*(matom(iom)-1)/2
! omp atomic update
                   s(k) = s(k) + om(iom) * om(jom) * points(i)%w
                enddo !! jom=1,iom
             enddo !! iom=1,ialfa
             do iom=ialfa+1,nom
                do jom=ialfa+1,iom
                   k = matom(jom) + matom(iom)*(matom(iom)-1)/2
! omp atomic update
                   s(k) = s(k) + om(iom) * om(jom) * points(i)%w
                enddo !! jom=ialfa+1,iom
             enddo !! iom=ialfa+1,nom

          enddo !! i=(ifrag-1)*nptos/10+1,ifrag*nptos/10
          call prwkdone(ifrag,nptos,ipart)
       enddo !! ifrag = 1,ipart
!$omp enddo
!$omp end parallel

    else    !!!!!!!!!!!!!!! RHF !!!!!!!!!!!!!!

!$omp parallel default(none) &
!$omp&         private(ifrag,i,j,vec,rrr,k,ggg,om,rhon,l,iom,jom, id ) &
!$omp&          shared(ipart,nptos,icalcat,points,nuclei,inddist,icmi  &
!$omp&                ,icma,nprim,primitives,nom,c,matom,pel         ) &
!$omp&         reduction(+:popul,poom,s)
!$omp do schedule(dynamic)
       do ifrag = 1,ipart
          do i=(ifrag-1)*nptos/ipart+1,ifrag*nptos/ipart
             do j=1,icalcat
                vec(1) = points(i)%x-nuclei(inddist(j))%x
                vec(2) = points(i)%y-nuclei(inddist(j))%y
                vec(3) = points(i)%z-nuclei(inddist(j))%z
                rrr    = dot_product(vec,vec)
                do k=icmi(inddist(j)),icma(inddist(j))
                   call gaussi(vec(1),vec(2),vec(3),rrr,k,ggg(k))
                enddo !! k=icmi(inddist(j)),icma(inddist(j))
             enddo !! j=1,icalcat
             om   = 0._dp
             rhon = 0._dp
             do j=1,nom
                do l=1,icalcat
                   do k=icmi(inddist(l)),icma(inddist(l))
                      om(j) = om(j) + c(matom(j),k)*ggg(k)
                   enddo !! k=icmi(inddist(l)),icma(inddist(l))
                enddo !! l=1,icalcat
                rhon    = rhon   + om(j)*om(j)*pel(matom(j))
                poom(j) = poom(j)+ om(j)*om(j)*points(i)%w
             enddo !! j=1,nom
             popul = popul+rhon*points(i)%w

! ... Restricted Atomic Overlap Matrix
             do iom=1,nom
                do jom=1,iom
                   k = matom(jom) + matom(iom)*(matom(iom)-1)/2
                   s(k) = s(k) + om(iom) * om(jom) * points(i)%w
                enddo !! jom=1,iom
             enddo !! iom=1,nom

          enddo !! i=(ifrag-1)*nptos/10+1,ifrag*nptos/10
          call prwkdone(ifrag,nptos,ipart)
       enddo !! ifrag = 1,ipart

!$omp enddo
!$omp end parallel

      endif !! (uhf) then
      deallocate(ggg,icma,icmi,inddist,points)
      call deallocate_wfn

  end subroutine maincalc


  subroutine maindens(ipart)

    use :: wfn ,                              only: nom,ialfa,nprim,primitives,nuclei,uhf,c,pel,matom &
                                                  , deallocate_wfn
    use :: grid,                              only: nptos,points

    implicit none

    real   (kind = dp)                          :: rrr,vec(3)
    integer(kind = i4)                          :: ipart
    integer(kind = i4)                          :: i,j,k,l,iom,jom,ifrag,id
    real   (kind = dp),allocatable,dimension(:) :: ggg

    allocate(poom(nom),ggg(nprim),om(nom))
    popul = 0._dp
!$omp parallel default(none) &
!$omp&          shared(poom)
!$omp workshare
    poom  = 0._dp
!$omp end workshare
!$omp end parallel

    if (uhf) then !!!!!!!!! UHF !!!!!!!!!!!!!!

!$omp parallel default(none) &
!$omp&         private(ifrag,i,j,vec,rrr,k,ggg,om,rhon,l,iom,jom, id ) &
!$omp&          shared(ipart,nptos,icalcat,points,nuclei,inddist,icmi  &
!$omp&                ,icma,nprim,primitives,nom,c,matom,pel,ialfa   ) &
!$omp&         reduction(+:popul,poom)
!$omp do schedule(dynamic)
       do ifrag = 1,ipart
          call prwkdone(ifrag,nptos,ipart)
          do i=(ifrag-1)*nptos/ipart+1,ifrag*nptos/ipart
             do j=1,icalcat
                vec(1) = points(i)%x-nuclei(inddist(j))%x
                vec(2) = points(i)%y-nuclei(inddist(j))%y
                vec(3) = points(i)%z-nuclei(inddist(j))%z
                rrr    = dot_product(vec,vec)
                do k=icmi(inddist(j)),icma(inddist(j))
                   call gaussi(vec(1),vec(2),vec(3),rrr,k,ggg(k))
                enddo !! k=icmi(inddist(j)),icma(inddist(j))
             enddo !! j=1,icalcat
             om   = 0._dp
             rhon = 0._dp
             do j=1,nom
                do l=1,icalcat
                   do k=icmi(inddist(l)),icma(inddist(l))
                      om(j)=om(j)+c(matom(j),k)*ggg(k)
                   enddo !! k=icmi(inddist(l)),icma(inddist(l))
                enddo !! l=1,icalcat
                rhon    = rhon   + om(j)*om(j)*pel(matom(j))
                poom(j) = poom(j)+ om(j)*om(j)*points(i)%w
             enddo !! j=1,nom
             popul = popul+rhon*points(i)%w

          enddo !! i=(ifrag-1)*nptos/10+1,ifrag*nptos/10
          call prwkdone(ifrag,nptos,ipart)
       enddo !! ifrag = 1,ipart
!$omp enddo
!$omp end parallel

    else    !!!!!!!!!!!!!!! RHF !!!!!!!!!!!!!!

!$omp parallel default(none) &
!$omp&         private(ifrag,i,j,vec,rrr,k,ggg,om,rhon,l,iom,jom, id ) &
!$omp&          shared(ipart,nptos,icalcat,points,nuclei,inddist,icmi  &
!$omp&                ,icma,nprim,primitives,nom,c,matom,pel         ) &
!$omp&         reduction(+:popul,poom)
!$omp do schedule(dynamic)
       do ifrag = 1,ipart
          do i=(ifrag-1)*nptos/ipart+1,ifrag*nptos/ipart
             do j=1,icalcat
                vec(1) = points(i)%x-nuclei(inddist(j))%x
                vec(2) = points(i)%y-nuclei(inddist(j))%y
                vec(3) = points(i)%z-nuclei(inddist(j))%z
                rrr    = dot_product(vec,vec)
                do k=icmi(inddist(j)),icma(inddist(j))
                   call gaussi(vec(1),vec(2),vec(3),rrr,k,ggg(k))
                enddo !! k=icmi(inddist(j)),icma(inddist(j))
             enddo !! j=1,icalcat
             om   = 0._dp
             rhon = 0._dp
             do j=1,nom
                do l=1,icalcat
                   do k=icmi(inddist(l)),icma(inddist(l))
                      om(j) = om(j) + c(matom(j),k)*ggg(k)
                   enddo !! k=icmi(inddist(l)),icma(inddist(l))
                enddo !! l=1,icalcat
                rhon    = rhon   + om(j)*om(j)*pel(matom(j))
                poom(j) = poom(j)+ om(j)*om(j)*points(i)%w
             enddo !! j=1,nom
             popul = popul+rhon*points(i)%w
          enddo !! i=(ifrag-1)*nptos/10+1,ifrag*nptos/10
          call prwkdone(ifrag,nptos,ipart)
       enddo !! ifrag = 1,ipart

!$omp enddo
!$omp end parallel

      endif !! (uhf) then
      deallocate(ggg,icma,icmi,inddist,points)
      call deallocate_wfn

  end subroutine maindens

  subroutine prwkdone(nfrag,ntot,ipart)

      use :: omp_lib        !! Enabling OpenMP functions !!
      use :: main    ,          only : print_second
      implicit none
      integer(kind = i4),intent(in) :: nfrag,ntot,ipart

      call print_second("CPU "                     &
                // int2str(omp_get_thread_num())   &
                // ": from "                       &
                // int2str((nfrag-1)*ntot/ipart+1) &
                // " to "                          &
                // int2str(nfrag*ntot/ipart)       &
                // " of "                          &
                // int2str(ntot) &
                       )

  end subroutine prwkdone

  subroutine gaussi(x,y,z,rrr,k,ggg)
!
!         calcula el valor de la gaussiana k en el punto x,y,z como ggg
!         utiliza rrr. x,y,z son coordenadas relativas al nucleo
!
      use :: wfn,                 only: primitives
      implicit none
      integer(kind = i4),intent(in ) :: k
      real   (kind = dp),intent(in ) :: x,y,z,rrr
      real   (kind = dp)             :: ggg

       if (rrr.LE.1.0d-100) then
          ggg=1._dp
       else
          ggg=dexp(-primitives(k)%expo*rrr)
       endif

       select case(primitives(k)%types)
         case( 1)
           return
         case( 2)
           ggg = ggg * x
           return
         case( 3)
           ggg = ggg * y
           return
         case( 4)
           ggg = ggg * z
         case( 5)
           ggg = ggg * x * x
         case( 6)
           ggg = ggg * y * y
         case( 7)
           ggg = ggg * z * z
         case( 8)
           ggg = ggg * x * y
         case( 9)
           ggg = ggg * x * z
         case(10)
           ggg = ggg * y * z
         case(11)
           ggg = ggg * x * x * x
         case(12)
           ggg = ggg * y * y * y
         case(13)
           ggg = ggg * z * z * z
         case(14)
           ggg = ggg * x * x * y
         case(15)
           ggg = ggg * x * x * z
         case(16)
           ggg = ggg * y * y * z
         case(17)
           ggg = ggg * x * y * y
         case(18)
           ggg = ggg * x * z * z
         case(19)
           ggg = ggg * y * z * z
         case(20)
           ggg = ggg * x * y * z
         case default
           write(ou,'(" ** This type assig. is unknown (",I3,") **")') primitives(k)%types
           stop
       end select

  end subroutine gaussi

end module computation
