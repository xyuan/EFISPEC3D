!=====================================================================================================================================
!                                      EFISPEC3D                                             !
!                            (Elements FInis SPECtraux 3D)                                   !
!                                                                                            !
!                               http://efispec.free.fr                                       !
!                                                                                            !
!                                                                                            !
!                            This file is part of EFISPEC3D                                  !
!              Please refer to http://efispec.free.fr if you use it or part of it            !
!                                                                                            !
!                                                                                            !
!1 ---> French License: CeCILL V2                                                            !
!                                                                                            !
!         Copyright BRGM 2009 "contributeurs : Florent  DE MARTIN (BRGM)                     !
!                                              David    MICHEA    (BRGM)                     !
!                                              Philippe THIERRY   (Intel)"                   !
!                                                                                            !
!         Contact: f.demartin at brgm.fr                                                     !
!                                                                                            !
!         Ce logiciel est un programme informatique servant a resoudre l'equation du         !
!         mouvement en trois dimensions via une methode des elements finis spectraux.        !
!                                                                                            !
!         Ce logiciel est regi par la licence CeCILL soumise au droit francais et            !
!         respectant les principes de diffusion des logiciels libres. Vous pouvez            !
!         utiliser, modifier et/ou redistribuer ce programme sous les conditions de la       !
!         licence CeCILL telle que diffusee par le CEA, le CNRS et l'INRIA sur le site       !
!         "http://www.cecill.info".                                                          !
!                                                                                            !
!         En contrepartie de l'accessibilite au code source et des droits de copie, de       !
!         modification et de redistribution accordes par cette licence, il n'est offert      !
!         aux utilisateurs qu'une garantie limitee. Pour les memes raisons, seule une        !
!         responsabilite restreinte pese sur l'auteur du programme, le titulaire des         !
!         droits patrimoniaux et les concedants successifs.                                  !
!                                                                                            !
!         A cet egard l'attention de l'utilisateur est attiree sur les risques associes      !
!         au chargement, a l'utilisation, a la modification et/ou au developpement et a      !
!         la reproduction du logiciel par l'utilisateur etant donne sa specificite de        !
!         logiciel libre, qui peut le rendre complexe a manipuler et qui le reserve donc     !
!         a des developpeurs et des professionnels avertis possedant des connaissances       !
!         informatiques approfondies. Les utilisateurs sont donc invites a charger et        !
!         tester l'adequation du logiciel a leurs besoins dans des conditions permettant     !
!         d'assurer la securite de leurs systemes et ou de leurs donnees et, plus            !
!         generalement, a l'utiliser et l'exploiter dans les memes conditions de             !
!         securite.                                                                          !
!                                                                                            !
!         Le fait que vous puissiez acceder a cet en-tete signifie que vous avez pris        !
!         connaissance de la licence CeCILL et que vous en avez accepte les termes.          !
!                                                                                            !
!                                                                                            !
!2 ---> International license: GNU GPL V3                                                    !
!                                                                                            !
!         EFISPEC3D is a computer program that solves the three-dimensional equations of     !
!         motion using a finite spectral-element method.                                     !
!                                                                                            !
!         Copyright (C) 2009 Florent DE MARTIN                                               !
!                                                                                            !
!         Contact: f.demartin at brgm.fr                                                     !
!                                                                                            !
!         This program is free software: you can redistribute it and/or modify it under      !
!         the terms of the GNU General Public License as published by the Free Software      !
!         Foundation, either version 3 of the License, or (at your option) any later         !
!         version.                                                                           !
!                                                                                            !
!         This program is distributed in the hope that it will be useful, but WITHOUT ANY    !
!         WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    !
!         PARTICULAR PURPOSE. See the GNU General Public License for more details.           !
!                                                                                            !
!         You should have received a copy of the GNU General Public License along with       !
!         this program. If not, see http://www.gnu.org/licenses/.                            !
!                                                                                            !
!                                                                                            !
!3 ---> Third party libraries                                                                !
!                                                                                            !
!         EFISPEC3D uses the following of third party libraries:                             !
!                                                                                            !
!           --> METIS 5.1.0 under the Apache License Version 2.0 (compatible with GNU GPL V3)!
!               see http://glaros.dtc.umn.edu/gkhome/metis/metis/overview                    !
!                                                                                            !
!           --> Lib_VTK_IO under GNU GPL V3 License                                          !
!               see S. Zaghi's website: https://github.com/szaghi/Lib_VTK_IO                 !
!                                                                                            !
!           --> INTERP_LINEAR under GNU GPL License                                          !
!               see J. Burkardt website: http://people.sc.fsu.edu/~jburkardt/                !
!                                                                                            !
!                                                                                            !
!4 ---> Related Articles                                                                     !
!                                                                                            !
!         De Martin, F., Matsushima, M., Kawase, H. (BSSA, 2013)                             !
!            Impact of geometric effects on near-surface Green's functions                   !
!            doi:10.1785/0120130039                                                          !
!                                                                                            !
!         Aochi, H., Ducellier, A., Dupros, F., Delatre, M., Ulrich, T., De Martin, F.,      !
!         Yoshimi, M., (Pure Appl. Geophys. 2013)                                            !
!            Finite Difference Simulations of Seismic Wave Propagation for the 2007 Mw 6.6   !
!            Niigata-ken Chuetsu-Oki Earthquake: Validity of Models and Reliable             !
!            Input Ground Motion in the Near-Field                                           !
!            doi:10.1007/s00024-011-0429-5                                                   !
!                                                                                            !
!         De Martin, F. (BSSA, 2011)                                                         !
!            Verification of a Spectral-Element Method Code for the Southern California      !
!            Earthquake Center LOH.3 Viscoelastic Case                                       !
!            doi:10.1785/0120100305                                                          !
!                                                                                            !
!=====================================================================================================================================

!>@file
!!This files contains a module to initialize MPI buffers between cpu myrank and its neighbor cpus.

!>@brief
!!This module contains subroutines to initialize MPI buffers between cpu myrank and its neighbor cpus.
module mod_init_mpibuffer

   implicit none

   private

   public  :: init_mpi_buffers
   private :: create_gll_buffer_recv
   private :: remove_duplicate
   
   contains

!
!
!>@brief
!!This subroutine searches for common GLLs between cpu myrank and neighbor cpus.
!!Send and receive buffers are filled only with unique GLLs (i.e., no duplicated GLL are present in MPI buffers).
!***********************************************************************************************************************************************************************************
   subroutine init_mpi_buffers()
!***********************************************************************************************************************************************************************************

       use mpi

       use mod_global_variables, only :&
                                       IG_NGLL&
                                      ,ig_myrank&
                                      ,ig_ncpu_neighbor&
                                      ,ig_nhexa_outer&
                                      ,ig_hexa_gll_glonum&
                                      ,tg_cpu_neighbor&
                                      ,ig_cpu_neighbor_info&
                                      ,ig_nneighbor_all_kind&
                                      ,rg_gll_coordinate&
                                      ,ig_mpi_buffer_offset&
                                      ,rg_mpi_buffer_send&
                                      ,rg_mpi_buffer_recv&
                                      ,IG_NDOF&
                                      ,ig_mpi_request_send&
                                      ,ig_mpi_request_recv&
                                      ,error_stop&
                                      ,cg_prefix&
                                      ,cg_myrank&
                                      ,get_newunit&
                                      ,LG_ASYNC_MPI_COMM&
                                      ,LG_OUTPUT_DEBUG_FILE&
                                      ,IG_LST_UNIT&
                                      ,ig_mpi_buffer_sizemax

       implicit none

       integer, parameter                   :: NFACE =  6
       integer, parameter                   :: NEDGE = 12
       integer, parameter                   :: NNODE =  8
                
       real   , dimension(:,:), allocatable :: gll_coord_send
       real   , dimension(:,:), allocatable :: gll_coord_recv
                
       integer, dimension(MPI_STATUS_SIZE)  :: statut
       integer                              :: i
       integer                              :: j
       integer                              :: ios
       integer                              :: myunit_debug
       integer                              :: icpu_neighbor
       integer                              :: icpu
       integer                              :: isurf
       integer                              :: ngll_duplicate
       integer                              :: ngll_send
       integer                              :: ngll_recv
       integer                              :: num_gll
       integer                              :: surf_num
       integer                              :: elt_num
       integer                              :: ngll_unique
       integer                              :: mpi_buffer_sizemax
       integer                              :: buffer_sizemax
       integer, dimension(:), allocatable   :: buffer_gll_duplicate

       character(len=  6)                   :: cl_rank
       character(len=255)                   :: info

       if (ig_myrank == 0) then
          write(IG_LST_UNIT,'(" ",/,a)') "creating mpi buffers between cpu myrank and its neighbors..."
          call flush(IG_LST_UNIT)
       endif

       buffer_sizemax = (6*IG_NGLL*IG_NGLL + 12*IG_NGLL + 8) * ig_nhexa_outer !DAVID: worst case, may be reduced
       allocate(buffer_gll_duplicate(buffer_sizemax),stat=ios)
       if (ios /= 0) then
          write(info,'(a)') "error in subroutine init_mpi_buffers while allocating buffer_gll_duplicate"
          call error_stop(info)
       endif
        
       if (LG_ASYNC_MPI_COMM) then
           allocate(ig_mpi_request_send(ig_ncpu_neighbor))
           allocate(ig_mpi_request_recv(ig_ncpu_neighbor))
           allocate(ig_mpi_buffer_offset(ig_ncpu_neighbor+1))
           ig_mpi_buffer_offset(:) = 0
           ig_mpi_request_send(:)  = 0
           ig_mpi_request_recv(:)  = 0
       endif
       
       mpi_buffer_sizemax = 0

       do icpu_neighbor = 1,ig_ncpu_neighbor

           icpu                    = tg_cpu_neighbor(icpu_neighbor)%icpu                  
           buffer_gll_duplicate(:) = 0
           ngll_duplicate          = 0

           do isurf = 1,ig_nneighbor_all_kind

               if (ig_cpu_neighbor_info(3,isurf) == icpu) then

                 surf_num = ig_cpu_neighbor_info(2,isurf)
                 elt_num  = ig_cpu_neighbor_info(1,isurf)
!
!----------------->put gll points from contact face in the buffer 
                   if (surf_num <= NFACE) then
                       do i = 1,IG_NGLL
                           do j = 1,IG_NGLL
                               select case(surf_num)
                                   case(1)
                                       num_gll = ig_hexa_gll_glonum(j,i,1,elt_num)
                                   case(2)
                                       num_gll = ig_hexa_gll_glonum(j,1,i,elt_num)
                                   case(3)
                                       num_gll = ig_hexa_gll_glonum(IG_NGLL,j,i,elt_num)
                                   case(4)
                                       num_gll = ig_hexa_gll_glonum(j,IG_NGLL,i,elt_num)
                                   case(5)
                                       num_gll = ig_hexa_gll_glonum(1,j,i,elt_num)
                                   case(6)
                                       num_gll = ig_hexa_gll_glonum(j,i,IG_NGLL,elt_num)
                               end select

                               ngll_duplicate = ngll_duplicate + 1

                               if (ngll_duplicate > buffer_sizemax) then
                                  write(info,'(a)') "error in subroutine init_mpi_buffers: face in contact : size of buffer_gll_duplicate too small"
                                  call error_stop(info)
                               endif

                               buffer_gll_duplicate(ngll_duplicate) = num_gll

                           enddo
                       enddo
!
!----------------->put gll points from contact edge in the buffer 
                   else if(surf_num <= NFACE+NEDGE) then
                       do i = 1,IG_NGLL
                           select case(surf_num - NFACE)
                               case(1)
                                   num_gll = ig_hexa_gll_glonum(i,1,1,elt_num)
                               case(2)
                                   num_gll = ig_hexa_gll_glonum(IG_NGLL,i,1,elt_num)
                               case(3)
                                   num_gll = ig_hexa_gll_glonum(i,IG_NGLL,1,elt_num)
                               case(4)
                                   num_gll = ig_hexa_gll_glonum(1,i,1,elt_num)
                               case(5)
                                   num_gll = ig_hexa_gll_glonum(i,1,IG_NGLL,elt_num)
                               case(6)
                                   num_gll = ig_hexa_gll_glonum(IG_NGLL,i,IG_NGLL,elt_num)
                               case(7)
                                   num_gll = ig_hexa_gll_glonum(i,IG_NGLL,IG_NGLL,elt_num)
                               case(8)
                                   num_gll = ig_hexa_gll_glonum(1,i,IG_NGLL,elt_num)
                               case(9)
                                   num_gll = ig_hexa_gll_glonum(1,1,i,elt_num)
                               case(10)
                                   num_gll = ig_hexa_gll_glonum(IG_NGLL,1,i,elt_num)
                               case(11)
                                   num_gll = ig_hexa_gll_glonum(IG_NGLL,IG_NGLL,i,elt_num)
                               case(12)
                                   num_gll = ig_hexa_gll_glonum(1,IG_NGLL,i,elt_num)
                           end select

                           ngll_duplicate = ngll_duplicate + 1

                           if (ngll_duplicate > buffer_sizemax) then
                              write(info,'(a)') "error in subroutine init_mpi_buffers: edge in contact : size of buffer_gll_duplicate too small"
                              call error_stop(info)
                           endif

                           buffer_gll_duplicate(ngll_duplicate) = num_gll

                       enddo
!
!----------------->put gll points from contact corner nodes in the buffer 
                   else
                       select case(surf_num - (NFACE+NEDGE))
                           case(1)
                               num_gll = ig_hexa_gll_glonum(1,1,1,elt_num)
                           case(2)
                               num_gll = ig_hexa_gll_glonum(IG_NGLL,1,1,elt_num)
                           case(3)
                               num_gll = ig_hexa_gll_glonum(IG_NGLL,IG_NGLL,1,elt_num)
                           case(4)
                               num_gll = ig_hexa_gll_glonum(1,IG_NGLL,1,elt_num)
                           case(5)
                               num_gll = ig_hexa_gll_glonum(1,1,IG_NGLL,elt_num)
                           case(6)
                               num_gll = ig_hexa_gll_glonum(IG_NGLL,1,IG_NGLL,elt_num)
                           case(7)
                               num_gll = ig_hexa_gll_glonum(IG_NGLL,IG_NGLL,IG_NGLL,elt_num)
                           case(8)
                               num_gll = ig_hexa_gll_glonum(1,IG_NGLL,IG_NGLL,elt_num)
                       end select

                       ngll_duplicate = ngll_duplicate + 1

                       if (ngll_duplicate > buffer_sizemax) then
                          write(info,'(a)') "error in subroutine init_mpi_buffers: node in contact : size of buffer_gll_duplicate too small"
                          call error_stop(info)
                       endif

                       buffer_gll_duplicate(ngll_duplicate) = num_gll

                   endif
               endif
           enddo
!
!--------->remove duplicates
           call remove_duplicate(buffer_gll_duplicate,ngll_duplicate,tg_cpu_neighbor(icpu_neighbor)%gll_send,ngll_unique)
           
           if (mpi_buffer_sizemax < ngll_unique) mpi_buffer_sizemax = ngll_unique
           tg_cpu_neighbor(icpu_neighbor)%ngll = ngll_unique
!
!--------->ig_mpi_buffer_offset(x) -> give offset in rg_mpi_buffer for xth connected proc (start from 0 !!!)
           if (LG_ASYNC_MPI_COMM) then
              ig_mpi_buffer_offset(icpu_neighbor+1) = ig_mpi_buffer_offset(icpu_neighbor) + tg_cpu_neighbor(icpu_neighbor)%ngll*IG_NDOF
           endif

           if (LG_OUTPUT_DEBUG_FILE)  then
               open(unit=get_newunit(myunit_debug),file="debug."//trim(cg_prefix)//".mpi.buffer.cpu."//trim(cg_myrank))
               write(unit=myunit_debug,fmt='(3(a,i6),a)') "from proc ",ig_myrank," to proc ",icpu," : ",ngll_unique," gll points in MPI buffer"
               do i = 1,ngll_unique
                   write(unit=myunit_debug,fmt='(i10,3(e14.7,1x))') tg_cpu_neighbor(icpu_neighbor)%gll_send(i)&
                                                                   ,rg_gll_coordinate(1,tg_cpu_neighbor(icpu_neighbor)%gll_send(i))&
                                                                   ,rg_gll_coordinate(2,tg_cpu_neighbor(icpu_neighbor)%gll_send(i))&
                                                                   ,rg_gll_coordinate(3,tg_cpu_neighbor(icpu_neighbor)%gll_send(i))
               enddo
               close(myunit_debug)
           endif

       enddo

!
!----->check if cpus connected together send and receive the same number of gll
       do icpu_neighbor = 1,ig_ncpu_neighbor

           icpu      = tg_cpu_neighbor(icpu_neighbor)%icpu
           ngll_send = tg_cpu_neighbor(icpu_neighbor)%ngll 

           call mpi_sendrecv(ngll_send,1,MPI_INTEGER,icpu,89,ngll_recv,1,MPI_INTEGER,icpu,89,MPI_COMM_WORLD,statut,ios)

           if (ngll_send /= ngll_recv) then
              write(cl_rank,'(i6.6)') icpu
              write(info,'(a)') "error in subroutine init_mpi_buffers: different number ngll_send and ngll_recv betwen cpu"//trim(cg_myrank)//" and "//trim(cl_rank)
              call error_stop(info)
           endif

       enddo
       
       ig_mpi_buffer_sizemax = mpi_buffer_sizemax
       
       allocate(gll_coord_send(IG_NDOF*mpi_buffer_sizemax,ig_ncpu_neighbor),stat=ios)
       if (ios /= 0) then
          write(info,'(a)') "error in subroutine init_mpi_buffers while allocating gll_coord_send"
          call error_stop(info)
       endif

       allocate(gll_coord_recv(IG_NDOF*mpi_buffer_sizemax,ig_ncpu_neighbor),stat=ios)
       if (ios /= 0) then
          write(info,'(a)') "error in subroutine init_mpi_buffers while allocating gll_coord_recv"
          call error_stop(info)
       endif
!
!----->send/recv 2D buffer's glls coordinates for linked procs
       do icpu_neighbor = 1,ig_ncpu_neighbor

           icpu = tg_cpu_neighbor(icpu_neighbor)%icpu    

           !
           !fill coords to send
           do i = 1, tg_cpu_neighbor(icpu_neighbor)%ngll
              gll_coord_send((i-1)*IG_NDOF+1,icpu_neighbor) = rg_gll_coordinate(1,tg_cpu_neighbor(icpu_neighbor)%gll_send(i))
              gll_coord_send((i-1)*IG_NDOF+2,icpu_neighbor) = rg_gll_coordinate(2,tg_cpu_neighbor(icpu_neighbor)%gll_send(i))
              gll_coord_send((i-1)*IG_NDOF+3,icpu_neighbor) = rg_gll_coordinate(3,tg_cpu_neighbor(icpu_neighbor)%gll_send(i))
           enddo
           
           call mpi_sendrecv(gll_coord_send(1,icpu_neighbor),IG_NDOF*tg_cpu_neighbor(icpu_neighbor)%ngll,MPI_REAL,icpu,90,gll_coord_recv(1,icpu_neighbor),IG_NDOF*tg_cpu_neighbor(icpu_neighbor)%ngll,MPI_REAL,icpu,90,MPI_COMM_WORLD,statut,ios)

           if (ios /= 0) then
              write(info,'(a)') "error in subroutine init_mpi_buffers while sending gll coordinates"
              call error_stop(info)
           endif

       enddo
!
!----->build eqres
       do icpu_neighbor = 1,ig_ncpu_neighbor
          call create_gll_buffer_recv(gll_coord_send(1,icpu_neighbor),gll_coord_recv(1,icpu_neighbor),tg_cpu_neighbor(icpu_neighbor)%ngll,tg_cpu_neighbor(icpu_neighbor)%gll_recv,tg_cpu_neighbor(icpu_neighbor)%gll_send,mpi_buffer_sizemax)
       enddo
!
!----->deallocate arrays
       deallocate(gll_coord_send)
       deallocate(gll_coord_recv)
       deallocate(ig_cpu_neighbor_info)
!
!----->allocate array for mpi buffers
       if (LG_ASYNC_MPI_COMM) then
          allocate(rg_mpi_buffer_send(ig_mpi_buffer_offset(ig_ncpu_neighbor+1)))
          allocate(rg_mpi_buffer_recv(ig_mpi_buffer_offset(ig_ncpu_neighbor+1)))
       endif

       if (ig_myrank == 0) then
          write(IG_LST_UNIT,'(a)') "done"
          call flush(IG_LST_UNIT)
       endif

!***********************************************************************************************************************************************************************************
   end subroutine init_mpi_buffers
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine creates GLLs buffers received by cpu myrank from other cpus.
!!Identical GLLs between connected cpus are found by computing minimal distance among GLLs. 
!!Two GLLs from different cpus with the smallest minimal distance are supposed to be identical.
!>@param gll_coord_send : @f$x,y,z@f$-coordinates of GLLs in MPI_send buffers
!>@param gll_coord_recv : @f$x,y,z@f$-coordinates of GLLs in MPI_recv buffers
!>@param ngll        : number of unique common GLLs between cpu myrank and its connected cpus
!>@param gll_recv    : list of GLLs received by cpu myrank from its connected cpus
!>@param gll_send    : list of GLLs sent by cpu myrank to its connected cpus
!>@param max_size    : maximum buffer size for cpu myrank
!***********************************************************************************************************************************************************************************
   subroutine create_gll_buffer_recv(gll_coord_send, gll_coord_recv, ngll, gll_recv, gll_send, max_size)
!***********************************************************************************************************************************************************************************

       use mpi
       
       use mod_global_variables, only :&
                                       ig_ncpu_neighbor&
                                      ,ig_myrank&
                                      ,IG_NDOF&
                                      ,error_stop

       implicit none

       integer, intent(in)                                :: ngll
       integer, intent(in)                                :: max_size
       integer, intent(in)              , dimension(ngll) :: gll_send
       integer, intent(out), allocatable, dimension(:)    :: gll_recv

       real   , intent(in)                                :: gll_coord_send(IG_NDOF*max_size)
       real   , intent(in)                                :: gll_coord_recv(IG_NDOF*max_size)
                                                          
       integer                                            :: i
       integer                                            :: j
       integer                                            :: myindex
       integer                                            :: ios
       integer                                            :: imin
                                                          
       real                                               :: xs
       real                                               :: ys
       real                                               :: zs
       real                                               :: xr
       real                                               :: yr
       real                                               :: zr
       real                                               :: dist
       real                                               :: mindist
                                                          
       character(len=255)                                 :: info

       allocate(gll_recv(ngll),stat=ios)
       if (ios /= 0) then
          write(info,'(a)') "Error in subroutine create_gll_buffer_recv while allocating gll_recv"
          call error_stop(info)
       endif

ext:   do i = 1,ngll

           mindist = huge(dist)
           imin    = 0
           xr      = gll_coord_recv((i-1)*IG_NDOF+1)
           yr      = gll_coord_recv((i-1)*IG_NDOF+2)
           zr      = gll_coord_recv((i-1)*IG_NDOF+3)

int:       do j = 1,ngll

               xs = gll_coord_send((j-1)*IG_NDOF+1)
               ys = gll_coord_send((j-1)*IG_NDOF+2)
               zs = gll_coord_send((j-1)*IG_NDOF+3)

               dist = sqrt((xs-xr)**2 + (ys-yr)**2 + (zs-zr)**2)

               if ( dist == 0.0 ) then
                   gll_recv(i) = gll_send(j)
                   cycle ext
                elseif (dist < mindist) then 
                    mindist = dist
                    myindex = j
               endif

           enddo int

           gll_recv(i) = gll_send(myindex)

       enddo ext

       return

!***********************************************************************************************************************************************************************************
   end subroutine create_gll_buffer_recv
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine removes duplicated GLLs between cpu myrank and its connected cpus.
!>@param x1 : array with duplicated GLLs
!>@param n  : size of array x1
!>@param x3 : array without duplicated GLLs
!>@param m  : size of array x3
!***********************************************************************************************************************************************************************************
   subroutine remove_duplicate(x1,n,x3,m)
!***********************************************************************************************************************************************************************************
       
       implicit none
       
       integer, intent( in)                            :: n  
       integer, intent( in)             , dimension(n) :: x1 
       integer, intent(out), allocatable, dimension(:) :: x3 
       integer, intent(out)                            :: m 
                                                      
       integer                          , dimension(n) :: x2 !array without duplicate (size n)
       integer                                         :: i
       integer                                         :: j

       m     = 1
       x2(1) = x1(1)

outer: do i=2,n
           do j=1,m
               if (x2(j) == x1(i)) then
!
!----------------->Found a match so start looking again
                   cycle outer
               endif
           enddo
!
!--------->No match found so add it to array x2
           m      = m + 1
           x2(m)  = x1(i)
       enddo outer

       allocate(x3(m))
       x3(1:m) = x2(1:m)

       return

!***********************************************************************************************************************************************************************************
   end subroutine remove_duplicate
!***********************************************************************************************************************************************************************************

end module mod_init_mpibuffer
