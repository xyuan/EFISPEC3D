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
!!This file contains a module to assemble forces between cpu myrank and its neighbor cpus.

!>@brief
!!This module contains subroutines to assemble forces between cpu myrank and its neighbor cpus by synchrone or asynchrone MPI communications.
module mod_mpi_assemble

   use mpi

   implicit none

   public :: assemble_force_sync_comm
   public :: assemble_force_async_comm_init
   public :: assemble_force_async_comm_end

   contains

!
!
!>@brief subroutine to assemble forces between connected cpus by MPI synchrone communications.
!>@return mod_global_variables::rg_gll_acceleration
!***********************************************************************************************************************************************************************************
   subroutine assemble_force_sync_comm()
!***********************************************************************************************************************************************************************************
   
      use mpi
      
      use mod_global_variables, only :&
                                      tg_cpu_neighbor&
                                     ,rg_gll_acceleration&
                                     ,IG_NDOF&
                                     ,ig_ngll_total&
                                     ,ig_myrank&
                                     ,ig_mpi_buffer_sizemax&
                                     ,ig_ncpu_neighbor
      
      implicit none
      
      real   , dimension(IG_NDOF,ig_mpi_buffer_sizemax)                  :: dummyr
      real   , dimension(IG_NDOF,ig_mpi_buffer_sizemax,ig_ncpu_neighbor) :: dlacce
      
      integer, dimension(MPI_STATUS_SIZE)                                :: statut
      integer                                                            :: igll
      integer                                                            :: icon
      integer                                                            :: ngll
      integer                                                            :: ios
      !
      !
      !******************************************
      !fill buffer dlacce with gll to be send
      !******************************************
      do icon = 1,ig_ncpu_neighbor
         do igll = 1,tg_cpu_neighbor(icon)%ngll
   
            dlacce(1,igll,icon) = rg_gll_acceleration(1,tg_cpu_neighbor(icon)%gll_send(igll))
            dlacce(2,igll,icon) = rg_gll_acceleration(2,tg_cpu_neighbor(icon)%gll_send(igll))
            dlacce(3,igll,icon) = rg_gll_acceleration(3,tg_cpu_neighbor(icon)%gll_send(igll))
   
         enddo
      enddo
      !
      !
      !**********************************************************
      !synchrone comm: exchange buffers between cpu in contact
      !**********************************************************
      do icon = 1,ig_ncpu_neighbor
      
         ngll = tg_cpu_neighbor(icon)%ngll
         
         call mpi_sendrecv(dlacce(1,1,icon),IG_NDOF*ngll,mpi_real,tg_cpu_neighbor(icon)%icpu,300,dummyr(1,1),IG_NDOF*ngll,mpi_real,tg_cpu_neighbor(icon)%icpu,300,mpi_comm_world,statut,ios)
         
         do igll = 1,ngll
   
            rg_gll_acceleration(1,tg_cpu_neighbor(icon)%gll_recv(igll)) = rg_gll_acceleration(1,tg_cpu_neighbor(icon)%gll_recv(igll)) + dummyr(1,igll)
            rg_gll_acceleration(2,tg_cpu_neighbor(icon)%gll_recv(igll)) = rg_gll_acceleration(2,tg_cpu_neighbor(icon)%gll_recv(igll)) + dummyr(2,igll)
            rg_gll_acceleration(3,tg_cpu_neighbor(icon)%gll_recv(igll)) = rg_gll_acceleration(3,tg_cpu_neighbor(icon)%gll_recv(igll)) + dummyr(3,igll)
   
         enddo
      
      enddo
   
      return
!***********************************************************************************************************************************************************************************
   end subroutine assemble_force_sync_comm
!***********************************************************************************************************************************************************************************

!
!
!>@brief subroutine to initialize forces assembly between connected cpus by MPI asynchrone communications.
!***********************************************************************************************************************************************************************************
   subroutine assemble_force_async_comm_init()
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only :&
                                      tg_cpu_neighbor&
                                     ,rg_gll_acceleration&
                                     ,IG_NDOF&
                                     ,ig_ngll_total&
                                     ,ig_myrank&
                                     ,ig_mpi_buffer_offset&
                                     ,rg_mpi_buffer_send&
                                     ,rg_mpi_buffer_recv&
                                     ,ig_mpi_request_send&
                                     ,ig_mpi_request_recv&
                                     ,ig_mpi_buffer_sizemax&
                                     ,ig_ncpu_neighbor
   
      implicit none
   
      real, dimension(IG_NDOF,ig_mpi_buffer_sizemax,ig_ncpu_neighbor) :: dlacce
      integer                                                         :: igll
      integer                                                         :: icon
      integer                                                         :: ngll
      integer                                                         :: ios
      integer                                                         :: dof_offset
      integer                                                         :: proc_offset   
      !
      !
      !*****************************************************************************************************************************************
      !fill buffer dlacce with gll to be send
      !*****************************************************************************************************************************************
      do icon = 1,ig_ncpu_neighbor
         do igll =1,tg_cpu_neighbor(icon)%ngll
   
             dlacce(1,igll,icon) = rg_gll_acceleration(1,tg_cpu_neighbor(icon)%gll_send(igll))
             dlacce(2,igll,icon) = rg_gll_acceleration(2,tg_cpu_neighbor(icon)%gll_send(igll))
             dlacce(3,igll,icon) = rg_gll_acceleration(3,tg_cpu_neighbor(icon)%gll_send(igll))
   
         enddo
      enddo
      !
      !
      !*****************************************************************************************************************************************
      !asynchrone comm: exchange buffers between cpu in contact
      !*****************************************************************************************************************************************
      do icon = 1,ig_ncpu_neighbor
   
         ngll        = tg_cpu_neighbor(icon)%ngll
         proc_offset = ig_mpi_buffer_offset(icon)
   
         !
         !fill buffer for cpu myrank
         do igll = 1,ngll
   
            !beginning of x buffer = 0
            dof_offset = 0
            rg_mpi_buffer_send(proc_offset+dof_offset+igll) = dlacce(1,igll,icon)
   
            !beginning of y buffer = ngll
            dof_offset = ngll
            rg_mpi_buffer_send(proc_offset+dof_offset+igll) = dlacce(2,igll,icon)
   
            !beginning of z buffer = 2*ngll
            dof_offset = 2*ngll
            rg_mpi_buffer_send(proc_offset+dof_offset+igll) = dlacce(3,igll,icon)
   
         enddo
      
         call mpi_Irecv(rg_mpi_buffer_recv(proc_offset+1),IG_NDOF*ngll,mpi_real,tg_cpu_neighbor(icon)%icpu,300,mpi_comm_world,ig_mpi_request_recv(icon),ios)
         call mpi_Isend(rg_mpi_buffer_send(proc_offset+1),IG_NDOF*ngll,mpi_real,tg_cpu_neighbor(icon)%icpu,300,mpi_comm_world,ig_mpi_request_send(icon),ios)
   
      enddo
   
      return
!***********************************************************************************************************************************************************************************
   end subroutine assemble_force_async_comm_init
!***********************************************************************************************************************************************************************************

!
!
!>@brief subroutine to finalize forces assembly between connected cpus by MPI asynchrone communications.
!>@return mod_global_variables::rg_gll_acceleration
!***********************************************************************************************************************************************************************************
   subroutine assemble_force_async_comm_end()
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only :&
                                      tg_cpu_neighbor&
                                     ,rg_gll_acceleration&
                                     ,IG_NDOF&
                                     ,ig_ngll_total&
                                     ,ig_myrank&
                                     ,ig_mpi_buffer_offset&
                                     ,rg_mpi_buffer_send&
                                     ,rg_mpi_buffer_recv&
                                     ,ig_mpi_request_send&
                                     ,ig_mpi_request_recv&
                                     ,ig_mpi_buffer_sizemax&
                                     ,ig_ncpu_neighbor
   
      implicit none
   
      integer                             :: igll
      integer                             :: ngll
      integer                             :: ios
      integer                             :: dof_offset
      integer                             :: proc_offset
      integer                             :: nb_msg
      integer                             :: r_index
      integer, dimension(MPI_STATUS_SIZE) :: statut
      
      nb_msg = 0
   
      do while(nb_msg < ig_ncpu_neighbor)
   !
   !---->wait for any received message
         call mpi_waitany(ig_ncpu_neighbor, ig_mpi_request_recv, r_index, statut, ios)
   
         nb_msg      = nb_msg + 1   
         ngll        = tg_cpu_neighbor(r_index)%ngll
         proc_offset = ig_mpi_buffer_offset(r_index)
   !
   !---->sum contributions for this proc
   
         !beginning of x buffer = 0
         dof_offset = 0
         do igll = 1,ngll
            rg_gll_acceleration(1,tg_cpu_neighbor(r_index)%gll_recv(igll)) = rg_gll_acceleration(1,tg_cpu_neighbor(r_index)%gll_recv(igll)) + rg_mpi_buffer_recv(proc_offset+dof_offset+igll)
         enddo
   
         !beginning of y buffer = ngll
         dof_offset = ngll
         do igll = 1,ngll
            rg_gll_acceleration(2,tg_cpu_neighbor(r_index)%gll_recv(igll)) = rg_gll_acceleration(2,tg_cpu_neighbor(r_index)%gll_recv(igll)) + rg_mpi_buffer_recv(proc_offset+dof_offset+igll)
         enddo
   
         !beginning of z buffer = 2*ngll
         dof_offset = 2*ngll
         do igll = 1,ngll
            rg_gll_acceleration(3,tg_cpu_neighbor(r_index)%gll_recv(igll)) = rg_gll_acceleration(3,tg_cpu_neighbor(r_index)%gll_recv(igll)) + rg_mpi_buffer_recv(proc_offset+dof_offset+igll)
         enddo
   
      enddo

      return
!***********************************************************************************************************************************************************************************
   end subroutine assemble_force_async_comm_end
!***********************************************************************************************************************************************************************************

end module mod_mpi_assemble
