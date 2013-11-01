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
!!This file contains a module to allocate arrays and to compute an approximation of the total RAM used by EFISPEC3D.

!>@brief
!!This module contains subroutines to allocate arrays and to compute an approximation of the total RAM used by EFISPEC3D.
module mod_init_memory

   implicit none

   private

   public :: init_array_real
   public :: init_array_int 
   public :: memory_consumption

!>@brief
!!Interface init_array_real to redirect allocation to n-rank arrays.
   interface init_array_real

      module procedure init_array_rank1_real
      module procedure init_array_rank2_real
      module procedure init_array_rank3_real
      module procedure init_array_rank4_real
      module procedure init_array_rank5_real

   end interface init_array_real


!>@brief
!!Interface init_array_int to redirect allocation to n-rank arrays.
   interface init_array_int 

      module procedure init_array_rank1_int 
      module procedure init_array_rank2_int 
      module procedure init_array_rank3_int 
      module procedure init_array_rank4_int 

   end interface init_array_int


   contains

!
!
!>@brief This subroutine allocates real-value array of rank 1
!>@param t     : real-value array to be allocated
!>@param n1    : size of the array
!>@param tname : name of the array
!***********************************************************************************************************************************************************************************
   function init_array_rank1_real(t,n1,tname) result(err)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : error_stop

      implicit none

      integer         , intent( in)                            :: n1
      real            , intent(out), allocatable, dimension(:) :: t
      character(len=*), intent( in)                            :: tname

      integer :: err
      integer :: i

      allocate(t(n1),stat=err)

      if (err /= 0) then

         call error_stop("error while allocating array "//trim(adjustl(tname)))

      else

         do i = 1,n1
            t(i) = 0.0
         enddo
         
      endif
 
      return
!***********************************************************************************************************************************************************************************
   end function init_array_rank1_real
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine allocates real-value array of rank 2
!>@param t  : real-value array to be allocated --> t(n2,n1)
!>@param n1 : size of the dimension 1 of the array
!>@param n2 : size of the dimension 2 of the array
!>@param tname : name of the array
!***********************************************************************************************************************************************************************************
   function init_array_rank2_real(t,n1,n2,tname) result(err)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : error_stop

      implicit none

      integer         , intent( in)                              :: n1
      integer         , intent( in)                              :: n2
      real            , intent(out), allocatable, dimension(:,:) :: t
      character(len=*), intent( in)                              :: tname

      integer :: err
      integer :: i
      integer :: j

      allocate(t(n2,n1),stat=err)

      if (err /= 0) then

         call error_stop("error while allocating array "//trim(adjustl(tname)))

      else

         do i = 1,n1
            do j = 1,n2
               t(j,i) = 0.0
            enddo
         enddo
         
      endif
 
      return
!***********************************************************************************************************************************************************************************
   end function init_array_rank2_real
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine allocates real-value array of rank 3
!>@param t  : real-value array to be allocated --> t(n3,n2,n1)
!>@param n1 : size of the dimension 1 of the array
!>@param n2 : size of the dimension 2 of the array
!>@param n3 : size of the dimension 3 of the array
!>@param tname : name of the array
!***********************************************************************************************************************************************************************************
   function init_array_rank3_real(t,n1,n2,n3,tname) result(err)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : error_stop

      implicit none

      integer         , intent( in)                                :: n1
      integer         , intent( in)                                :: n2
      integer         , intent( in)                                :: n3
      real            , intent(out), allocatable, dimension(:,:,:) :: t
      character(len=*), intent( in)                                :: tname

      integer :: err
      integer :: i
      integer :: j
      integer :: k

      allocate(t(n3,n2,n1),stat=err)

      if (err /= 0) then

         call error_stop("error while allocating array "//trim(adjustl(tname)))

      else

         do i = 1,n1
            do j = 1,n2
               do k = 1,n3
                  t(k,j,i) = 0.0
               enddo
            enddo
         enddo
         
      endif
 
      return
!***********************************************************************************************************************************************************************************
   end function init_array_rank3_real
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine allocates real-value array of order 4
!>@param t  : real-value array to be allocated
!>@param n1 : size of the dimension 1 of the array
!>@param n2 : size of the dimension 2 of the array
!>@param n3 : size of the dimension 3 of the array
!>@param n4 : size of the dimension 4 of the array
!>@param tname : name of the array
!***********************************************************************************************************************************************************************************
   function init_array_rank4_real(t,n1,n2,n3,n4,tname) result(err)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : error_stop

      implicit none

      integer         , intent( in)                                  :: n1
      integer         , intent( in)                                  :: n2
      integer         , intent( in)                                  :: n3
      integer         , intent( in)                                  :: n4
      real            , intent(out), allocatable, dimension(:,:,:,:) :: t
      character(len=*), intent(in)                                   :: tname

      integer :: err
      integer :: i
      integer :: j
      integer :: k
      integer :: l

      allocate(t(n4,n3,n2,n1),stat=err)

      if (err /= 0) then

         call error_stop("error while allocating array "//trim(adjustl(tname)))

      else

         do i = 1,n1
            do j = 1,n2
               do k = 1,n3
                  do l = 1,n4
                     t(l,k,j,i) = 0.0
                  enddo
               enddo
            enddo
         enddo
         
      endif
 
      return
!***********************************************************************************************************************************************************************************
   end function init_array_rank4_real
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine allocates real-value array of order 5
!>@param t  : real-value array to be allocated
!>@param n1 : size of the dimension 1 of the array
!>@param n2 : size of the dimension 2 of the array
!>@param n3 : size of the dimension 3 of the array
!>@param n4 : size of the dimension 4 of the array
!>@param n5 : size of the dimension 5 of the array
!>@param tname : name of the array
!***********************************************************************************************************************************************************************************
   function init_array_rank5_real(t,n1,n2,n3,n4,n5,tname) result(err)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : error_stop

      implicit none

      integer         , intent( in)                                    :: n1
      integer         , intent( in)                                    :: n2
      integer         , intent( in)                                    :: n3
      integer         , intent( in)                                    :: n4
      integer         , intent( in)                                    :: n5
      real            , intent(out), allocatable, dimension(:,:,:,:,:) :: t
      character(len=*), intent(in)                                     :: tname

      integer :: err
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: m

      allocate(t(n5,n4,n3,n2,n1),stat=err)

      if (err /= 0) then

         call error_stop("error while allocating array "//trim(adjustl(tname)))

      else

         do i = 1,n1
            do j = 1,n2
               do k = 1,n3
                  do l = 1,n4
                     do m = 1,n5
                        t(m,l,k,j,i) = 0.0
                     enddo
                  enddo
               enddo
            enddo
         enddo
         
      endif
 
      return
!***********************************************************************************************************************************************************************************
   end function init_array_rank5_real
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine allocates integer-value array of rank 1
!>@param t     : integer-value array to be allocated
!>@param n1    : size of the array
!>@param tname : name of the array
!***********************************************************************************************************************************************************************************
   function init_array_rank1_int(t,n1,tname) result(err)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : error_stop

      implicit none

      integer         , intent( in)                            :: n1
      integer         , intent(out), allocatable, dimension(:) :: t
      character(len=*), intent( in)                            :: tname

      integer :: err
      integer :: i

      allocate(t(n1),stat=err)

      if (err /= 0) then

         call error_stop("error while allocating array "//trim(adjustl(tname)))

      else

         do i = 1,n1
            t(i) = 0
         enddo
         
      endif
 
      return
!***********************************************************************************************************************************************************************************
   end function init_array_rank1_int
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine allocates integer-value array of rank 2
!>@param t     : integer-value array to be allocated
!>@param n1    : size of dimension 1 of the array
!>@param n2    : size of dimension 2 of the array
!>@param tname : name of the array
!***********************************************************************************************************************************************************************************
   function init_array_rank2_int(t,n1,n2,tname) result(err)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : error_stop

      implicit none

      integer         , intent( in)                              :: n1
      integer         , intent( in)                              :: n2
      integer         , intent(out), allocatable, dimension(:,:) :: t
      character(len=*), intent( in)                              :: tname

      integer :: err
      integer :: i
      integer :: j

      allocate(t(n2,n1),stat=err)

      if (err /= 0) then

         call error_stop("error while allocating array "//trim(adjustl(tname)))

      else

         do i = 1,n1
            do j = 1,n2
               t(j,i) = 0
            enddo
         enddo
         
      endif
 
      return
!***********************************************************************************************************************************************************************************
   end function init_array_rank2_int
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine allocates int-value array of rank 3
!>@param t  : int-value array to be allocated --> t(n3,n2,n1)
!>@param n1 : size of the dimension 1 of the array
!>@param n2 : size of the dimension 2 of the array
!>@param n3 : size of the dimension 3 of the array
!>@param tname : name of the array
!***********************************************************************************************************************************************************************************
   function init_array_rank3_int(t,n1,n2,n3,tname) result(err)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : error_stop

      implicit none

      integer         , intent( in)                               :: n1
      integer         , intent( in)                               :: n2
      integer         , intent( in)                               :: n3
      integer        , intent(out), allocatable, dimension(:,:,:) :: t
      character(len=*), intent( in)                               :: tname

      integer :: err
      integer :: i
      integer :: j
      integer :: k

      allocate(t(n3,n2,n1),stat=err)

      if (err /= 0) then

         call error_stop("error while allocating array "//trim(adjustl(tname)))

      else

         do i = 1,n1
            do j = 1,n2
               do k = 1,n3
                  t(k,j,i) = 0
               enddo
            enddo
         enddo
         
      endif
 
      return
!***********************************************************************************************************************************************************************************
   end function init_array_rank3_int
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine allocates int-value array of order 4
!>@param t  : int-value array to be allocated
!>@param n1 : size of the dimension 1 of the array
!>@param n2 : size of the dimension 2 of the array
!>@param n3 : size of the dimension 3 of the array
!>@param n4 : size of the dimension 4 of the array
!>@param tname : name of the array
!***********************************************************************************************************************************************************************************
   function init_array_rank4_int(t,n1,n2,n3,n4,tname) result(err)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : error_stop

      implicit none

      integer         , intent( in)                                 :: n1
      integer         , intent( in)                                 :: n2
      integer         , intent( in)                                 :: n3
      integer         , intent( in)                                 :: n4
      integer        , intent(out), allocatable, dimension(:,:,:,:) :: t
      character(len=*), intent(in)                                  :: tname

      integer :: err
      integer :: i
      integer :: j
      integer :: k
      integer :: l

      allocate(t(n4,n3,n2,n1),stat=err)

      if (err /= 0) then

         call error_stop("error while allocating array "//trim(adjustl(tname)))

      else

         do i = 1,n1
            do j = 1,n2
               do k = 1,n3
                  do l = 1,n4
                     t(l,k,j,i) = 0
                  enddo
               enddo
            enddo
         enddo
         
      endif
 
      return
!***********************************************************************************************************************************************************************************
   end function init_array_rank4_int
!***********************************************************************************************************************************************************************************

!
!
!>@brief Subroutine to compute an approximation of total RAM used by global variables of EFISPEC3D.
!!See module mod_global_variables
!***********************************************************************************************************************************************************************************
   subroutine memory_consumption()
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables
      
      implicit none
      
      real, parameter :: MO  = (1024.0)**2
      real            :: total_memory
      real            :: total_memory_all_cpu
      real            :: size_cpu_world(ig_ncpu)
      integer         :: ios
      integer         :: i
      
      total_memory         = 0.0
      total_memory_all_cpu = 0.0

      total_memory = total_memory + real(sizeof(LG_ASYNC_MPI_COMM))/MO
      total_memory = total_memory + real(sizeof(LG_OUTPUT_CPUTIME))/MO
      total_memory = total_memory + real(sizeof(LG_OUTPUT_DEBUG_FILE))/MO
      total_memory = total_memory + real(sizeof(LG_SNAPSHOT_VTK))/MO
      total_memory = total_memory + real(sizeof(LG_SNAPSHOT_GMT))/MO
      total_memory = total_memory + real(sizeof(LG_OUTPUT_MEDIUM_VTK))/MO
      total_memory = total_memory + real(sizeof(LG_VISCO))/MO
      total_memory = total_memory + real(sizeof(IG_LAGRANGE_ORDER))/MO
      total_memory = total_memory + real(sizeof(IG_NGLL))/MO
      total_memory = total_memory + real(sizeof(IG_LST_UNIT))/MO
      total_memory = total_memory + real(sizeof(IG_NDOF))/MO
      total_memory = total_memory + real(sizeof(IG_NRELAX))/MO
      total_memory = total_memory + real(sizeof(RG_NEWMARK_GAMMA))/MO
      total_memory = total_memory + real(sizeof(RG_PI))/MO
      
      total_memory = total_memory + real(sizeof(rg_dcsource_gll_force))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_jacobian_det))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_dxidx))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_dxidy))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_dxidz))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_detdx))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_detdy))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_detdz))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_dzedx))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_dzedy))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_dzedz))/MO
      total_memory = total_memory + real(sizeof(rg_quadp_gll_normal))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_rho))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_rhovs2))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_rhovp2))/MO
      total_memory = total_memory + real(sizeof(rg_quadp_gll_rhovs))/MO
      total_memory = total_memory + real(sizeof(rg_quadp_gll_rhovp))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_wkqs))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_wkqp))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_ksixx))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_ksiyy))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_ksizz))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_ksixy))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_ksixz))/MO
      total_memory = total_memory + real(sizeof(rg_hexa_gll_ksiyz))/MO
      total_memory = total_memory + real(sizeof(RG_RELAX_COEFF))/MO
      total_memory = total_memory + real(sizeof(rg_mem_var_exp))/MO
      total_memory = total_memory + real(sizeof(rg_quadp_gll_jaco_det))/MO
      total_memory = total_memory + real(sizeof(rg_gll_displacement))/MO
      total_memory = total_memory + real(sizeof(rg_gll_velocity))/MO
      total_memory = total_memory + real(sizeof(rg_gll_acceleration))/MO
      total_memory = total_memory + real(sizeof(rg_gll_acctmp))/MO
      total_memory = total_memory + real(sizeof(rg_gnode_abscissa))/MO
      total_memory = total_memory + real(sizeof(rg_gnode_abscissa_dist))/MO
      total_memory = total_memory + real(sizeof(rg_gll_mass_matrix))/MO
      total_memory = total_memory + real(sizeof(rg_gnode_x))/MO
      total_memory = total_memory + real(sizeof(rg_gnode_y))/MO
      total_memory = total_memory + real(sizeof(rg_gnode_z))/MO
      total_memory = total_memory + real(sizeof(rg_receiver_snapshot_z))/MO
      total_memory = total_memory + real(sizeof(ig_receiver_snapshot_locnum))/MO
      total_memory = total_memory + real(sizeof(rg_dcsource_user_func))/MO
      total_memory = total_memory + real(sizeof(rg_sfsource_user_func))/MO
      total_memory = total_memory + real(sizeof(rg_mpi_buffer_send))/MO
      total_memory = total_memory + real(sizeof(rg_mpi_buffer_recv))/MO
      total_memory = total_memory + real(sizeof(rg_dt))/MO
      total_memory = total_memory + real(sizeof(rg_dt2))/MO
      total_memory = total_memory + real(sizeof(rg_simu_current_time))/MO
      total_memory = total_memory + real(sizeof(rg_simu_total_time))/MO
      total_memory = total_memory + real(sizeof(rg_mesh_xmax))/MO
      total_memory = total_memory + real(sizeof(rg_mesh_xmin))/MO
      total_memory = total_memory + real(sizeof(rg_mesh_ymax))/MO
      total_memory = total_memory + real(sizeof(rg_mesh_ymin))/MO
      total_memory = total_memory + real(sizeof(rg_mesh_zmax))/MO
      total_memory = total_memory + real(sizeof(rg_mesh_zmin))/MO
      total_memory = total_memory + real(sizeof(rg_receiver_snapshot_dxdy))/MO
      total_memory = total_memory + real(sizeof(rg_receiver_snapshot_nx))/MO
      total_memory = total_memory + real(sizeof(rg_receiver_snapshot_ny))/MO
      total_memory = total_memory + real(sizeof(ig_hexa_gll_glonum))/MO
      total_memory = total_memory + real(sizeof(ig_quadp_gll_glonum))/MO
      total_memory = total_memory + real(sizeof(ig_quadf_gll_glonum))/MO
      total_memory = total_memory + real(sizeof(ig_hexa_gnode_glonum))/MO
      total_memory = total_memory + real(sizeof(ig_quadp_gnode_glonum))/MO
      total_memory = total_memory + real(sizeof(ig_quadf_gnode_glonum))/MO
      total_memory = total_memory + real(sizeof(ig_hexa_gnode_xiloc))/MO
      total_memory = total_memory + real(sizeof(ig_hexa_gnode_etloc))/MO
      total_memory = total_memory + real(sizeof(ig_hexa_gnode_zeloc))/MO
      total_memory = total_memory + real(sizeof(ig_quad_gnode_xiloc))/MO
      total_memory = total_memory + real(sizeof(ig_quad_gnode_etloc))/MO
      total_memory = total_memory + real(sizeof(ig_hexa_receiver_unit))/MO
      total_memory = total_memory + real(sizeof(ig_mpi_request_send))/MO
      total_memory = total_memory + real(sizeof(ig_mpi_request_recv))/MO
      total_memory = total_memory + real(sizeof(ig_mpi_buffer_offset))/MO
      total_memory = total_memory + real(sizeof(ig_hexa_material_number))/MO
      total_memory = total_memory + real(sizeof(ig_material_type))/MO
      total_memory = total_memory + real(sizeof(ig_receiver_snapshot_glonum))/MO
      total_memory = total_memory + real(sizeof(ig_receiver_snapshot_mpi_shift))/MO
      total_memory = total_memory + real(sizeof(ig_receiver_snapshot_total_number))/MO
      total_memory = total_memory + real(sizeof(ig_ncpu))/MO
      total_memory = total_memory + real(sizeof(ig_myrank))/MO
      total_memory = total_memory + real(sizeof(ig_ncpu_neighbor))/MO
      total_memory = total_memory + real(sizeof(ig_nhexa))/MO
      total_memory = total_memory + real(sizeof(ig_nhexa_outer))/MO
      total_memory = total_memory + real(sizeof(ig_nhexa_inner))/MO
      total_memory = total_memory + real(sizeof(ig_nquad_parax))/MO
      total_memory = total_memory + real(sizeof(ig_nquad_fsurf))/MO
      total_memory = total_memory + real(sizeof(ig_mesh_nnode))/MO
      total_memory = total_memory + real(sizeof(ig_hexa_nnode))/MO
      total_memory = total_memory + real(sizeof(ig_quad_nnode))/MO
      total_memory = total_memory + real(sizeof(ig_line_nnode))/MO
      total_memory = total_memory + real(sizeof(ig_hexa_node2gll))/MO
      total_memory = total_memory + real(sizeof(ig_ndt))/MO
      total_memory = total_memory + real(sizeof(ig_idt))/MO
      total_memory = total_memory + real(sizeof(ig_receiver_saving_incr))/MO
      total_memory = total_memory + real(sizeof(ig_snapshot_saving_incr))/MO
      total_memory = total_memory + real(sizeof(lg_boundary_absorption))/MO
      total_memory = total_memory + real(sizeof(ig_ndcsource))/MO
      total_memory = total_memory + real(sizeof(ig_nsfsource))/MO
      total_memory = total_memory + real(sizeof(ig_nreceiver_hexa))/MO
      total_memory = total_memory + real(sizeof(ig_nreceiver_snapshot))/MO
      total_memory = total_memory + real(sizeof(ig_medium_type))/MO
      total_memory = total_memory + real(sizeof(ig_nmaterial))/MO
      total_memory = total_memory + real(sizeof(ig_mpi_buffer_sizemax))/MO
      total_memory = total_memory + real(sizeof(ig_quadp_neighbor_hexa))/MO
      total_memory = total_memory + real(sizeof(ig_quadp_neighbor_hexaface))/MO
      total_memory = total_memory + real(sizeof(ig_cpu_neighbor_info))/MO
      total_memory = total_memory + real(sizeof(ig_nneighbor_all_kind))/MO
      total_memory = total_memory + real(sizeof(ig_cpu_name_len))/MO
      total_memory = total_memory + real(sizeof(cg_cpu_name))/MO
      total_memory = total_memory + real(sizeof(cg_prefix))/MO
      total_memory = total_memory + real(sizeof(cg_myrank))/MO
      total_memory = total_memory + real(sizeof(cg_ncpu))/MO
      total_memory = total_memory + real(sizeof(lg_snapshot))/MO
      total_memory = total_memory + real(sizeof(lg_snapshot_displacement))/MO
      total_memory = total_memory + real(sizeof(lg_snapshot_velocity))/MO
      total_memory = total_memory + real(sizeof(lg_snapshot_acceleration))/MO

      total_memory = total_memory + real(sizeof(rg_gll_weight))/MO
      total_memory = total_memory + real(sizeof(rg_gll_abscissa))/MO
      total_memory = total_memory + real(sizeof(rg_gll_lagrange_deriv))/MO
      total_memory = total_memory + real(sizeof(rg_gll_abscissa_dist))/MO
      total_memory = total_memory + real(sizeof(ig_ngll_total))/MO
      total_memory = total_memory + real(sizeof(tg_sfsource))/MO
      total_memory = total_memory + real(sizeof(tg_dcsource))/MO
      total_memory = total_memory + real(sizeof(tg_receiver_hexa))/MO
      total_memory = total_memory + real(sizeof(tg_receiver_snapshot_quad))/MO
      total_memory = total_memory + real(sizeof(rg_gll_coordinate))/MO
      total_memory = total_memory + real(sizeof(tg_cpu_neighbor))/MO
      total_memory = total_memory + real(sizeof(tg_elastic_material))/MO
      total_memory = total_memory + real(sizeof(tg_viscoelastic_material))/MO

      !
      !cpu 0 gather all info
      call mpi_gather(total_memory,1,mpi_real,size_cpu_world,1,mpi_real,0,mpi_comm_world,ios)

      !
      !cpu 0 write info in *.lst
      if (ig_myrank == 0) then
      
         write(IG_LST_UNIT,'(" ",/,a)') "Memory consumption approximation"
      
         do i = 1,ig_ncpu
            write(IG_LST_UNIT,'(a,i8,a,f15.3,a)') "cpu ",i-1," : ",size_cpu_world(i)," Mo"
            total_memory_all_cpu = total_memory_all_cpu + size_cpu_world(i)
         enddo
      
         write(IG_LST_UNIT,'(a,9x,f15.3,a)') "Total ",total_memory_all_cpu," Mo"
      
      endif
      
      return

!***********************************************************************************************************************************************************************************
   end subroutine memory_consumption
!***********************************************************************************************************************************************************************************

end module mod_init_memory
