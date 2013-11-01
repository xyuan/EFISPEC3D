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
!!This file contains a module write snapshots readable by GMT or ParaView.

!>@brief
!!This module contains subroutines to compute and write snapshots of the free surface (either in GMT or VTK xml formats).
module mod_snapshot

   use mpi

   implicit none

   public  :: init_snapshot
   public  :: write_snapshot
   private :: write_snapshot_gmt_native_format
   private :: write_header_gmt_native_format
   private :: write_snapshot_vtk
   public  :: write_collection_vtk
   private :: collection_vtk
   public  :: write_peak_ground_motion

   contains

!
!
!>@brief
!!This subroutine generates a structured grid of receivers on the free surface.
!>@return mod_global_variables::tg_receiver_snapshot_quad : information about receivers' used for snapshots in cpu myrank. See mod_global_variables::type_receiver_quad
!>@return mod_global_variables::rg_receiver_snapshot_z : @f$z@f$-coordinate of receivers' used for snapshots
!***********************************************************************************************************************************************************************************
   subroutine init_snapshot()
!***********************************************************************************************************************************************************************************

      use mpi
      
      use mod_global_variables, only :&
                                      tg_receiver_snapshot_quad&
                                     ,ig_quadf_gll_glonum&
                                     ,ig_nreceiver_snapshot&
                                     ,ig_receiver_snapshot_glonum&
                                     ,rg_receiver_snapshot_z&
                                     ,ig_receiver_snapshot_locnum&
                                     ,ig_receiver_snapshot_total_number&
                                     ,ig_receiver_snapshot_mpi_shift&
                                     ,IG_LST_UNIT&
                                     ,IG_NGLL&
                                     ,ig_myrank&
                                     ,error_stop&
                                     ,ig_ncpu&
                                     ,rg_receiver_snapshot_nx&
                                     ,rg_receiver_snapshot_ny&
                                     ,rg_mesh_xmin&
                                     ,rg_mesh_xmax&
                                     ,rg_mesh_ymin&
                                     ,rg_mesh_ymax&
                                     ,rg_receiver_snapshot_dxdy
     
      use mod_receiver        , only :&
                                      compute_info_quad_receiver&
                                     ,search_closest_quad_gll

      implicit none
      
      real   , dimension(ig_ncpu) :: rldum1
      real   , allocatable, dimension(:)   :: rldum2
      real   , allocatable, dimension(:)   :: rldum3
      real                                 :: rec_x
      real                                 :: rec_y
      real                                 :: rec_dmin

                                           
      integer, allocatable, dimension(:)   :: ildum1
      integer             , dimension(1)   :: min_loc
      integer                              :: nx
      integer                              :: ny
      integer                              :: ix
      integer                              :: iy
      integer                              :: iz
   
      integer                              :: rec_iel
      integer                              :: rec_lgll
      integer                              :: rec_mgll
      integer                              :: irec
      integer                              :: jrec
      integer                              :: nrec
   
      integer                              :: icpu
      integer                              :: ios
   
      logical, allocatable, dimension(:)   :: is_rec_in_cpu
   
      character(len=255)                   :: info
   
      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(" ",/,a)') "searching receivers for snapshot..."
         call flush(IG_LST_UNIT)
      endif
! 
!---->set up the regular grid of receivers at the free surface
      nx = int((rg_mesh_xmax - rg_mesh_xmin)/rg_receiver_snapshot_dxdy)+1
      ny = int((rg_mesh_ymax - rg_mesh_ymin)/rg_receiver_snapshot_dxdy)+1
   
      rg_receiver_snapshot_nx = nx
      rg_receiver_snapshot_ny = ny
   
      allocate(is_rec_in_cpu(nx*ny),stat=ios)
      if (ios /= 0) then
         write(info,'(a)') "error in subroutine init_snapshot while allocating is_rec_in_cpu"
         call error_stop(info)
      else
         do ix = 1,nx*ny
            is_rec_in_cpu(ix) = .false.
         enddo
      endif
   
      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(3(a,i0))') " -->number of receivers for snapshot = ",nx,"*",ny," = ",nx * ny
         call flush(IG_LST_UNIT)
      endif
! 
!---->first pass on nx*ny receivers to count the number of receivers that belong to cpu 'myrank'
      irec = 0
      nrec = 0
   
      do iy = 1,ny
     
         do ix = 1,nx
    
            irec  = irec + 1
            rec_x = rg_mesh_xmin + real(ix-1)*rg_receiver_snapshot_dxdy
            rec_y = rg_mesh_ymax - real(iy-1)*rg_receiver_snapshot_dxdy
! 
!---------->find the closest gll point to the receiver irec and find to which free surface quadrangle element the receiver belongs
            call search_closest_quad_gll(rec_x                       &
                                                     ,rec_y                       &
                                                     ,ig_quadf_gll_glonum &
                                                     ,rec_dmin                    &
                                                     ,rec_lgll                    &
                                                     ,rec_mgll                    &
                                                     ,rec_iel                     )
   
            call mpi_allgather(rec_dmin,1,mpi_real,rldum1,1,mpi_real,mpi_comm_world,ios)
            min_loc = minloc(rldum1(1:ig_ncpu))
   
            if (min_loc(1) == ig_myrank+1) then
               nrec                = nrec + 1
               is_rec_in_cpu(irec) = .true.
            endif
  
            if ( (ig_myrank == 0) .and. ( (irec == 1)  .or. (mod(irec,5000) == 0) .or. (irec == nx*ny) ) ) then
               write(IG_LST_UNIT,'(i10,a)') irec," receivers found"
               call flush(IG_LST_UNIT)
            endif
   
         enddo
     
      enddo
! 
!---->allocate array structure tg_receiver_snapshot_quad
      ig_nreceiver_snapshot = nrec
   
      allocate(tg_receiver_snapshot_quad(ig_nreceiver_snapshot),stat=ios)
      if (ios /= 0) then
         write(info,'(a)') "error in subroutine init_snapshot while allocating tg_receiver_snapshot_quad"
         call error_stop(info)
      endif
! 
!---->second pass on receiver that belong to cpu 'myrank' and fill array structure tg_receiver_snapshot_quad
      irec = 0
      jrec = 0
   
      do iy = 1,ny
         do ix = 1,nx
   
            irec  = irec + 1
   
            if (is_rec_in_cpu(irec)) then
   
               jrec  = jrec + 1
               rec_x = rg_mesh_xmin + real(ix-1)*rg_receiver_snapshot_dxdy
               rec_y = rg_mesh_ymax - real(iy-1)*rg_receiver_snapshot_dxdy
   
               call search_closest_quad_gll(rec_x                       &
                                                        ,rec_y                       &
                                                        ,ig_quadf_gll_glonum &
                                                        ,rec_dmin                    &
                                                        ,rec_lgll                    &
                                                        ,rec_mgll                    &
                                                        ,rec_iel                     )
   
               tg_receiver_snapshot_quad(jrec)%x    = rec_x
               tg_receiver_snapshot_quad(jrec)%y    = rec_y
               tg_receiver_snapshot_quad(jrec)%dmin = rec_dmin
               tg_receiver_snapshot_quad(jrec)%iel  = rec_iel 
               tg_receiver_snapshot_quad(jrec)%lgll = rec_lgll
               tg_receiver_snapshot_quad(jrec)%mgll = rec_mgll
               tg_receiver_snapshot_quad(jrec)%rglo = irec
   
            endif
         
         enddo
      enddo
      deallocate(is_rec_in_cpu)
!   
!---->compute receivers local coordinates and initialize PGx
      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(a)') " -->computing local coordinates of receivers"
         call flush(IG_LST_UNIT)
      endif
   
      do irec = 1,ig_nreceiver_snapshot
   
         call compute_info_quad_receiver(tg_receiver_snapshot_quad(irec))
   
         tg_receiver_snapshot_quad(irec)%pgd_x = 0.0
         tg_receiver_snapshot_quad(irec)%pgd_y = 0.0
         tg_receiver_snapshot_quad(irec)%pgd_z = 0.0
   
         tg_receiver_snapshot_quad(irec)%pgv_x = 0.0
         tg_receiver_snapshot_quad(irec)%pgv_y = 0.0
         tg_receiver_snapshot_quad(irec)%pgv_z = 0.0
   
         tg_receiver_snapshot_quad(irec)%pga_x = 0.0
         tg_receiver_snapshot_quad(irec)%pga_y = 0.0
         tg_receiver_snapshot_quad(irec)%pga_z = 0.0

      enddo
!
!---->cpu0 gathers the number of receiver of all cpus. Allocation is done over all cpus to avoid seg_fault with option -check all
      allocate(ig_receiver_snapshot_total_number(ig_ncpu))
      allocate(ig_receiver_snapshot_mpi_shift(ig_ncpu))
      call mpi_gather(ig_nreceiver_snapshot,1,mpi_integer,ig_receiver_snapshot_total_number,1,mpi_integer,0,mpi_comm_world,ios)

      if (ig_myrank == 0) then
!
!------->cpu0 prepares array ig_receiver_snapshot_mpi_shift for mpi_gatherv
         ig_receiver_snapshot_mpi_shift(1) = 0
         do icpu = 2,ig_ncpu
            ig_receiver_snapshot_mpi_shift(icpu) = ig_receiver_snapshot_mpi_shift(icpu-1) + ig_receiver_snapshot_total_number(icpu-1)
         enddo

      endif
! 
!---->cpu0 gathers global number of all snapshot receivers of all cpus (because cpu0 writes the *.grd files alone)
      allocate(ig_receiver_snapshot_glonum(nx*ny),stat=ios)
      if (ios /= 0) then
         write(info,'(a)') "error in subroutine init_snapshot while allocating ig_receiver_snapshot_glonum"
         call error_stop(info)
      endif
! 
!---->cpu0 gathers z-coordinate of all snapshot receivers of all cpus (because cpu0 writes the VTK files alone)
      allocate(rldum3(nx*ny),stat=ios)
      if (ios /= 0) then
         write(info,'(a)') "error in subroutine init_snapshot while allocating rldum3"
         call error_stop(info)
      endif
!
!---->cpu0 gathers global number and z-coordinate of all snapshot receivers of all cpus 
      allocate(ildum1(ig_nreceiver_snapshot),stat=ios)
      allocate(rldum2(ig_nreceiver_snapshot),stat=ios)
      if (ios /= 0) then
         write(info,'(a)') "error in subroutine init_snapshot while allocating *dum*"
         call error_stop(info)
      else
         do irec = 1,ig_nreceiver_snapshot
            ildum1(irec) = tg_receiver_snapshot_quad(irec)%rglo
            rldum2(irec) = tg_receiver_snapshot_quad(irec)%z
         enddo
      endif
   
      call mpi_gatherv(ildum1                             &
                      ,ig_nreceiver_snapshot &
                      ,mpi_integer                        &
                      ,ig_receiver_snapshot_glonum    &
                      ,ig_receiver_snapshot_total_number    &
                      ,ig_receiver_snapshot_mpi_shift         &
                      ,mpi_integer                        &
                      ,0                                  &
                      ,mpi_comm_world                     &
                      ,ios)

      call mpi_gatherv(rldum2                             &
                      ,ig_nreceiver_snapshot &
                      ,mpi_real                           &
                      ,rldum3                             &
                      ,ig_receiver_snapshot_total_number    &
                      ,ig_receiver_snapshot_mpi_shift         &
                      ,mpi_real                           &
                      ,0                                  &
                      ,mpi_comm_world                     &
                      ,ios)

      deallocate(rldum2)

      if (ig_myrank == 0) then

         allocate(rg_receiver_snapshot_z(nx*ny),stat=ios)
         if (ios /= 0) then
            write(info,'(a)') "error in subroutine init_snapshot while allocating rg_receiver_snapshot_z"
            call error_stop(info)
         endif

         allocate(ig_receiver_snapshot_locnum(rg_receiver_snapshot_nx*rg_receiver_snapshot_ny),stat=ios)
         if (ios /= 0) then
            write(info,'(a)') "error in subroutine init_snapshot while allocating ig_receiver_snapshot_locnum"
            call error_stop(info)
         endif

         irec = 0
 
         do iy = 1,rg_receiver_snapshot_ny
            do ix = 1,rg_receiver_snapshot_nx

               irec = irec + 1

               do iz = 1,rg_receiver_snapshot_nx*rg_receiver_snapshot_ny

                  if (ig_receiver_snapshot_glonum(iz) == irec) then

                     ig_receiver_snapshot_locnum  (irec) = iz
                     rg_receiver_snapshot_z(irec) = rldum3(iz)

                  endif

               enddo
         
            enddo
         enddo

      endif
 
      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(a)') "done"
         call flush(IG_LST_UNIT)
      endif
  
      return

!***********************************************************************************************************************************************************************************
   end subroutine init_snapshot
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine computes and writes @f$x,y,z@f$-displacements, velocities and accelerations at receivers used for free surface snapshots.
!>@return Snapshots are written in binary native GMT format and/or in binary VTK xml format.
!***********************************************************************************************************************************************************************************
   subroutine write_snapshot()
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only :&
                                      LG_SNAPSHOT_GMT&
                                     ,LG_SNAPSHOT_VTK&
                                     ,IG_NDOF&
                                     ,IG_NGLL&
                                     ,ig_idt&
                                     ,cg_prefix&
                                     ,lg_snapshot_displacement&
                                     ,lg_snapshot_velocity&
                                     ,lg_snapshot_acceleration&
                                     ,ig_nreceiver_snapshot&
                                     ,ig_snapshot_saving_incr&
                                     ,tg_receiver_snapshot_quad&
                                     ,rg_gll_displacement&
                                     ,rg_gll_velocity&
                                     ,rg_gll_acceleration&
                                     ,ig_myrank

      use mod_gll_value       , only : get_quad_gll_value

      use mod_lagrange        , only : quad_lagrange_interp
     
      implicit none

      real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL)    :: gll_dis
      real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL)    :: gll_vel
      real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL)    :: gll_acc

      real, dimension(ig_nreceiver_snapshot) :: ux
      real, dimension(ig_nreceiver_snapshot) :: uy
      real, dimension(ig_nreceiver_snapshot) :: uz
      real, dimension(ig_nreceiver_snapshot) :: vx
      real, dimension(ig_nreceiver_snapshot) :: vy
      real, dimension(ig_nreceiver_snapshot) :: vz
      real, dimension(ig_nreceiver_snapshot) :: ax
      real, dimension(ig_nreceiver_snapshot) :: ay
      real, dimension(ig_nreceiver_snapshot) :: az

      integer                                             :: irec
                                                          
      character(len=255)                                  :: fname
      character(len=  6)                                  :: csnapshot
      character(len=  4)                                  :: vname
!
!
!*********************************************************************************************************************************
!---->before writing snapshots, displacement, velocity and acceleration are interpolated at the snapshot_receiver x,y,z location
!*********************************************************************************************************************************

      do irec = 1,ig_nreceiver_snapshot

!
!------->get GLL nodes displacement, velocity and acceleration xyz-values from quadrangle element which contains the snapshot receiver
         call get_quad_gll_value(rg_gll_displacement,tg_receiver_snapshot_quad(irec)%gll,gll_dis)
         call get_quad_gll_value(rg_gll_velocity    ,tg_receiver_snapshot_quad(irec)%gll,gll_vel)
         call get_quad_gll_value(rg_gll_acceleration,tg_receiver_snapshot_quad(irec)%gll,gll_acc)

!
!------->interpolate displacement at the snapshot receiver
         call quad_lagrange_interp(gll_dis,tg_receiver_snapshot_quad(irec)%lag,ux(irec),uy(irec),uz(irec))
!
!------->compute PGD in x, y and z directions
         if ( abs(ux(irec)) > tg_receiver_snapshot_quad(irec)%pgd_x ) tg_receiver_snapshot_quad(irec)%pgd_x = abs(ux(irec))
         if ( abs(uy(irec)) > tg_receiver_snapshot_quad(irec)%pgd_y ) tg_receiver_snapshot_quad(irec)%pgd_y = abs(uy(irec))
         if ( abs(uz(irec)) > tg_receiver_snapshot_quad(irec)%pgd_z ) tg_receiver_snapshot_quad(irec)%pgd_z = abs(uz(irec))
!
!------->interpolate velocity at the snapshot receiver
         call quad_lagrange_interp(gll_vel,tg_receiver_snapshot_quad(irec)%lag,vx(irec),vy(irec),vz(irec))
!
!------->compute PGV in x, y and z directions
         if ( abs(vx(irec)) > tg_receiver_snapshot_quad(irec)%pgv_x ) tg_receiver_snapshot_quad(irec)%pgv_x = abs(vx(irec))
         if ( abs(vy(irec)) > tg_receiver_snapshot_quad(irec)%pgv_y ) tg_receiver_snapshot_quad(irec)%pgv_y = abs(vy(irec))
         if ( abs(vz(irec)) > tg_receiver_snapshot_quad(irec)%pgv_z ) tg_receiver_snapshot_quad(irec)%pgv_z = abs(vz(irec))
         if ( sqrt(vx(irec)**2 + vy(irec)**2 + vz(irec)**2)  > tg_receiver_snapshot_quad(irec)%pgv_xyz ) tg_receiver_snapshot_quad(irec)%pgv_xyz = sqrt(vx(irec)**2 + vy(irec)**2 + vz(irec)**2)
!
!------->interpolate acceleration at the snapshot receiver
         call quad_lagrange_interp(gll_acc,tg_receiver_snapshot_quad(irec)%lag,ax(irec),ay(irec),az(irec))
!
!------->compute PGA in x, y and z directions
         if ( abs(ax(irec)) > tg_receiver_snapshot_quad(irec)%pga_x ) tg_receiver_snapshot_quad(irec)%pga_x = abs(ax(irec))
         if ( abs(ay(irec)) > tg_receiver_snapshot_quad(irec)%pga_y ) tg_receiver_snapshot_quad(irec)%pga_y = abs(ay(irec))
         if ( abs(az(irec)) > tg_receiver_snapshot_quad(irec)%pga_z ) tg_receiver_snapshot_quad(irec)%pga_z = abs(az(irec))

      enddo
!
!
!*********************************************************************************************************************************
!---->write snapshot in gmt or VTK format
!*********************************************************************************************************************************
      write(csnapshot,'(i6.6)') ig_idt
!
!
!---->snapshot displacement
      if ( lg_snapshot_displacement .and. (mod(ig_idt-1,ig_snapshot_saving_incr) == 0) ) then

         if (LG_SNAPSHOT_GMT) then

            fname = trim(cg_prefix)//".snapshot.ux."//trim(csnapshot)//".grd"
            call write_snapshot_gmt_native_format(fname,ux)

            fname = trim(cg_prefix)//".snapshot.uy."//trim(csnapshot)//".grd"
            call write_snapshot_gmt_native_format(fname,uy)

            fname = trim(cg_prefix)//".snapshot.uz."//trim(csnapshot)//".grd"
            call write_snapshot_gmt_native_format(fname,uz)

         endif

         if (LG_SNAPSHOT_VTK) then

            fname = trim(cg_prefix)//".snapshot.uxyz."//trim(csnapshot)//".vts"
            vname = "disp"
            call write_snapshot_vtk(fname,vname,ux,uy,uz)

         endif

      endif
!
!
!---->snapshot velocity
      if ( lg_snapshot_velocity .and. (mod(ig_idt-1,ig_snapshot_saving_incr) == 0) ) then

         if (LG_SNAPSHOT_GMT) then

            fname = trim(cg_prefix)//".snapshot.vx."//trim(csnapshot)//".grd"
            call write_snapshot_gmt_native_format(fname,vx)
          
            fname = trim(cg_prefix)//".snapshot.vy."//trim(csnapshot)//".grd"
            call write_snapshot_gmt_native_format(fname,vy)
          
            fname = trim(cg_prefix)//".snapshot.vz."//trim(csnapshot)//".grd"
            call write_snapshot_gmt_native_format(fname,vz)

         endif

         if (LG_SNAPSHOT_VTK) then

            fname = trim(cg_prefix)//".snapshot.vxyz."//trim(csnapshot)//".vts"
            vname = "velo"
            call write_snapshot_vtk(fname,vname,vx,vy,vz)

         endif

      endif
!
!
!---->snapshot acceleration
      if ( lg_snapshot_acceleration .and. (mod(ig_idt-1,ig_snapshot_saving_incr) == 0) )  then

         if (LG_SNAPSHOT_GMT) then

            fname = trim(cg_prefix)//".snapshot.ax."//trim(csnapshot)//".grd"
            call write_snapshot_gmt_native_format(fname,ax)
          
            fname = trim(cg_prefix)//".snapshot.ay."//trim(csnapshot)//".grd"
            call write_snapshot_gmt_native_format(fname,ay)
          
            fname = trim(cg_prefix)//".snapshot.az."//trim(csnapshot)//".grd"
            call write_snapshot_gmt_native_format(fname,az)

         endif

         if (LG_SNAPSHOT_VTK) then

            fname = trim(cg_prefix)//".snapshot.axyz."//trim(csnapshot)//".vts"
            vname = "acce"
            call write_snapshot_vtk(fname,vname,ax,ay,az)

         endif

      endif

      return
     
!***********************************************************************************************************************************************************************************
   end subroutine write_snapshot
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine writes structured grid snapshot in binary native GMT format
!>@param fname : file name of the snapshot
!>@param val   : values to be written
!***********************************************************************************************************************************************************************************
   subroutine write_snapshot_gmt_native_format(fname,val)
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only :&
                                      ig_nreceiver_snapshot&
                                     ,tg_receiver_snapshot_quad&
                                     ,ig_receiver_snapshot_glonum&
                                     ,ig_receiver_snapshot_total_number&
                                     ,ig_receiver_snapshot_mpi_shift&
                                     ,rg_receiver_snapshot_dxdy&
                                     ,rg_mesh_xmax&
                                     ,rg_mesh_ymax&
                                     ,rg_mesh_xmin&
                                     ,rg_mesh_ymin&
                                     ,rg_receiver_snapshot_nx&
                                     ,rg_receiver_snapshot_ny&
                                     ,ig_myrank&
                                     ,ig_ncpu

      implicit none

      character(len=255), intent(in)      :: fname
      real              , intent(in)      :: val(ig_nreceiver_snapshot)

      integer           , parameter       :: HEADER_SIZE = 892

      real                                :: val_max
      real                                :: val_min
      real                                :: val_max_all_cpu
      real                                :: val_min_all_cpu
      real, allocatable                   :: val_all_cpu(:)

      integer                             :: nboctet_real
      integer                             :: irec
      integer                             :: rglo
      integer                             :: desc
      integer                             :: ios
      integer(kind=mpi_offset_kind)       :: pos
      integer, dimension(mpi_status_size) :: statut

!
!---->compute min and max. cpu 0 gather data for writing gmt header
      if (ig_nreceiver_snapshot > 0) then
         val_max = maxval(val)
         val_min = minval(val)
      else
         val_min = 0.0
         val_max = 0.0
      endif
      call mpi_reduce(val_max,val_max_all_cpu,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ios)
      call mpi_reduce(val_min,val_min_all_cpu,1,MPI_REAL,MPI_MIN,0,MPI_COMM_WORLD,ios)
!
!---->open file 'fname' for snapshot
      call mpi_file_open(mpi_comm_world,fname,mpi_mode_wronly + mpi_mode_create,mpi_info_null,desc,ios)
      call mpi_type_size(mpi_real,nboctet_real,ios)

      if (ig_myrank == 0) then
!
!------>cpu 0 writes gmt grd file header
        call write_header_gmt_native_format(desc                                &
                                           ,rg_receiver_snapshot_nx   ,rg_receiver_snapshot_ny    &
                                           ,rg_mesh_xmin   ,rg_mesh_xmax    &
                                           ,rg_mesh_ymin   ,rg_mesh_ymax    &
                                           ,val_min_all_cpu  ,val_max_all_cpu   &
                                           ,rg_receiver_snapshot_dxdy,rg_receiver_snapshot_dxdy )

      endif
!
!---->cpu0 allocates array val_all_cpu to gather val of all cpus
!     allocation is done over all cpus to avoid seg_fault with -check all compiler option
      allocate(val_all_cpu(rg_receiver_snapshot_nx*rg_receiver_snapshot_ny))

!
!---->cpu0 gathers array val to write gmt *.grd file
      call mpi_gatherv(val                                &
                      ,ig_nreceiver_snapshot &
                      ,mpi_real                           &
                      ,val_all_cpu                        &
                      ,ig_receiver_snapshot_total_number    &
                      ,ig_receiver_snapshot_mpi_shift         &
                      ,mpi_real                           &
                      ,0                                  &
                      ,mpi_comm_world                     &
                      ,ios                                )
!
!---->cpu0 writes all receivers of all cpus using global numbering of snapshot receiver
      if (ig_myrank == 0) then
         do irec = 1,rg_receiver_snapshot_nx*rg_receiver_snapshot_ny
 
            rglo =  ig_receiver_snapshot_glonum(irec)

            pos  = HEADER_SIZE + (rglo-1)*nboctet_real
       
            call mpi_file_write_at(desc,pos,val_all_cpu(irec),1,mpi_real,statut,ios)
       
         enddo
      endif
!
!---->close files
      call mpi_file_close(desc,ios)

      return

!***********************************************************************************************************************************************************************************
   end subroutine write_snapshot_gmt_native_format
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine writes header of binary native GMT files.
!>@param desc  : mpi_unit of the file
!>@param nx    : number of receivers along @f$x@f$-direction
!>@param ny    : number of receivers along @f$y@f$-direction
!>@param x_min : minimum @f$x@f$-coordinate
!>@param x_max : maximum @f$x@f$-coordinate
!>@param y_min : minimum @f$y@f$-coordinate
!>@param y_max : maximum @f$y@f$-coordinate
!>@param z_min : minimum @f$z@f$-coordinate
!>@param z_max : maximum @f$z@f$-coordinate
!>@param x_inc : @f$x@f$-increment between receivers
!>@param y_inc : @f$y@f$-increment between receivers
!***********************************************************************************************************************************************************************************
   subroutine write_header_gmt_native_format(desc,nx,ny,x_min,x_max,y_min,y_max,z_min,z_max,x_inc,y_inc)
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only : ig_myrank

      implicit none

      integer, intent(in)                 :: desc
      integer, intent(in)                 :: nx
      integer, intent(in)                 :: ny
      real   , intent(in)                 :: x_min
      real   , intent(in)                 :: x_max
      real   , intent(in)                 :: y_min
      real   , intent(in)                 :: y_max
      real   , intent(in)                 :: z_min
      real   , intent(in)                 :: z_max
      real   , intent(in)                 :: x_inc
      real   , intent(in)                 :: y_inc

      integer, parameter                  :: GRD_COMMAND_LEN = 320
      integer, parameter                  :: GRD_REMARK_LEN  = 160
      integer, parameter                  :: GRD_TITLE_LEN   =  80
      integer, parameter                  :: GRD_UNIT_LEN    =  80
      integer, parameter                  :: GRD_VARNAME_LEN =  80
      integer, parameter                  :: NODE_OFFSET     =   0
      real   , parameter                  :: Z_SCALE_FACTOR  = 1.0 !grd values must be multiplied by this
      real   , parameter                  :: Z_ADD_OFFSET    = 0.0 !After scaling, add this

      integer                             :: nboctet_double_precision
      integer                             :: nboctet_char 
      integer                             :: nboctet_int
      integer                             :: ios
      integer(kind=mpi_offset_kind)       :: pos
      integer, dimension(mpi_status_size) :: statut

      character(len=GRD_UNIT_LEN)         :: x_units !units in x-direction
      character(len=GRD_UNIT_LEN)         :: y_units !units in y-direction
      character(len=GRD_UNIT_LEN)         :: z_units !grid value units
      character(len=GRD_TITLE_LEN)        :: title   !name of data set
      character(len=GRD_COMMAND_LEN)      :: command !name of generating command
      character(len=GRD_REMARK_LEN)       :: remark  !comments on this data set


      call mpi_type_size(mpi_double_precision,nboctet_double_precision,ios)
      call mpi_type_size(mpi_character       ,nboctet_char            ,ios)
      call mpi_type_size(mpi_integer         ,nboctet_int             ,ios)

      x_units = "m"//CHAR(0)
      y_units = "m"//CHAR(0)
      z_units = " "//CHAR(0)
      title   = "EFISPEC3D snapshot"//CHAR(0)
      command = "mpi_file_write_at()"//CHAR(0)
      remark  = "create by EFISPEC3D: http://efispec.free.fr"//CHAR(0)

      pos = 0
      call mpi_file_write_at(desc,pos,nx,1,mpi_integer,statut,ios)

      pos = pos + nboctet_int
      call mpi_file_write_at(desc,pos,ny,1,mpi_integer,statut,ios)

      pos = pos + nboctet_int
      call mpi_file_write_at(desc,pos,NODE_OFFSET,1,mpi_integer,statut,ios)

      pos = pos + nboctet_int
      call mpi_file_write_at(desc,pos,dble(x_min),1,mpi_double_precision,statut,ios)

      pos = pos + nboctet_double_precision
      call mpi_file_write_at(desc,pos,dble(x_max),1,mpi_double_precision,statut,ios)

      pos = pos + nboctet_double_precision
      call mpi_file_write_at(desc,pos,dble(y_min),1,mpi_double_precision,statut,ios)

      pos = pos + nboctet_double_precision
      call mpi_file_write_at(desc,pos,dble(y_max),1,mpi_double_precision,statut,ios)

      pos = pos + nboctet_double_precision
      call mpi_file_write_at(desc,pos,dble(z_min),1,mpi_double_precision,statut,ios)

      pos = pos + nboctet_double_precision
      call mpi_file_write_at(desc,pos,dble(z_max),1,mpi_double_precision,statut,ios)

      pos = pos + nboctet_double_precision
      call mpi_file_write_at(desc,pos,dble(x_inc),1,mpi_double_precision,statut,ios)

      pos = pos + nboctet_double_precision
      call mpi_file_write_at(desc,pos,dble(y_inc),1,mpi_double_precision,statut,ios)

      pos = pos + nboctet_double_precision
      call mpi_file_write_at(desc,pos,dble(Z_SCALE_FACTOR),1,mpi_double_precision,statut,ios)

      pos = pos + nboctet_double_precision
      call mpi_file_write_at(desc,pos,dble(Z_ADD_OFFSET),1,mpi_double_precision,statut,ios)

      pos = pos + nboctet_double_precision
      call mpi_file_write_at(desc,pos,x_units,GRD_UNIT_LEN,mpi_character,statut,ios)

      pos = pos + nboctet_char*GRD_UNIT_LEN
      call mpi_file_write_at(desc,pos,y_units,GRD_UNIT_LEN,mpi_character,statut,ios)

      pos = pos + nboctet_char*GRD_UNIT_LEN
      call mpi_file_write_at(desc,pos,z_units,GRD_UNIT_LEN,mpi_character,statut,ios)

      pos = pos + nboctet_char*GRD_UNIT_LEN
      call mpi_file_write_at(desc,pos,title,GRD_TITLE_LEN,mpi_character,statut,ios)

      pos = pos + nboctet_char*GRD_TITLE_LEN
      call mpi_file_write_at(desc,pos,command,GRD_COMMAND_LEN,mpi_character,statut,ios)

      pos = pos + nboctet_char*GRD_COMMAND_LEN
      call mpi_file_write_at(desc,pos,remark,GRD_REMARK_LEN,mpi_character,statut,ios)

      return

!***********************************************************************************************************************************************************************************
   end subroutine write_header_gmt_native_format
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine writes structured grid snapshot in binary VTK xml format.
!>@param fname   : file name of the snapshot
!>@param varname : variable name
!>@param valx    : @f$x@f$-direction values
!>@param valy    : @f$y@f$-direction values
!>@param valz    : @f$z@f$-direction values
!***********************************************************************************************************************************************************************************
   subroutine write_snapshot_vtk(fname,varname,valx,valy,valz)
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only :&
                                      ig_nreceiver_snapshot&
                                     ,ig_receiver_snapshot_total_number&
                                     ,ig_receiver_snapshot_mpi_shift&
                                     ,ig_receiver_snapshot_locnum&
                                     ,rg_receiver_snapshot_z&
                                     ,ig_myrank&
                                     ,error_stop&
                                     ,rg_receiver_snapshot_nx&
                                     ,rg_receiver_snapshot_ny&
                                     ,rg_receiver_snapshot_dxdy&
                                     ,rg_mesh_xmin&
                                     ,rg_mesh_ymax

      use mod_vtk_io

      implicit none

      character(len=255), intent(in) :: fname
      character(len=  4), intent(in) :: varname

      real              , intent(in) :: valx(ig_nreceiver_snapshot)
      real              , intent(in) :: valy(ig_nreceiver_snapshot)
      real              , intent(in) :: valz(ig_nreceiver_snapshot)
                        
      real                           :: val_all_cpu(rg_receiver_snapshot_nx*rg_receiver_snapshot_ny)
      real                           :: val_all_cpu_order(rg_receiver_snapshot_nx*rg_receiver_snapshot_ny)
      real                           :: xrec(rg_receiver_snapshot_nx*rg_receiver_snapshot_ny)
      real                           :: yrec(rg_receiver_snapshot_nx*rg_receiver_snapshot_ny)
                        
      integer                        :: ix
      integer                        :: iy
      integer                        :: irec
      integer                        :: irec_gatherv
      integer                        :: ios

      character(len=255)             :: info

!
!---->cpu0 gathers array valx to be store in VTK file
      call mpi_gatherv(valx                               &
                      ,ig_nreceiver_snapshot &
                      ,mpi_real                           &
                      ,val_all_cpu                        &
                      ,ig_receiver_snapshot_total_number    &
                      ,ig_receiver_snapshot_mpi_shift         &
                      ,mpi_real                           &
                      ,0                                  &
                      ,mpi_comm_world                     &
                      ,ios                                )
!
!---->cpu0 writes all receivers of all cpus using global numbering of snapshot receiver
      if (ig_myrank == 0) then

!
!------->compute x, y and z coordinates of receiver snapshot 
         irec = 0

         do iy = 1,rg_receiver_snapshot_ny
            do ix = 1,rg_receiver_snapshot_nx

               irec                    = irec + 1
               irec_gatherv            = ig_receiver_snapshot_locnum(irec)
                                       
               xrec(irec)              = rg_mesh_xmin + real(ix-1)*rg_receiver_snapshot_dxdy
               yrec(irec)              = rg_mesh_ymax - real(iy-1)*rg_receiver_snapshot_dxdy !y for VTK coordinate system

               val_all_cpu_order(irec) = val_all_cpu(irec_gatherv)

            enddo
         enddo 

!
!------->initialize VTK file
         ios = VTK_INI_XML(                                 &
                           output_format = 'BINARY'         &
                          ,filename      = trim(fname)      &
                          ,mesh_topology = 'StructuredGrid' &
                          ,nx1           = 1                &
                          ,nx2           = rg_receiver_snapshot_nx   &
                          ,ny1           = 1                &
                          ,ny2           = rg_receiver_snapshot_ny   &
                          ,nz1           = 1                &
                          ,nz2           = 1                )

         if (ios /= 0) then
            write(info,'(a)') "error in subroutine write_snapshot_vtk: VTK_INI_XML"
            call error_stop(info)
         endif

!
!------->store free surface geometry
         ios = VTK_GEO_XML(                                              &
                           nx1           = 1                             &
                          ,nx2           = rg_receiver_snapshot_nx                &
                          ,ny1           = 1                             &
                          ,ny2           = rg_receiver_snapshot_ny                &
                          ,nz1           = 1                             &
                          ,nz2           = 1                             &
                          ,NN            = rg_receiver_snapshot_nx*rg_receiver_snapshot_ny &
                          ,X             = xrec                          &
                          ,Y             = yrec                          &
                          ,Z             = rg_receiver_snapshot_z  )

         if (ios /= 0) then
            write(info,'(a)') "error in subroutine write_snapshot_vtk: VTK_GEO_XML (init)"
            call error_stop(info)
         endif
         
!
!------->store x-component
         ios = VTK_DAT_XML(                                              &
                           var_location     = 'NODE'                     &
                          ,var_block_action = 'OPEN'                     )

         if (ios /= 0) then
            write(info,'(a)') "error in subroutine write_snapshot_vtk: VTK_DAT_XML (open)"
            call error_stop(info)
         endif
         
         ios = VTK_VAR_XML(                                        &
                           NC_NN   = rg_receiver_snapshot_nx*rg_receiver_snapshot_ny &
                          ,varname = "X-"//trim(varname)           &
                          ,var     = val_all_cpu_order             &
                                                                   )

         if (ios /= 0) then
            write(info,'(a)') "error in subroutine write_snapshot_vtk: VTK_VAR_XML X-component"
            call error_stop(info)
         endif

      endif

!
!---->cpu0 gathers array valy to be store in VTK file
      call mpi_gatherv(valy                               &
                      ,ig_nreceiver_snapshot &
                      ,mpi_real                           &
                      ,val_all_cpu                        &
                      ,ig_receiver_snapshot_total_number    &
                      ,ig_receiver_snapshot_mpi_shift         &
                      ,mpi_real                           &
                      ,0                                  &
                      ,mpi_comm_world                     &
                      ,ios                                )


      if (ig_myrank == 0) then

!
!------->order data from gatherv
         irec = 0
         do iy = 1,rg_receiver_snapshot_ny
            do ix = 1,rg_receiver_snapshot_nx

               irec                    = irec + 1
               irec_gatherv            = ig_receiver_snapshot_locnum(irec)
               val_all_cpu_order(irec) = val_all_cpu(irec_gatherv)

            enddo
         enddo 

!
!------->store y-component
         ios = VTK_VAR_XML(                                        &
                           NC_NN   = rg_receiver_snapshot_nx*rg_receiver_snapshot_ny &
                          ,varname = "Y-"//trim(varname)           &
                          ,var     = val_all_cpu_order             &
                                                                   )
         if (ios /= 0) then
            write(info,'(a)') "error in subroutine write_snapshot_vtk: VTK_VAR_XML Y-component"
            call error_stop(info)
         endif

      endif

!
!---->cpu0 gathers array valz to be store in VTK file
      call mpi_gatherv(valz                               &
                      ,ig_nreceiver_snapshot &
                      ,mpi_real                           &
                      ,val_all_cpu                        &
                      ,ig_receiver_snapshot_total_number    &
                      ,ig_receiver_snapshot_mpi_shift         &
                      ,mpi_real                           &
                      ,0                                  &
                      ,mpi_comm_world                     &
                      ,ios                                )

      if (ig_myrank == 0) then

!
!------->order data from gatherv
         irec = 0
         do iy = 1,rg_receiver_snapshot_ny
            do ix = 1,rg_receiver_snapshot_nx

               irec                    = irec + 1
               irec_gatherv            = ig_receiver_snapshot_locnum(irec)
               val_all_cpu_order(irec) = val_all_cpu(irec_gatherv)

            enddo
         enddo
!
!------->store z-component
         ios = VTK_VAR_XML(                                        &
                           NC_NN   = rg_receiver_snapshot_nx*rg_receiver_snapshot_ny &
                          ,varname = "Z-"//trim(varname)           &
                          ,var     = val_all_cpu_order             &
                                                                   )

         if (ios /= 0) then
            write(info,'(a)') "error in subroutine write_snapshot_vtk: VTK_VAR_XML Z-component"
            call error_stop(info)
         endif
         
!
!------->finalize VTK file
         ios = VTK_DAT_XML(                                              &
                           var_location     = 'NODE'                     &
                          ,var_block_action = 'CLOSE'                    )

         if (ios /= 0) then
            write(info,'(a)') "error in subroutine write_snapshot_vtk: VTK_DAT_XML (close)"
            call error_stop(info)
         endif
         
         ios = VTK_GEO_XML()

         if (ios /= 0) then
            write(info,'(a)') "error in subroutine write_snapshot_vtk: VTK_GEO_XML (finalize)"
            call error_stop(info)
         endif

         ios = VTK_END_XML()

         if (ios /= 0) then
            write(info,'(a)') "error in subroutine write_snapshot_vtk: VTK_END_XML"
            call error_stop(info)
         endif

      endif

      return
!***********************************************************************************************************************************************************************************
   end subroutine write_snapshot_vtk
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine selects which ParaView collection files should be written (displacement, velocity and/or acceleration collections)
!!depending on value of variables mod_global_variables::lg_snapshot_displacement, mod_global_variables::lg_snapshot_velocity and mod_global_variables::lg_snapshot_acceleration.
!***********************************************************************************************************************************************************************************
   subroutine write_collection_vtk()
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only :&
                                      lg_snapshot_displacement&
                                     ,lg_snapshot_velocity&
                                     ,lg_snapshot_acceleration
               
      implicit none

      if (lg_snapshot_displacement) then
         call collection_vtk("snapshot.uxyz")
      endif

      if (lg_snapshot_velocity) then
         call collection_vtk("snapshot.vxyz")
      endif

      if (lg_snapshot_acceleration) then
         call collection_vtk("snapshot.axyz")
      endif

      return
!***********************************************************************************************************************************************************************************
   end subroutine write_collection_vtk
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine write Paraview collection file *.pvd of VTK xml snapshot files *.vts.
!>@param varname : variable name
!***********************************************************************************************************************************************************************************
   subroutine collection_vtk(varname)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only :&
                                      get_newunit&
                                     ,cg_prefix&
                                     ,rg_dt&
                                     ,ig_ndt&
                                     ,ig_snapshot_saving_incr
               
      implicit none

      character(len=*), intent(in) :: varname

      real                         :: time

      integer                      :: istep
      integer                      :: myunit

      character(len=255)           :: fname
      character(len=6  )           :: csnapshot

      open(unit=get_newunit(myunit),file=trim(cg_prefix)//".collection."//trim(varname)//".pvd")

      write(unit=myunit,fmt='(a)') "<?xml version=""1.0""?>"
      write(unit=myunit,fmt='(a)') "<VTKFile type=""Collection"" version=""0.1"" byte_order=""BigEndian"">"
      write(unit=myunit,fmt='(a)') "  <Collection>"

      do istep = 1,ig_ndt,ig_snapshot_saving_incr

         write(csnapshot,'(I6.6)') istep

         time  = real(istep-1)*rg_dt

         fname = trim(cg_prefix)//"."//trim(varname)//"."//trim(csnapshot)//".vts"

         write(unit=myunit,fmt='(a,E14.7,3a)') "    <DataSet timestep=""",time,""" group="""" part=""0"" file=""",trim(fname),"""/>"

      enddo

      write(unit=myunit,fmt='(a)') "  </Collection>"
      write(unit=myunit,fmt='(a)') "</VTKFile>"

      close(myunit)

      return
!***********************************************************************************************************************************************************************************
   end subroutine collection_vtk
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine writes peak ground motions files in GMT or VTK formats.
!***********************************************************************************************************************************************************************************
   subroutine write_peak_ground_motion()
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only :&
                                      tg_receiver_snapshot_quad&
                                     ,ig_myrank&
                                     ,ig_nreceiver_snapshot&
                                     ,cg_prefix&
                                     ,LG_SNAPSHOT_GMT&
                                     ,LG_SNAPSHOT_VTK&
                                     ,IG_LST_UNIT

      implicit none

      real, dimension(ig_nreceiver_snapshot) :: pgd_x
      real, dimension(ig_nreceiver_snapshot) :: pgd_y
      real, dimension(ig_nreceiver_snapshot) :: pgd_z
      real, dimension(ig_nreceiver_snapshot) :: pgv_x
      real, dimension(ig_nreceiver_snapshot) :: pgv_y
      real, dimension(ig_nreceiver_snapshot) :: pgv_z
      real, dimension(ig_nreceiver_snapshot) :: pga_x
      real, dimension(ig_nreceiver_snapshot) :: pga_y
      real, dimension(ig_nreceiver_snapshot) :: pga_z

      integer                                             :: irec

      character(len=255)                                  :: fname

      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(" ",/,a)') "writing snapshots of peak ground motion..."
         call flush(IG_LST_UNIT)
      endif

      do irec = 1,ig_nreceiver_snapshot
         pgd_x(irec) = tg_receiver_snapshot_quad(irec)%pgd_x
      enddo

      do irec = 1,ig_nreceiver_snapshot
         pgd_y(irec) = tg_receiver_snapshot_quad(irec)%pgd_y
      enddo

      do irec = 1,ig_nreceiver_snapshot
         pgd_z(irec) = tg_receiver_snapshot_quad(irec)%pgd_z
      enddo

      do irec = 1,ig_nreceiver_snapshot
         pgv_x(irec) = tg_receiver_snapshot_quad(irec)%pgv_x
      enddo

      do irec = 1,ig_nreceiver_snapshot
         pgv_y(irec) = tg_receiver_snapshot_quad(irec)%pgv_y
      enddo

      do irec = 1,ig_nreceiver_snapshot
         pgv_z(irec) = tg_receiver_snapshot_quad(irec)%pgv_z
      enddo

      do irec = 1,ig_nreceiver_snapshot
         pga_x(irec) = tg_receiver_snapshot_quad(irec)%pga_x
      enddo

      do irec = 1,ig_nreceiver_snapshot
         pga_y(irec) = tg_receiver_snapshot_quad(irec)%pga_y
      enddo

      do irec = 1,ig_nreceiver_snapshot
         pga_z(irec) = tg_receiver_snapshot_quad(irec)%pga_z
      enddo

      if (LG_SNAPSHOT_GMT) then

         fname = trim(cg_prefix)//".snapshot.PGDx.grd"
         call write_snapshot_gmt_native_format(fname,pgd_x)
         
         fname = trim(cg_prefix)//".snapshot.PGDy.grd"
         call write_snapshot_gmt_native_format(fname,pgd_y)
         
         fname = trim(cg_prefix)//".snapshot.PGDz.grd"
         call write_snapshot_gmt_native_format(fname,pgd_z)
         
         fname = trim(cg_prefix)//".snapshot.PGVx.grd"
         call write_snapshot_gmt_native_format(fname,pgv_x)
         
         fname = trim(cg_prefix)//".snapshot.PGVy.grd"
         call write_snapshot_gmt_native_format(fname,pgv_y)
         
         fname = trim(cg_prefix)//".snapshot.PGVz.grd"
         call write_snapshot_gmt_native_format(fname,pgv_z)
         
         fname = trim(cg_prefix)//".snapshot.PGAx.grd"
         call write_snapshot_gmt_native_format(fname,pga_x)
         
         fname = trim(cg_prefix)//".snapshot.PGAy.grd"
         call write_snapshot_gmt_native_format(fname,pga_y)
         
         fname = trim(cg_prefix)//".snapshot.PGAz.grd"
         call write_snapshot_gmt_native_format(fname,pga_z)

      endif

      if (LG_SNAPSHOT_VTK) then

         fname = trim(cg_prefix)//".snapshot.PGDxyz.vts"
         call write_snapshot_vtk(fname,"displacement",pgd_x,pgd_y,pgd_z)
       
         fname = trim(cg_prefix)//".snapshot.PGVxyz.vts"
         call write_snapshot_vtk(fname,"velocity"    ,pgv_x,pgv_y,pgv_z)
       
         fname = trim(cg_prefix)//".snapshot.PGAxyz.vts"
         call write_snapshot_vtk(fname,"acceleration",pga_x,pga_y,pga_z)

      endif

      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(a)') "done"
         call flush(IG_LST_UNIT)
      endif

      return

!***********************************************************************************************************************************************************************************
   end subroutine write_peak_ground_motion
!***********************************************************************************************************************************************************************************

end module mod_snapshot
