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
!!<b>This file contains the main spectral finite element program: EFISPEC3D</b>.
!!
!>@brief
!!This program is composed by :
!!- an initialization phase;
!!- a time loop;
!!- a finalization phase.

!>@return void
program EFISPEC3D

   use mpi
   
   use mod_global_variables, only :&
                                    ig_ndt&
                                   ,rg_simu_current_time&
                                   ,ig_idt&
                                   ,rg_dt&
                                   ,IG_LST_UNIT&
                                   ,LG_OUTPUT_CPUTIME&
                                   ,LG_SNAPSHOT_VTK&
                                   ,IG_LAGRANGE_ORDER&
                                   ,lg_snapshot&
                                   ,ig_nreceiver_snapshot&
                                   ,ig_receiver_saving_incr&
                                   ,ig_ncpu&
                                   ,ig_nhexa&
                                   ,ig_myrank&
                                   ,cg_myrank&
                                   ,cg_cpu_name&
                                   ,ig_cpu_name_len&
                                   ,cg_ncpu&
                                   ,cg_prefix&
                                   ,get_newunit&
                                   ,error_stop&
                                   ,lg_boundary_absorption&
                                   ,ig_medium_type&
                                   ,ig_nhexa_outer&
                                   ,LG_ASYNC_MPI_COMM

   use mod_init_efi
   
   use mod_init_mesh
   
   use mod_init_mpibuffer
   
   use mod_solver
   
   use mod_init_medium

   use mod_snapshot          , only :&
                                     init_snapshot&
                                    ,write_snapshot&
                                    ,write_peak_ground_motion&
                                    ,write_collection_vtk

   use mod_init_memory

   use mod_mpi_assemble

   use mod_receiver          , only :&
                                     write_receiver_output&
                                    ,init_hexa_receiver&
                                    ,init_quad_receiver

   use mod_source            , only :&
                                     init_double_couple_source&
                                    ,init_single_force_source&
                                    ,compute_double_couple_source
   
   implicit none
   
   double precision, allocatable, dimension(:) :: buffer_double
   double precision, allocatable, dimension(:) :: buffer_double_all_cpu

   double precision                            :: start_time
   double precision                            :: end_time
                                               
   double precision                            :: dltim1
   double precision                            :: dltim2
                                               
   double precision                            :: time_init_mesh
   double precision                            :: time_init_gll
   double precision                            :: time_init_mpi_buffers
   double precision                            :: time_init_medium
   double precision                            :: time_init_source
   double precision                            :: time_init_receiver
   double precision                            :: time_init_jacobian
   double precision                            :: time_init_mass_matrix
   double precision                            :: time_compute_dcsource
   double precision                            :: time_init_snapshot
   double precision                            :: time_memory_consumption
                                               
   double precision                            :: t_tloop
   double precision                            :: t_newmark
   double precision                            :: t_forext
   double precision                            :: t_forabs
   double precision                            :: t_forint_inner
   double precision                            :: t_forint_outer
   double precision                            :: t_resol
   double precision                            :: t_snafsu
   double precision                            :: t_outavd
   double precision                            :: t_linkfo
   double precision                            :: t_link_init
                                               
   integer, parameter                          :: NTIME_INIT=11 
   integer, allocatable, dimension(:)          :: buffer_integer
   integer                                     :: ios
   integer                                     :: icpu
   integer                                     :: itime
   integer                                     :: unit_time
                                               
   character(len= 95)                          :: myfmt
   character(len=255)                          :: info

!
!
!*************************************************************************************************************************
!*************************************************************************************************************************
!phase 1: initialization
!*************************************************************************************************************************
!*************************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   call init_mpi()
   start_time = mpi_wtime()
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   call init_input_variables()
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   time_init_mesh = mpi_wtime()

   call init_mesh()

   time_init_mesh = mpi_wtime() - time_init_mesh
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   time_init_gll = mpi_wtime()

   call init_gll_nodes()
   call init_gll_nodes_coordinates()

   time_init_gll = mpi_wtime() - time_init_gll
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   time_init_mpi_buffers = mpi_wtime()

   call init_mpi_buffers()

   time_init_mpi_buffers = mpi_wtime() - time_init_mpi_buffers
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   time_init_medium = mpi_wtime()

   if (ig_medium_type == 0) then

      if (ig_myrank == 0) write(IG_LST_UNIT,'(" ",/,a)') "generating medium with constant mechanical properties per hexa"
      call init_hexa_medium()

   else

      write(info,'(a)') "error in efispec: illegal medium type"
      call error_stop(info)

   endif

   call init_quadp_medium()

   time_init_medium = mpi_wtime() - time_init_medium
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   time_init_source = mpi_wtime()

   call init_double_couple_source()
   call init_single_force_source()
   
   time_init_source = mpi_wtime() - time_init_source
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   time_init_receiver = mpi_wtime()

   call init_hexa_receiver()
   call init_quad_receiver()

   time_init_receiver = mpi_wtime() - time_init_receiver
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   time_init_jacobian = mpi_wtime()

   call init_jacobian_matrix_hexa()
   call init_jacobian_matrix_quad()
   
   time_init_jacobian = mpi_wtime() - time_init_jacobian
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   time_init_mass_matrix = mpi_wtime()

   call init_mass_matrix()
   
   time_init_mass_matrix = mpi_wtime() - time_init_mass_matrix
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   time_compute_dcsource = mpi_wtime()
   
   call compute_double_couple_source()
   
   time_compute_dcsource = mpi_wtime() - time_compute_dcsource
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   time_init_snapshot = mpi_wtime()
   
   if (lg_snapshot) call init_snapshot()
   
   time_init_snapshot = mpi_wtime() - time_init_snapshot
   !**********************************************************************************************************************


   !
   !
   !**********************************************************************************************************************
   time_memory_consumption = mpi_wtime()

   call memory_consumption()

   time_memory_consumption = mpi_wtime() - time_memory_consumption
   !**********************************************************************************************************************


   end_time = (mpi_wtime() - start_time)
   call mpi_reduce(end_time,dltim1,1,mpi_double_precision,mpi_max,0,mpi_comm_world,ios)


   !
   !
   !**********************************************************************************************************************
   !summarize elapsed time in subroutines for initialization
   !**********************************************************************************************************************
   allocate(buffer_double(NTIME_INIT))

   allocate(buffer_double_all_cpu(NTIME_INIT*ig_ncpu))

   buffer_double( 1) = time_init_mesh
   buffer_double( 2) = time_init_gll
   buffer_double( 3) = time_init_mpi_buffers
   buffer_double( 4) = time_init_medium
   buffer_double( 5) = time_init_source
   buffer_double( 6) = time_init_receiver
   buffer_double( 7) = time_init_jacobian
   buffer_double( 8) = time_init_mass_matrix
   buffer_double( 9) = time_compute_dcsource
   buffer_double(10) = time_init_snapshot
   buffer_double(11) = time_memory_consumption

   call mpi_gather(buffer_double,NTIME_INIT,mpi_double_precision,buffer_double_all_cpu,NTIME_INIT,mpi_double_precision,0,mpi_comm_world,ios)

   if (ig_myrank == 0) then

      write(myfmt,'(i)') NTIME_INIT
      myfmt = "(a,i6,1x,"//trim(adjustl(myfmt))//"(e14.7,1x))"

      write(unit=IG_LST_UNIT,fmt='(" ",/,a)') "elapsed time for initialization"
      write(unit=IG_LST_UNIT,fmt='(      a)') "             init_mesh    init_gll_nodes  init_mpi_buff  init_medium    init_source    init_receiver    init_jaco     init_mass      init_dcs      init_snapshot   init_memory"

      do icpu = 1,ig_ncpu
         write(unit=IG_LST_UNIT,fmt=trim(myfmt)) "cpu ",icpu,(buffer_double_all_cpu((icpu-1)*NTIME_INIT+itime),itime=1,NTIME_INIT)
      enddo

   endif

   deallocate(buffer_double)
   deallocate(buffer_double_all_cpu)

!
!
!***********************************************************************************************************************
!***********************************************************************************************************************
!phase 2: time loop
!***********************************************************************************************************************
!***********************************************************************************************************************
   call mpi_barrier(mpi_comm_world,ios)

   start_time = mpi_wtime()

   if (LG_OUTPUT_CPUTIME) then 
      open(unit=get_newunit(unit_time),file=trim(cg_prefix)//".time.cpu."//trim(cg_myrank),status='replace')
   endif

   if (ig_myrank == 0) then
      write(IG_LST_UNIT,'(" ",/,a)') "starting time loop"
      write(IG_LST_UNIT,'(" ",/,a)') " -->time of simulation"
   endif
   
   do ig_idt = 1,ig_ndt

      !
      !
      !*****************************************************************************************************************
      !->time of the simulation
      !*****************************************************************************************************************
      rg_simu_current_time = (ig_idt-1)*rg_dt

      !
      !
      !*****************************************************************************************************************
      !->output time step in file *.lst
      !*****************************************************************************************************************
      if ( (ig_myrank == 0) .and. mod(ig_idt-1,1000) == 0 ) then
         write(IG_LST_UNIT,'(7X,E14.7)') rg_simu_current_time
         call flush(IG_LST_UNIT)
      endif

      !
      !
      !*****************************************************************************************************************
      !->compute cpu time for each cpu for each time steps
      !*****************************************************************************************************************
      if (LG_OUTPUT_CPUTIME) then 
        t_tloop = mpi_wtime()
      endif

      !
      !
      !*****************************************************************************************************************
      !->compute displacements at step n+1
      !*****************************************************************************************************************
      if (LG_OUTPUT_CPUTIME) then 
        t_newmark = mpi_wtime()
      endif
   
      call newmark_ini()
   
      if (LG_OUTPUT_CPUTIME) then
         t_newmark = mpi_wtime() - t_newmark
      endif

      !
      !
      !*****************************************************************************************************************
      !->compute external forces
      !*****************************************************************************************************************
      if (LG_OUTPUT_CPUTIME) then
        t_forext = mpi_wtime()
      endif
   
      call compute_external_force()
   
      if (LG_OUTPUT_CPUTIME) then
        t_forext = mpi_wtime() - t_forext
      endif      

      !
      !
      !*****************************************************************************************************************
      !->compute boundary absorption forces
      !*****************************************************************************************************************
      if (LG_OUTPUT_CPUTIME) then 
        t_forabs = mpi_wtime()
      endif
      
      if (lg_boundary_absorption) call compute_absorption_forces()
      
      if (LG_OUTPUT_CPUTIME) then
         t_forabs = mpi_wtime() - t_forabs
      endif

      !
      !
      !*****************************************************************************************************************
      !->compute internal forces for inner and outer elements
      !*****************************************************************************************************************
      if (LG_OUTPUT_CPUTIME) then
        t_forint_outer = mpi_wtime()
      endif
     
      if (IG_LAGRANGE_ORDER == 4) call compute_internal_forces_order4(1,ig_nhexa_outer)
      if (IG_LAGRANGE_ORDER == 5) call compute_internal_forces_order5(1,ig_nhexa_outer)
      if (IG_LAGRANGE_ORDER == 6) call compute_internal_forces_order6(1,ig_nhexa_outer)
   
      if (LG_OUTPUT_CPUTIME) then
         t_forint_outer = mpi_wtime() - t_forint_outer
      endif
     
      if (LG_OUTPUT_CPUTIME) then
         t_link_init = mpi_wtime()
      endif
    
      if (LG_ASYNC_MPI_COMM) call assemble_force_async_comm_init()
   
      if (LG_OUTPUT_CPUTIME) then
         t_link_init = mpi_wtime() - t_link_init
      endif 
      
      if (LG_OUTPUT_CPUTIME) then
         t_forint_inner = mpi_wtime()
      endif
      
      if (IG_LAGRANGE_ORDER == 4) call compute_internal_forces_order4(ig_nhexa_outer+1,ig_nhexa)
      if (IG_LAGRANGE_ORDER == 5) call compute_internal_forces_order5(ig_nhexa_outer+1,ig_nhexa)
      if (IG_LAGRANGE_ORDER == 6) call compute_internal_forces_order6(ig_nhexa_outer+1,ig_nhexa)
   
      if (LG_OUTPUT_CPUTIME) then
         t_forint_inner = mpi_wtime() - t_forint_inner
      endif

      !
      !
      !*****************************************************************************************************************
      !->assemble forces for linked cpu
      !*****************************************************************************************************************
      if (LG_OUTPUT_CPUTIME) then 
         t_linkfo = mpi_wtime()
      endif
      
      if (LG_ASYNC_MPI_COMM) then
         call assemble_force_async_comm_end()
      else
         call assemble_force_sync_comm()
      endif  
      
      if (LG_OUTPUT_CPUTIME) then 
         t_linkfo = mpi_wtime() - t_linkfo
      endif

      !      
      !      
      !*****************************************************************************************************************
      !->solve for acceleration and velocity at step n+1
      !*****************************************************************************************************************
      if (LG_OUTPUT_CPUTIME) then 
         t_resol = mpi_wtime()
      endif
   
      call newmark_end()
   
      if (LG_OUTPUT_CPUTIME) then 
         t_resol = mpi_wtime() - t_resol
      endif

      !
      !
      !*****************************************************************************************************************
      !->output snapshot of the surface of the domain + compute PGD, PGV and PGA
      !*****************************************************************************************************************
      if (LG_OUTPUT_CPUTIME) then 
         t_snafsu = mpi_wtime()
      endif

      if (lg_snapshot) call write_snapshot()

      if (LG_OUTPUT_CPUTIME) then
          t_snafsu = mpi_wtime() - t_snafsu
      endif

      !
      !
      !*****************************************************************************************************************
      !->output acceleration, velocity and displacement time history of receivers
      !*****************************************************************************************************************
      if (LG_OUTPUT_CPUTIME) then 
         t_outavd = mpi_wtime()
      endif
      if ( (ig_idt == 1) .or. mod((ig_idt-1),ig_receiver_saving_incr) == 0 ) then
         call write_receiver_output()
      endif
      if (LG_OUTPUT_CPUTIME) then 
         t_outavd = mpi_wtime() - t_outavd
      endif

      !
      !
      !*****************************************************************************************************************
      !->output cpu time for each cpu for each time step + detail of different subroutines
      !*****************************************************************************************************************
      if (LG_OUTPUT_CPUTIME) then

         t_tloop = mpi_wtime() - t_tloop

         write(unit_time,'(i10,1x,11(e22.15,1x))') ig_idt  & !1  time step
                                                  ,t_tloop               & !2  computation of entire time loop
                                                  ,t_newmark             & !3  computation of displacement (step n+1)
                                                  ,t_forext              & !4  computation of external forces
                                                  ,t_forabs              & !5  computation of boundary absorption forces
                                                  ,t_forint_outer        & !6  computation of internal forces of outer elements
                                                  ,t_forint_inner        & !7  computation of internal forces of inner elements
                                                  ,t_linkfo              & !8  computation of assemble forces
                                                  ,t_resol               & !9  computation of acceleration and velocity (step n+1)
                                                  ,t_snafsu              & !10 computation of free surface snapshot
                                                  ,t_outavd              & !11 computation of receiver time history
                                                  ,t_link_init             !12 assemble_force_async_comm_init

         call flush(unit_time)

      endif
   
   enddo!loop on time step

!
!
!*********************************************************************************************************************
!*********************************************************************************************************************
!phase 3: output elapsed time info and finalize computation
!*********************************************************************************************************************
!*********************************************************************************************************************

   if (lg_snapshot) then

      !***************************************************************************************************************
      !write peak ground motion at the free surface
      !***************************************************************************************************************
      call write_peak_ground_motion()

      !***************************************************************************************************************
      !write collection file if snapshot are saved in VTK format
      !***************************************************************************************************************
      if (LG_SNAPSHOT_VTK) call write_collection_vtk()

   endif

   if (LG_OUTPUT_CPUTIME) then 
      close(unit_time)
   endif

   !
   !
   !***************************************************************************************************************
   !compute elapsed time for time loop for all cpu
   !***************************************************************************************************************
   end_time = (mpi_wtime() - start_time)
   call mpi_reduce(end_time,dltim2,1,mpi_double_precision,mpi_max,0,mpi_comm_world,ios)

   !
   !
   !***************************************************************************************************************
   !cpu 0 gather elapsed time for time loop of each cpu
   !***************************************************************************************************************
   allocate(buffer_double_all_cpu(ig_ncpu))
   call mpi_gather(end_time,1,mpi_double_precision,buffer_double_all_cpu,1,mpi_double_precision,0,mpi_comm_world,ios)

   !
   !
   !***************************************************************************************************************
   !cpu 0 gather number of hexa of each cpu
   !***************************************************************************************************************
   allocate(buffer_integer(ig_ncpu))
   call mpi_gather(ig_nhexa,1,mpi_integer,buffer_integer,1,mpi_integer,0,mpi_comm_world,ios)

   !
   !
   !***************************************************************************************************************
   !output results in file *.lst
   !***************************************************************************************************************
   if (ig_myrank == 0) then
      write(IG_LST_UNIT,'(" ",/,a,e15.7,a)') "elapsed time for initialization        = ",dltim1," s"
      write(IG_LST_UNIT,'(a,e15.7,a)')       "elapsed time for time loop computation = ",dltim2," s"
      write(IG_LST_UNIT,'(a,e15.7,a)')       "total elapsed time for computation     = ",dltim1+dltim2," s"
   
      write(IG_LST_UNIT,'("",/,a,i0,a)'    ) "average time per time step and per hexa (order ",IG_LAGRANGE_ORDER,") for the simulation"
      do icpu = 1,ig_ncpu
         write(IG_LST_UNIT,'(a,i8,1x,e15.7,a)') " -->cpu ",icpu-1,buffer_double_all_cpu(icpu)/(dble(ig_ndt)*dble(buffer_integer(icpu)))," s"
      enddo
   endif
   
   call mpi_finalize(ios)
   stop

end program EFISPEC3D
