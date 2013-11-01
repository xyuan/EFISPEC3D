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
!!This file contains a module to initialize some global variable vectors and matrices.

!>@brief
!!This module contains subroutines to initialize some global variable vectors and matrices.
module mod_init_efi

   implicit none

   private

   public :: init_mpi
   public :: init_input_variables
   public :: init_gll_nodes
   public :: init_gll_nodes_coordinates
   public :: init_mass_matrix
   public :: init_jacobian_matrix_hexa
   public :: init_jacobian_matrix_quad

   contains

!>@brief
!!This subroutine initializes MPI
!>@return mod_global_variables::ig_ncpu
!>@return mod_global_variables::ig_myrank
!>@return mod_global_variables::cg_cpu_name
!>@return mod_global_variables::ig_cpu_name_len
!***********************************************************************************************************************************************************************************
   subroutine init_mpi()
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only :&
                                      error_stop&
                                     ,ig_ncpu&
                                     ,cg_ncpu&
                                     ,ig_myrank&
                                     ,cg_myrank&
                                     ,cg_cpu_name&
                                     ,ig_cpu_name_len
      
      implicit none

      integer            :: ios
      character(len=255) :: info
      
      call mpi_init(ios)
      if (ios /= 0) then
         write(info,'(a)') "error while initializing mpi"
         call error_stop(info)
      endif
      
      call mpi_comm_size(mpi_comm_world,ig_ncpu,ios)
      write(cg_ncpu,'(i6.6)') ig_ncpu
      
      call mpi_comm_rank(mpi_comm_world,ig_myrank,ios)
      write(cg_myrank,'(i6.6)') ig_myrank
      
      call mpi_get_processor_name(cg_cpu_name,ig_cpu_name_len,ios)
      
      return

!***********************************************************************************************************************************************************************************
   end subroutine init_mpi
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine initializes the simulation by writing header of listing file *.lst and reading variables of configuration file *.cfg
!>@return global variables read from file *.cfg needed to initialize a simulation
!***********************************************************************************************************************************************************************************
   subroutine init_input_variables()
!***********************************************************************************************************************************************************************************
      
      use mpi

      use mod_global_variables
      
      implicit none
  
      integer :: unit_input_file
      integer :: ios
      integer :: icpu
      integer :: pos
      integer :: imat
      logical :: is_material_defined = .false.
      logical :: is_file_exist       = .false.
      logical :: is_activate_dcs     = .false. 
      logical :: is_activate_sfs     = .false.

      character(len=MPI_MAX_PROCESSOR_NAME) :: cpu_name(ig_ncpu)
      character(len=100)                    :: buffer_line
      character(len=100)                    :: buffer_value
      character(len=100)                    :: keyword
      character(len=255)                    :: info
   
   !
   !->Check if order of polynomial is correct
      if (IG_LAGRANGE_ORDER < 4 .or. IG_LAGRANGE_ORDER > 6) then
         write(info,'(a)') "error while checking order of lagrange polynomial"
         call error_stop(info)
      endif
   
   !
   !->read prefix
      open(unit=100,file="prefix",iostat=ios)
      if(ios /= 0) then
         write(info,'(a)') "error while opening file prefix"
         call error_stop(info)
      endif
   
      read(unit=100,fmt=*),cg_prefix
      cg_prefix = trim(adjustl(cg_prefix))
      close(100)

   !
   !->cpu0 gather all the name of the other cpu
      cpu_name(:) = " "
      call mpi_gather(cg_cpu_name,MPI_MAX_PROCESSOR_NAME,mpi_character,cpu_name,MPI_MAX_PROCESSOR_NAME,mpi_character,0,mpi_comm_world,ios)
      if (ios /= 0) then
         write(info,'(a)') "error while gathering cpus' name"
         call error_stop(info)
      endif
   
      if (ig_myrank == 0) then

         open(unit=IG_LST_UNIT,file=trim(cg_prefix)//".lst") !IG_LST_UNIT = 10
         write(IG_LST_UNIT,'(a)') "************************************************************"
         write(IG_LST_UNIT,'(a)') "***                    EFISPEC3D                         ***"
         write(IG_LST_UNIT,'(a)') "***         (Elements FInis SPECtraux 3D)                ***"
         write(IG_LST_UNIT,'(a)') "***                                                      ***"
         write(IG_LST_UNIT,'(a)') "***                   version 1.0                        ***"
         write(IG_LST_UNIT,'(a)') "***                                                      ***"
         write(IG_LST_UNIT,'(a)') "***               http://efispec.free.fr                 ***"
         write(IG_LST_UNIT,'(a)') "***                                                      ***"
         write(IG_LST_UNIT,'(a)') "************************************************************"

   !
   !---->cpu name used for computation
         write(IG_LST_UNIT,'("",/,a)') "name of the cpus used for computation"
         do icpu = 1,ig_ncpu
            write(IG_LST_UNIT,'(2a)')   " --> ",trim(adjustl(cpu_name(icpu)))
         enddo
   !
   !---->communication type
         if (LG_ASYNC_MPI_COMM) then
            write(IG_LST_UNIT,'("",/,a)') "communication between cpus are asynchrone"
         else
            write(IG_LST_UNIT,'("",/,a)') "communication between cpus are synchrone"
         endif
   !
   !---->spectral element's polynomial order
         write(IG_LST_UNIT,'("",/,a,I0,a)') "elements use order ",IG_LAGRANGE_ORDER," polynomial"
         write(IG_LST_UNIT,'(     a,I0,a)') " --> ",IG_NGLL**3," GLL points per hexa"

      endif


      !
      !
      !******************************
      !open configuration file *.cfg
      !******************************
      open(unit=get_newunit(unit_input_file),file=trim(cg_prefix)//".cfg",status='old',iostat=ios)
      if (ios /= 0) then
         write(info,'(a)') "error in subroutine init_input_variables while opening file *.cfg"
         call error_stop(info)
      endif


      !
      !
      !*********************************
      !search for keyword in file *.cfg
      !*********************************
      do
      
         read(unit=unit_input_file,fmt='(a)',iostat=ios) buffer_line
         if (ios /= 0) exit
         
         buffer_line = trim(adjustl(buffer_line))
         pos         = scan(buffer_line,'=')
         if (pos >= 2) then
            keyword      = strlowcase(buffer_line(1:pos-1))
            buffer_value = buffer_line(pos+1:100)
         else
            keyword      = "na"
         endif
      
         selectcase(trim(sweep_blanks(keyword)))
      
            case("durationofsimulation")
               read(buffer_value,*) rg_simu_total_time
      
            case("timestep")
               read(buffer_value,*) rg_dt
               rg_dt2 = rg_dt**2
      
            case("receiversavingincrement")
               read(buffer_value,*) ig_receiver_saving_incr
      
            case("snapshotsavingincrement")
               read(buffer_value,*) ig_snapshot_saving_incr
      
            case("snapshotspaceincrement")
               read(buffer_value,*) rg_receiver_snapshot_dxdy
      
            case("snapshotdisplacement") 
               read(buffer_value,*) lg_snapshot_displacement
      
            case("snapshotvelocity")
               read(buffer_value,*) lg_snapshot_velocity
      
            case("snapshotacceleration")
               read(buffer_value,*) lg_snapshot_acceleration
      
            case("boundaryabsorption")
               read(buffer_value,*) lg_boundary_absorption
      
            case("numberofmaterial")
               read(buffer_value,*) ig_nmaterial
               if (.not.allocated(tg_elastic_material)) allocate(tg_elastic_material(ig_nmaterial))
               if (.not.allocated(ig_material_type))    allocate(ig_material_type(ig_nmaterial))
               if (.not.allocated(tg_viscoelastic_material) .and. LG_VISCO) allocate(tg_viscoelastic_material(ig_nmaterial))
   
            case("material")
               read(buffer_value,*) imat
               is_material_defined = .true.
               if (.not.LG_VISCO) then
                  ig_material_type(imat) = 1
               else
                  ig_material_type(imat) = 2
               endif
      
            case("vs")
               read(buffer_value,*) tg_elastic_material(imat)%svel
      
            case("vp")
               read(buffer_value,*) tg_elastic_material(imat)%pvel
      
            case("rho")
               read(buffer_value,*) tg_elastic_material(imat)%dens
      
            case("qf")
               if (LG_VISCO) then
                  read(buffer_value,*) tg_viscoelastic_material(imat)%freq
               endif
      
            case("qs")
               if (LG_VISCO) then
                  read(buffer_value,*) tg_viscoelastic_material(imat)%qs
               endif
      
            case("qp")
               if (LG_VISCO) then
                  read(buffer_value,*) tg_viscoelastic_material(imat)%qp
               endif
      
            case default
      
         endselect
      
      enddo
   
      close(unit_input_file)


      !
      !
      !*************************************
      !set default value of input variables
      !*************************************
      ig_ndt = int(rg_simu_total_time/rg_dt) + 1
      
      if (ig_receiver_saving_incr == 0) then

         ig_receiver_saving_incr = ig_ndt + 10

      endif
      
      if ( (ig_snapshot_saving_incr == 0) .or. (.not.lg_snapshot_displacement .and. .not.lg_snapshot_velocity .and. .not.lg_snapshot_acceleration) ) then

         ig_snapshot_saving_incr  = ig_ndt + 10
         lg_snapshot              = .false.
         lg_snapshot_displacement = .false.
         lg_snapshot_velocity     = .false.
         lg_snapshot_acceleration = .false.

      else

         lg_snapshot  = .true.

      endif
   
      inquire(file=trim(cg_prefix)//".dcs", exist=is_file_exist)
      if(is_file_exist) then
         is_activate_dcs = .true.
      endif
   
      inquire(file=trim(cg_prefix)//".sfs", exist=is_file_exist)
      if(is_file_exist) then
         is_activate_sfs = .true.
      endif
   
      if (.not.(is_activate_dcs) .and. .not.(is_activate_sfs)) then
         write(info,'(a)') "error in subroutine init_input_variables: no source defined"
         call error_stop(info)
      endif


      !
      !
      !********************************************
      !check for illegal value of input variables
      !********************************************
      if (rg_simu_total_time <= 1.0e-10) then
         write(info,'(a)') "error in subroutine init_input_variables: duration of simulation must be greater than zero"
         call error_stop(info)
      endif
      
      if (rg_dt <= 1.0e-10) then
         write(info,'(a)') "error in subroutine init_input_variables: time step must be greater than zero"
         call error_stop(info)
      endif

      do imat = 1,ig_nmaterial
         if (tg_elastic_material(imat)%svel <= 1.0e-10) then
            write(info,'(a)') "error in subroutine init_input_variables: S-wave velocity must be greater than zero"
            call error_stop(info)
         endif
         if (tg_elastic_material(imat)%pvel <= 1.0e-10) then
            write(info,'(a)') "error in subroutine init_input_variables: P-wave velocity must be greater than zero"
            call error_stop(info)
         endif
         if (tg_elastic_material(imat)%dens <= 1.0e-10) then
            write(info,'(a)') "error in subroutine init_input_variables: density must be greater than zero"
            call error_stop(info)
         endif
      enddo
      
      if (.not.(is_material_defined)) then
         write(info,'(a)') "error in subroutine init_input_variables: no material defined"
         call error_stop(info)
      endif



      !
      !
      !*******************************************************************
      !write in listing file if the simulation is elastic or viscoelastic
      !*******************************************************************
      if (.not.LG_VISCO) then
         if (ig_myrank == 0) then
            write(IG_LST_UNIT,'(" ",/,a)') "simulation is elastic"
         endif
      else
         if (ig_myrank == 0) then
            write(IG_LST_UNIT,'(" ",/,a,I3,a)') "simulation is viscoelastic using ",IG_NRELAX," memory variables"
         endif
      endif
      !
      !
      !*************************************************
      !print value of input variables in listing file
      !*************************************************
      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(" ",/,a,i10)')                  "number of time step    = ",ig_ndt
         write(IG_LST_UNIT,'(a,f10.6)')                      "size   of time step    = ",rg_dt
         write(IG_LST_UNIT,'(a,f10.6)')                      "duration of simulation = ",rg_simu_total_time
         write(IG_LST_UNIT,'(a)')                            "                       "
         if (ig_receiver_saving_incr < ig_ndt) then
                              write(IG_LST_UNIT,'(a,i8,a)' ) "receivers' time history saved every ",ig_receiver_saving_incr," time steps"
         else
                              write(IG_LST_UNIT,'(a,i8,a)' ) "receivers' time history are not saved"
         endif
         if (ig_snapshot_saving_incr < ig_ndt) then
                              write(IG_LST_UNIT,'(a,i8,a)' ) "surface snapshot saved every        ",ig_snapshot_saving_incr," time steps"
         else
                              write(IG_LST_UNIT,'(a,i8,a)' ) "surface snapshot are not saved"
         endif
         if (lg_snapshot_displacement) then
                              write(IG_LST_UNIT,'(a)')       " --> displacement snapshot are     saved"
         else
                              write(IG_LST_UNIT,'(a)')       " --> displacement snapshot are not saved"
         endif
         if (lg_snapshot_velocity) then
                              write(IG_LST_UNIT,'(a)')       " --> velocity     snapshot are     saved"
         else
                              write(IG_LST_UNIT,'(a)')       " --> velocity     snapshot are not saved"
         endif
         if (lg_snapshot_acceleration) then 
                              write(IG_LST_UNIT,'(a)')       " --> acceleration snapshot are     saved"
         else
                              write(IG_LST_UNIT,'(a)')       " --> acceleration snapshot are not saved"
         endif
      endif
 
      return
!***********************************************************************************************************************************************************************************
   end subroutine init_input_variables
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine computes GLL nodes abscissa and weight in the reference domain [-1:1] as well as derivative of Lagrange polynomials at GLL nodes.
!!
!>@reference C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
!!
!>@authors Initially written for Matlab by G. von Winckel - 04/17/2004. Translated in FORTRAN for EFISPEC3D by F. De Martin
!!
!>@return mod_global_variables::rg_gll_abscissa       : GLL node coordinates @f$\xi,\eta,\zeta@f$ in the reference domain [-1:1]
!>@return mod_global_variables::rg_gll_weight         : GLL node weights
!>@return mod_global_variables::rg_gll_lagrange_deriv : derivative of Lagrange polynomials at GLL nodes
!>@return mod_global_variables::rg_gll_abscissa_dist  : inverse of distance between GLL nodes
!***********************************************************************************************************************************************************************************
   subroutine init_gll_nodes()
!***********************************************************************************************************************************************************************************
   
      use mpi
      
      use mod_global_variables, only : &
                                       cg_prefix&
                                      ,rg_gll_abscissa_dist&
                                      ,rg_gll_abscissa&
                                      ,rg_gll_weight&
                                      ,rg_gll_lagrange_deriv&
                                      ,IG_NGLL&
                                      ,IG_LST_UNIT&
                                      ,RG_PI&
                                      ,ig_myrank&
                                      ,IG_LAGRANGE_ORDER
      
      use mod_lagrange        , only : lagrad
      
      implicit none

      real, parameter                  :: EPS = 1.0E-7
      real, allocatable ,dimension(:)  :: xold
      real, dimension(IG_NGLL,IG_NGLL) :: vandermonde_matrix
      integer                          :: ngll
      integer                          :: i
      integer                          :: j
      integer                          :: k
      integer                          :: n
      
      !
      !
      !higher integration order among different group is taken as the integration order
      n = IG_LAGRANGE_ORDER

 
      !
      !
      !ngll is the number of integration node in 1d along the segment [-1;1]
      ngll = n + 1


      !
      !
      !*******************************
      !find gll abscissa and weight
      !*******************************
      !
      !use the chebyshev-gauss-lobatto nodes as the first guess
      do i = 1,ngll
         rg_gll_abscissa(i) = cos(RG_PI*real(i-1)/real(n))
      enddo
      !
      !initialize the legendre vandermonde matrix
      vandermonde_matrix(:,:) = 0.0
      !
      !compute p_{n} using the recursion relation. compute its first and second derivatives and update rg_gll_abscissa using the newton-raphson method
      allocate(xold(ngll))
      xold(:) = 2.0
      do while (maxval(abs(rg_gll_abscissa-xold)) > EPS)
      
         xold(:)                 = rg_gll_abscissa(:)
         vandermonde_matrix(:,1) = 1.0
         vandermonde_matrix(:,2) = rg_gll_abscissa(:)
      
         do k=2,n
            vandermonde_matrix(:,k+1) = ( (2.0*real(k)-1.0)*rg_gll_abscissa(:)*vandermonde_matrix(:,k)-(real(k)-1.0)*vandermonde_matrix(:,k-1) )/real(k)
         enddo
      
         rg_gll_abscissa(:) = xold(:)-( rg_gll_abscissa(:)*vandermonde_matrix(:,ngll)-vandermonde_matrix(:,n) )/( ngll*vandermonde_matrix(:,ngll) )
      
      enddo
      !
      !compute corresponding weight
      rg_gll_weight(:) = 2.0/(n*ngll*vandermonde_matrix(:,ngll)**2)
      
      !
      !
      !*******************************************************************************
      !write results rg_gll_abscissa,rg_gll_weight in file *.lst
      !*******************************************************************************
      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(" ",/,a)') " --> abscissa of GLLs - weight of GLLs"
         do i = 1,ngll
            write(IG_LST_UNIT,fmt='(2(7X,f10.7))') rg_gll_abscissa(i),rg_gll_weight(i)
         enddo
      endif
      
      !
      !
      !*******************************************************************************
      !compute subtraction between gll nodes (used in subroutines lagrap and lagrad)
      !*******************************************************************************
      do i = 1,IG_NGLL
         do j = 1,IG_NGLL
            rg_gll_abscissa_dist(j,i) = 0.0
            if (i /= j) rg_gll_abscissa_dist(j,i) =1.0/(rg_gll_abscissa(i) - rg_gll_abscissa(j))
         enddo
      enddo
      
      !
      !
      !***************************************************
      !compute derivative of lagrange poly at GLL nodes
      !***************************************************
      !
      !notation convention: rg_gll_lagrange_deriv(j,i) = h'_j(xgll_i)
      do i = 1,IG_NGLL
         do j = 1,IG_NGLL
           rg_gll_lagrange_deriv(j,i) = lagrad(j,rg_gll_abscissa(i),IG_NGLL)
         enddo
      enddo
      
      return

!***********************************************************************************************************************************************************************************
   end subroutine init_gll_nodes
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine computes GLL nodes @f$x,y,z@f$-coordinates in the physical domain, cartesian right-handed coordinate system.
!>@return mod_global_variables::rg_gll_coordinate : GLL nodes @f$x,y,z@f$-coordinates in cpu myrank
!***********************************************************************************************************************************************************************************
   subroutine init_gll_nodes_coordinates()
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only :&
                                      rg_gll_coordinate&
                                     ,ig_ngll_total&
                                     ,ig_nhexa&
                                     ,IG_NGLL&
                                     ,IG_NDOF&
                                     ,ig_hexa_gll_glonum&
                                     ,ig_hexa_nnode&
                                     ,ig_hexa_gnode_xiloc&
                                     ,ig_hexa_gnode_etloc&
                                     ,ig_hexa_gnode_zeloc&
                                     ,rg_gll_abscissa&
                                     ,ig_line_nnode&
                                     ,ig_hexa_gnode_glonum&
                                     ,rg_gnode_x&
                                     ,rg_gnode_y&
                                     ,rg_gnode_z&
                                     ,cg_prefix&
                                     ,cg_myrank&
                                     ,LG_OUTPUT_DEBUG_FILE&
                                     ,error_stop

      use mod_init_memory
      
      use mod_lagrange, only : lagrap_geom
      
      implicit none

      real                                      :: lag_xi
      real                                      :: lag_eta
      real                                      :: lag_zeta
                                                
      real                                      :: xgnode
      real                                      :: ygnode
      real                                      :: zgnode
                                                
      integer                                   :: ihexa
      integer                                   :: k
      integer                                   :: l
      integer                                   :: m
      integer                                   :: n
      integer                                   :: igll
      integer                                   :: ios
                                              
      integer(kind=1), dimension(ig_ngll_total) :: maskeq
                                              
      character(len=255)                        :: info 

      !
      !
      !*************************************************************************************
      !compute global coordinate of gll nodes once and for all
      !*************************************************************************************

      ios = init_array_real(rg_gll_coordinate,ig_ngll_total,IG_NDOF,"rg_gll_coordinate")
   
      maskeq(:) = 0 !-->coordinate not computed yet
      
      do ihexa = 1,ig_nhexa
   !
   !---->compute global x,y,z coordinate of gll nodes
         do k = 1,IG_NGLL        !zeta
            do l = 1,IG_NGLL     !eta
               do m = 1,IG_NGLL  !xi
   
                  igll = ig_hexa_gll_glonum(m,l,k,ihexa)
   
                  if (maskeq(igll) == 0) then
   
                     do n = 1,ig_hexa_nnode
   
                        lag_xi   = lagrap_geom(ig_hexa_gnode_xiloc(n),rg_gll_abscissa(m),ig_line_nnode)
                        lag_eta  = lagrap_geom(ig_hexa_gnode_etloc(n),rg_gll_abscissa(l),ig_line_nnode)
                        lag_zeta = lagrap_geom(ig_hexa_gnode_zeloc(n),rg_gll_abscissa(k),ig_line_nnode)
   
                        xgnode   = rg_gnode_x(ig_hexa_gnode_glonum(n,ihexa))
                        ygnode   = rg_gnode_y(ig_hexa_gnode_glonum(n,ihexa))
                        zgnode   = rg_gnode_z(ig_hexa_gnode_glonum(n,ihexa))
      
                        rg_gll_coordinate(1,igll) = rg_gll_coordinate(1,igll) + lag_xi*lag_eta*lag_zeta*xgnode
            
                        rg_gll_coordinate(2,igll) = rg_gll_coordinate(2,igll) + lag_xi*lag_eta*lag_zeta*ygnode
   
                        rg_gll_coordinate(3,igll) = rg_gll_coordinate(3,igll) + lag_xi*lag_eta*lag_zeta*zgnode
   
                     enddo
   
                  endif
   
                  maskeq(igll) = 1
   
               enddo
            enddo
         enddo
      enddo
   
      if (LG_OUTPUT_DEBUG_FILE) then
         open(unit=100,file="debug."//trim(cg_prefix)//".global.element.gll.coordinates.cpu."//trim(cg_myrank))
   
         do ihexa = 1,ig_nhexa
            do k = 1,IG_NGLL
               do l = 1,IG_NGLL
                  do m = 1,IG_NGLL
   
                     igll = ig_hexa_gll_glonum(m,l,k,ihexa)
   
                     write(100,'(3(f15.3,1x),i10)') rg_gll_coordinate(1,igll)&
                                                   ,rg_gll_coordinate(2,igll)&
                                                   ,rg_gll_coordinate(3,igll)&
                                                   ,ig_hexa_gll_glonum(m,l,k,ihexa)
                  enddo
               enddo
            enddo
         enddo
   
         close(100)
      endif
   
      return
   
!***********************************************************************************************************************************************************************************
   end subroutine init_gll_nodes_coordinates
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine computes and assembles the mass matrix @f$M@f$ of the system @f$ M\ddot{U} + C\dot{U} + KU = F^{ext} @f$.
!>@return mod_global_variables::rg_gll_mass_matrix : diagonal mass matrix in cpu myrank
!***********************************************************************************************************************************************************************************
   subroutine init_mass_matrix()
!***********************************************************************************************************************************************************************************
      
      use mpi
      
      use mod_global_variables, only : &
                                       IG_NGLL&
                                      ,ig_ngll_total&
                                      ,rg_gll_mass_matrix&
                                      ,rg_gll_weight&
                                      ,ig_hexa_gll_glonum&
                                      ,ig_nhexa&
                                      ,ig_hexa_material_number&
                                      ,rg_hexa_gll_jacobian_det&
                                      ,ig_ncpu&
                                      ,tg_cpu_neighbor&
                                      ,ig_myrank&
                                      ,rg_hexa_gll_rho&
                                      ,error_stop&
                                      ,ig_mpi_buffer_sizemax&
                                      ,ig_ncpu_neighbor&
                                      ,IG_LST_UNIT

      use mod_init_memory
      
      implicit none

      real   , dimension(ig_mpi_buffer_sizemax)                  :: dummyr
      real   , dimension(ig_mpi_buffer_sizemax,ig_ncpu_neighbor) :: dlxmas
                                                 
      integer, dimension(MPI_STATUS_SIZE)                        :: statut
      integer                                                    :: ios
      integer                                                    :: ngll
      integer                                                    :: icon
      integer                                                    :: ihexa
      integer                                                    :: igll
      integer                                                    :: kgll
      integer                                                    :: lgll
      integer                                                    :: mgll
                                                                 
      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(" ",/,a)') "computing mass matrix..."
         call flush(IG_LST_UNIT)
      endif

      !
      !
      !*****************************************************************************************************************************************
      !allocate mass matrix array
      !*****************************************************************************************************************************************
      ios = init_array_real(rg_gll_mass_matrix,ig_ngll_total,"rg_gll_mass_matrix")

      !
      !
      !*****************************************************************************************************************************************
      !compute mass matrix for cpu 'myrank'
      !*****************************************************************************************************************************************
      do ihexa = 1,ig_nhexa
         do kgll = 1,IG_NGLL        !zeta
            do lgll = 1,IG_NGLL     !eta
               do mgll = 1,IG_NGLL  !xi
   
                  igll = ig_hexa_gll_glonum(mgll,lgll,kgll,ihexa)
   
                  rg_gll_mass_matrix(igll) = rg_gll_mass_matrix(igll)&
                                           + rg_hexa_gll_jacobian_det(mgll,lgll,kgll,ihexa)*rg_hexa_gll_rho(mgll,lgll,kgll,ihexa)&
                                           * rg_gll_weight(mgll)*rg_gll_weight(lgll)*rg_gll_weight(kgll)
               enddo
            enddo
         enddo
      enddo
      !
      !
      !*****************************************************************************************************************************************
      !fill buffer dlxmas with gll to be send
      !*****************************************************************************************************************************************
      do icon = 1,ig_ncpu_neighbor
         do igll = 1,tg_cpu_neighbor(icon)%ngll
   
            kgll = tg_cpu_neighbor(icon)%gll_send(igll)
            dlxmas(igll,icon) = rg_gll_mass_matrix(kgll)
   
         enddo
      enddo
      !
      !
      !*****************************************************************************************************************************************
      !assemble mass matrix on all cpus once and for all
      !*****************************************************************************************************************************************
      do icon = 1,ig_ncpu_neighbor
      
         ngll = tg_cpu_neighbor(icon)%ngll
      
         call mpi_sendrecv(dlxmas(1,icon),ngll,mpi_real,tg_cpu_neighbor(icon)%icpu,100,dummyr,ngll,mpi_real,tg_cpu_neighbor(icon)%icpu,100,mpi_comm_world,statut,ios)
      
         do igll = 1,ngll
   
            kgll = tg_cpu_neighbor(icon)%gll_recv(igll)
            rg_gll_mass_matrix(kgll) = rg_gll_mass_matrix(kgll) + dummyr(igll)
   
         enddo
      
      enddo
      !
      !
      !
      !******************************************************************************************************************
      !Compute the inverse of mass matrix
      !******************************************************************************************************************
      rg_gll_mass_matrix(:) = 1.0/rg_gll_mass_matrix(:)
      
      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(a)') "done"
         call flush(IG_LST_UNIT)
      endif
      
      deallocate(rg_hexa_gll_rho)

      return
      
!***********************************************************************************************************************************************************************************
   end subroutine init_mass_matrix
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine computes Jacobian matrix and its determinant for hexahedron elements.
!>@return rg_hexa_gll_jacobian_det : determinant of jacobian matrix at GLL nodes of hexahedron elements in cpu myrank
!>@return rg_hexa_gll_dxidx        : @f$ \frac{\partial \xi  }{\partial x} @f$ at GLL nodes of hexahedron elements
!>@return rg_hexa_gll_dxidy        : @f$ \frac{\partial \xi  }{\partial y} @f$ at GLL nodes of hexahedron elements
!>@return rg_hexa_gll_dxidz        : @f$ \frac{\partial \xi  }{\partial z} @f$ at GLL nodes of hexahedron elements
!>@return rg_hexa_gll_detdx        : @f$ \frac{\partial \eta }{\partial x} @f$ at GLL nodes of hexahedron elements 
!>@return rg_hexa_gll_detdy        : @f$ \frac{\partial \eta }{\partial y} @f$ at GLL nodes of hexahedron elements
!>@return rg_hexa_gll_detdz        : @f$ \frac{\partial \eta }{\partial z} @f$ at GLL nodes of hexahedron elements
!>@return rg_hexa_gll_dzedx        : @f$ \frac{\partial \zeta}{\partial x} @f$ at GLL nodes of hexahedron elements  
!>@return rg_hexa_gll_dzedy        : @f$ \frac{\partial \zeta}{\partial y} @f$ at GLL nodes of hexahedron elements
!>@return rg_hexa_gll_dzedz        : @f$ \frac{\partial \zeta}{\partial z} @f$ at GLL nodes of hexahedron elements
!***********************************************************************************************************************************************************************************
   subroutine init_jacobian_matrix_hexa()
!***********************************************************************************************************************************************************************************

      use mpi
   
      use mod_global_variables, only : IG_NGLL&
                                      ,rg_gll_abscissa&
                                      ,rg_gll_weight&
                                      ,rg_gnode_x&
                                      ,rg_gnode_y&
                                      ,rg_gnode_z&
                                      ,ig_hexa_gnode_xiloc&
                                      ,ig_hexa_gnode_etloc&
                                      ,ig_hexa_gnode_zeloc&
                                      ,ig_nhexa&
                                      ,ig_line_nnode&
                                      ,rg_hexa_gll_jacobian_det&
                                      ,rg_hexa_gll_dxidx&
                                      ,rg_hexa_gll_dxidy&
                                      ,rg_hexa_gll_dxidz&
                                      ,rg_hexa_gll_detdx&
                                      ,rg_hexa_gll_detdy&
                                      ,rg_hexa_gll_detdz&
                                      ,rg_hexa_gll_dzedx&
                                      ,rg_hexa_gll_dzedy&
                                      ,rg_hexa_gll_dzedz&
                                      ,ig_hexa_nnode&
                                      ,ig_myrank&
                                      ,error_stop&
                                      ,IG_LST_UNIT

      use mod_init_memory

      use mod_jacobian

      use mod_lagrange      , only :&
                                    lagrap_geom&
                                   ,lagrad_geom

      implicit none
 
      real ,parameter    :: EPS = 1.0e-8
      real               :: xi
      real               :: eta
      real               :: zeta
      real               :: dxidx 
      real               :: dxidy
      real               :: dxidz
      real               :: detdx
      real               :: detdy
      real               :: detdz
      real               :: dzedx
      real               :: dzedy
      real               :: dzedz
      real               :: det
      
      integer            :: ihexa
      integer            :: k
      integer            :: l
      integer            :: m
      integer            :: ios

      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(" ",/,a)') "computing jacobian matrix in hexa..."
         call flush(IG_LST_UNIT)
      endif

      !
      !
      !********************************************************************************************************************
      !initialize memory
      !********************************************************************************************************************
      ios = init_array_real(rg_hexa_gll_jacobian_det,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_jacobian_det")

      ios = init_array_real(rg_hexa_gll_dxidx,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_dxidx")

      ios = init_array_real(rg_hexa_gll_dxidy,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_dxidy")

      ios = init_array_real(rg_hexa_gll_dxidz,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_dxidz")

      ios = init_array_real(rg_hexa_gll_detdx,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_detdx")

      ios = init_array_real(rg_hexa_gll_detdy,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_detdy")

      ios = init_array_real(rg_hexa_gll_detdz,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_detdz")

      ios = init_array_real(rg_hexa_gll_dzedx,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_dzedx")

      ios = init_array_real(rg_hexa_gll_dzedy,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_dzedy")

      ios = init_array_real(rg_hexa_gll_dzedz,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_dzedz")


      !
      !
      !********************************************************************************************************************
      !initialize jacobian matrix at each GLL nodes of each hexahedron elements
      !********************************************************************************************************************
      do ihexa = 1,ig_nhexa
   
         do k = 1,IG_NGLL
            do l = 1,IG_NGLL
               do m = 1,IG_NGLL
   
                  xi   = rg_gll_abscissa(m)
                  eta  = rg_gll_abscissa(l)
                  zeta = rg_gll_abscissa(k)
   
                  call compute_hexa_jacobian(ihexa,xi,eta,zeta,dxidx,dxidy,dxidz,detdx,detdy,detdz,dzedx,dzedy,dzedz,det)
   
                  rg_hexa_gll_jacobian_det(m,l,k,ihexa) = det

                  rg_hexa_gll_dxidx(m,l,k,ihexa) = dxidx 
                  rg_hexa_gll_dxidy(m,l,k,ihexa) = dxidy
                  rg_hexa_gll_dxidz(m,l,k,ihexa) = dxidz

                  rg_hexa_gll_detdx(m,l,k,ihexa) = detdx
                  rg_hexa_gll_detdy(m,l,k,ihexa) = detdy
                  rg_hexa_gll_detdz(m,l,k,ihexa) = detdz

                  rg_hexa_gll_dzedx(m,l,k,ihexa) = dzedx
                  rg_hexa_gll_dzedy(m,l,k,ihexa) = dzedy
                  rg_hexa_gll_dzedz(m,l,k,ihexa) = dzedz
   
               enddo
            enddo
         enddo
   
      enddo

      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(a)') "done"
         call flush(IG_LST_UNIT)
      endif
   
      return
   
!***********************************************************************************************************************************************************************************
   end subroutine init_jacobian_matrix_hexa
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine computes the determinant of Jacobian matrix and normal unit vector of quadrangle elements.
!>@return rg_quadp_gll_jaco_det : determinant of jacobian matrix at GLL nodes of quadrangle elements in cpu myrank
!>@return rg_quadp_gll_normal   : normal unit vector at GLL nodes of quadrangle elements
!***********************************************************************************************************************************************************************************
   subroutine init_jacobian_matrix_quad()
!***********************************************************************************************************************************************************************************

      use mpi
   
      use mod_global_variables, only : IG_NGLL&
                                      ,rg_gll_abscissa&
                                      ,rg_gll_weight&
                                      ,rg_gnode_x&
                                      ,rg_gnode_y&
                                      ,rg_gnode_z&
                                      ,ig_line_nnode&
                                      ,ig_quad_gnode_xiloc&
                                      ,ig_quad_gnode_etloc&
                                      ,ig_nquad_parax&
                                      ,rg_quadp_gll_jaco_det&
                                      ,rg_quadp_gll_normal&
                                      ,ig_quad_nnode&
                                      ,ig_myrank&
                                      ,error_stop&
                                      ,ig_quadp_gnode_glonum&
                                      ,IG_LST_UNIT&
                                      ,lg_boundary_absorption

      use mod_init_memory

      use mod_lagrange      , only :&
                                    lagrap_geom&
                                   ,lagrad_geom

      implicit none
 
      real ,parameter    :: EPS = 1.0e-8
      real               :: xxi
      real               :: xet
      real               :: yxi
      real               :: yet
      real               :: zxi
      real               :: zet
      real               :: det
      real               :: inv_det
      real               :: crprx
      real               :: crpry
      real               :: crprz
      real               :: vn(3)
      real               :: lag_deriv_xi
      real               :: lag_deriv_et
      real               :: lag_xi
      real               :: lag_et
      real               :: xgnode
      real               :: ygnode
      real               :: zgnode
      
      integer            :: iquad
      integer            :: k
      integer            :: l
      integer            :: m
      integer            :: nas
      integer            :: ios

      !
      !
      !*********************************************************************************************************************
      !if there are quad (eg, paraxial elt): compute jacobian and derivative of local coordinate with respect to global ones
      !*********************************************************************************************************************
   
      if (lg_boundary_absorption .and. ig_nquad_parax > 0) then

         if (ig_myrank == 0) then
            write(IG_LST_UNIT,'(" ",/,a)') "computing jacobian matrix in quad..."
            call flush(IG_LST_UNIT)
         endif

         !
         !
         !********************************************************************************************************************
         !initialize memory
         !********************************************************************************************************************
         ios = init_array_real(rg_quadp_gll_jaco_det,ig_nquad_parax,IG_NGLL,IG_NGLL,"rg_quadp_gll_jaco_det")

         ios = init_array_real(rg_quadp_gll_normal,ig_nquad_parax,IG_NGLL,IG_NGLL,3,"rg_quadp_gll_normal")
      

         !
         !
         !********************************************************************************************************************
         !initialize jacobian matrix and normal vector at each GLL nodes of each quadrangle elements
         !********************************************************************************************************************
         do iquad = 1,ig_nquad_parax
            do k = 1,IG_NGLL
               do l = 1,IG_NGLL
   
                  xxi = 0.0
                  xet = 0.0
                  yxi = 0.0
                  yet = 0.0
                  zxi = 0.0
                  zet = 0.0
   
                  do m = 1,ig_quad_nnode
   
                     lag_deriv_xi = lagrad_geom(ig_quad_gnode_xiloc(m),rg_gll_abscissa(l),ig_line_nnode)
                     lag_deriv_et = lagrad_geom(ig_quad_gnode_etloc(m),rg_gll_abscissa(k),ig_line_nnode)
   
                     lag_xi       = lagrap_geom(ig_quad_gnode_xiloc(m),rg_gll_abscissa(l),ig_line_nnode)
                     lag_et       = lagrap_geom(ig_quad_gnode_etloc(m),rg_gll_abscissa(k),ig_line_nnode)
   
                     xgnode       = rg_gnode_x(ig_quadp_gnode_glonum(m,iquad))
                     ygnode       = rg_gnode_y(ig_quadp_gnode_glonum(m,iquad))
                     zgnode       = rg_gnode_z(ig_quadp_gnode_glonum(m,iquad))
   
                     xxi          = xxi + lag_deriv_xi*lag_et*xgnode
                     yxi          = yxi + lag_deriv_xi*lag_et*ygnode
                     zxi          = zxi + lag_deriv_xi*lag_et*zgnode
                                  
                     xet          = xet + lag_xi*lag_deriv_et*xgnode
                     yet          = yet + lag_xi*lag_deriv_et*ygnode
                     zet          = zet + lag_xi*lag_deriv_et*zgnode
   
                  enddo
   
                  crprx   = yxi*zet - zxi*yet
                  crpry   = zxi*xet - xxi*zet
                  crprz   = xxi*yet - yxi*xet
   
                  det     = sqrt(crprx**2 + crpry**2 + crprz**2) !the determinant of the quad element is the magnitude of the cross-product
                  inv_det = 1.0/det
   
                  if (det < EPS) then
                     write(*,*) "***error in jacobi quad: det null or negative, det = ",det
                     write(*,*) "***element ",iquad, "cpu ",ig_myrank
                     call mpi_abort(mpi_comm_world,100,ios)
                     stop
                  endif
   !!
   !!------------>determinant of the quad element
                  rg_quadp_gll_jaco_det(l,k,iquad) = det
   !!
   !!------------>unit vector normal to the quad at the node k,l
                  vn(1) = crprx*inv_det
                  vn(2) = crpry*inv_det
                  vn(3) = crprz*inv_det
   
                  rg_quadp_gll_normal(1,l,k,iquad) = vn(1) !x
                  rg_quadp_gll_normal(2,l,k,iquad) = vn(2) !y
                  rg_quadp_gll_normal(3,l,k,iquad) = vn(3) !z
      
               enddo
            enddo
         enddo

         if (ig_myrank == 0) then
            write(IG_LST_UNIT,'(a)') "done"
            call flush(IG_LST_UNIT)
         endif

      endif
   
      return
   
!***********************************************************************************************************************************************************************************
   end subroutine init_jacobian_matrix_quad
!***********************************************************************************************************************************************************************************

end module mod_init_efi
