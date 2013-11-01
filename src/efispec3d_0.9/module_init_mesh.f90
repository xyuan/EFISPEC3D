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
!!This file contains a module to read mesh and create GLL nodes global numbering.

!>@brief
!!This module contains subroutines to read mesh files and creates GLL nodes global numbering in cpu myrank.
module mod_init_mesh
   
   implicit none

   private

!>size of buffer array to read file *.dat
   integer, parameter :: BUFFER_READ_SIZE = 8192                                                     
!>buffer array to read file *.dat
   integer, dimension(BUFFER_READ_SIZE) :: buffer 

   public  :: init_mesh
   private :: init_gll_number
   private :: propagate_gll_nodes_quad
   private :: propagate_gll_nodes_face
   private :: propagate_gll_nodes_edge
   private :: propagate_gll_nodes_corner
   private :: init_element
   private :: read_bin_int
   
   contains

!
!
!>@brief
!>This subroutine reads mesh files *.dat for cpu myrank and creates GLL numbering of hexahedron and quadrangle elements for cpu myrank
!***********************************************************************************************************************************************************************************
    subroutine init_mesh()
!***********************************************************************************************************************************************************************************
 
      use mpi
 
      use mod_global_variables, only : &
                                       ig_nhexa&
                                      ,ig_nhexa_outer&
                                      ,ig_nhexa_inner&
                                      ,ig_nquad_parax&
                                      ,ig_nquad_fsurf&
                                      ,ig_mesh_nnode&
                                      ,ig_hexa_nnode&
                                      ,ig_quad_nnode&
                                      ,ig_hexa_gnode_glonum&
                                      ,ig_quadp_gnode_glonum&
                                      ,ig_quadf_gnode_glonum&
                                      ,rg_gnode_x&
                                      ,rg_gnode_y&
                                      ,rg_gnode_z&
                                      ,rg_mesh_xmax&
                                      ,rg_mesh_ymax&
                                      ,rg_mesh_zmax&
                                      ,rg_mesh_xmin&
                                      ,rg_mesh_ymin&
                                      ,rg_mesh_zmin&
                                      ,rg_gll_displacement&
                                      ,rg_gll_velocity&
                                      ,rg_gll_acceleration&
                                      ,rg_gll_acctmp&
                                      ,cg_prefix&
                                      ,cg_myrank&
                                      ,cg_ncpu&
                                      ,ig_myrank&
                                      ,ig_ncpu_neighbor&
                                      ,tg_cpu_neighbor&
                                      ,ig_medium_type&
                                      ,ig_hexa_gll_glonum&
                                      ,ig_quadp_gll_glonum&
                                      ,ig_quadf_gll_glonum&
                                      ,ig_hexa_material_number&
                                      ,ig_quadp_neighbor_hexa&
                                      ,ig_quadp_neighbor_hexaface&
                                      ,IG_NGLL&
                                      ,ig_ngll_total&
                                      ,get_newunit&
                                      ,error_stop&
                                      ,ig_cpu_neighbor_info&
                                      ,ig_nneighbor_all_kind&
                                      ,lg_boundary_absorption&
                                      ,ig_ncpu&
                                      ,IG_LST_UNIT&
                                      ,IG_NDOF&
                                      ,info_all_cpu&
                                      ,LG_OUTPUT_DEBUG_FILE

      use mod_init_memory

      implicit none

      integer ,parameter :: NFACE=6
      integer ,parameter :: NEDGE=12
      integer ,parameter :: NNODE=8

      real               :: rl_max_coord_x
      real               :: rl_min_coord_x
      real               :: rl_max_coord_y
      real               :: rl_min_coord_y
      real               :: rl_max_coord_z
      real               :: rl_min_coord_z

      integer            :: ios
      integer            :: size_real_t
      integer            :: icpu
      integer            :: ihexa
      integer            :: iface
      integer            :: iedge
      integer            :: inode
      integer            :: ielt_type
      integer            :: ielt_number
      integer            :: ielt_face
      integer            :: ielt_edge
      integer            :: ielt_corner
      integer            :: ielt_coty
      integer            :: ielt_cpu
      integer            :: myunit_debug
      integer            :: myunit_dat
      integer            :: ifsurf
      integer            :: iparax
      character(len=255) :: fname
      character(len=255) :: info 
      
      integer            :: irec
      integer            :: igll
      integer            :: jgll
      integer            :: kgll
      integer            :: ineigh
      integer            :: nbneigh


      !
      !
      !********************************************************************************************************************************
      !open file that contains mesh information (outer and inner hexa, adjacency, etc.)
      !********************************************************************************************************************************
      fname = trim(cg_prefix)//"."//trim(cg_ncpu)//".elements.cpu."//trim(cg_myrank)//".dat"

      open(unit       = get_newunit(myunit_dat)&
          ,file       = trim(fname)&
          ,status     = 'old'&
          ,action     = 'read'&
          ,access     = 'sequential'&
          ,form       = 'unformatted'&
          ,recordtype = 'stream'&
          ,iostat     = ios)

      if (ios /= 0) then
         write(info,'(a)') "Error in subroutine init_mesh while opening hexa file"
         call error_stop(info)
      endif


      !
      !
      !********************************************************************************************************************************
      !read geometric nodes coordinates (right handed coordinate system: x-->est, y-->north, z-->upward)
      !********************************************************************************************************************************
      read(unit=myunit_dat) size_real_t,ig_mesh_nnode

      write(info,'(a)') "geometric nodes"
      call info_all_cpu(ig_mesh_nnode,info)

      ios = init_array_real(rg_gnode_x,ig_mesh_nnode,"rg_gnode_x")

      ios = init_array_real(rg_gnode_y,ig_mesh_nnode,"rg_gnode_y")

      ios = init_array_real(rg_gnode_z,ig_mesh_nnode,"rg_gnode_z")
      
      read(unit=myunit_dat) (rg_gnode_x(inode),rg_gnode_y(inode),rg_gnode_z(inode),inode=1,ig_mesh_nnode)


      !
      !
      !********************************************************************************************************************************
      !compute min/max (x,y,z) coordinates of the entire domain (i.e., for all cpu)
      !********************************************************************************************************************************
      rl_max_coord_x = maxval(rg_gnode_x(:))
      rl_max_coord_y = maxval(rg_gnode_y(:))
      rl_max_coord_z = maxval(rg_gnode_z(:))
      rl_min_coord_x = minval(rg_gnode_x(:))
      rl_min_coord_y = minval(rg_gnode_y(:))
      rl_min_coord_z = minval(rg_gnode_z(:))

      call mpi_allreduce(rl_max_coord_x,rg_mesh_xmax,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ios)
      call mpi_allreduce(rl_max_coord_y,rg_mesh_ymax,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ios)
      call mpi_allreduce(rl_max_coord_z,rg_mesh_zmax,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ios)
      call mpi_allreduce(rl_min_coord_x,rg_mesh_xmin,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ios)
      call mpi_allreduce(rl_min_coord_y,rg_mesh_ymin,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ios)
      call mpi_allreduce(rl_min_coord_z,rg_mesh_zmin,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ios)

      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(" "   ,/,a      )') "Boundaries of the entire domain"
         write(IG_LST_UNIT,'(" -->",  a,E14.7)') " xmin = ",rg_mesh_xmin
         write(IG_LST_UNIT,'(" -->",  a,E14.7)') " xmax = ",rg_mesh_xmax
         write(IG_LST_UNIT,'(" -->",  a,E14.7)') " ymin = ",rg_mesh_ymin
         write(IG_LST_UNIT,'(" -->",  a,E14.7)') " ymax = ",rg_mesh_ymax
         write(IG_LST_UNIT,'(" -->",  a,E14.7)') " zmin = ",rg_mesh_zmin
         write(IG_LST_UNIT,'(" -->",  a,E14.7)') " zmax = ",rg_mesh_zmax
      endif


      !
      !
      !********************************************************************************************************************************
      !read information about cpu connected to cpu myrank
      !********************************************************************************************************************************
      irec  = 1
      call read_bin_int(myunit_dat,irec,ig_ncpu_neighbor)
      write(info,'(a)') "connected cpu"
      call info_all_cpu(ig_ncpu_neighbor,info)
 
      allocate(tg_cpu_neighbor(ig_ncpu_neighbor),stat=ios)
      if (ios /= 0) then
         write(info,'(a)') "Error in subroutine init_mesh while allocating tg_cpu_neighbor"
         call error_stop(info)
      endif

      do icpu = 1,ig_ncpu_neighbor
         irec = irec + 1
         call read_bin_int(myunit_dat,irec,tg_cpu_neighbor(icpu)%icpu)
      enddo


      !
      !
      !********************************************************************************************************************************
      !read info about hexa and quad (number, number of geometric nodes per hexa, etc.)
      !********************************************************************************************************************************
      irec = irec + 1
      call read_bin_int(myunit_dat,irec,ig_nhexa)
      write(info,'(a)') " hexahedra"
      call info_all_cpu(ig_nhexa,info)
 
      irec = irec + 1
      call read_bin_int(myunit_dat,irec,ig_nhexa_outer)
      write(info,'(a)') " outer hexahedra"
      call info_all_cpu(ig_nhexa_outer,info)

      ig_nhexa_inner = ig_nhexa - ig_nhexa_outer
      write(info,'(a)') " inner hexahedra"
      call info_all_cpu(ig_nhexa_inner,info)
 
      irec = irec + 1
      call read_bin_int(myunit_dat,irec,ig_nquad_parax)
      write(info,'(a)') "paraxial quadrangles"
      call info_all_cpu(ig_nquad_parax,info)
 
      irec = irec + 1
      call read_bin_int(myunit_dat,irec,ig_nquad_fsurf)
      write(info,'(a)') "free surface quadrangles"
      call info_all_cpu(ig_nquad_fsurf,info)
 
      irec = irec + 1
      call read_bin_int(myunit_dat,irec,ig_hexa_nnode)
 
      if (ig_hexa_nnode == 8) then

         ig_quad_nnode = 4

      elseif (ig_hexa_nnode == 27) then

         ig_quad_nnode = 9

      else

         write(info,'(a)') "Error in subroutine init_mesh: illegal number of local nodes for hexa"
         call error_stop(info)

      endif


      !
      !
      !********************************************************************************************************************************
      !initialize convention numbering of finite element
      !********************************************************************************************************************************
      call init_element(ig_hexa_nnode,ig_quad_nnode)


      !
      !
      !********************************************************************************************************************************
      !initialize some arrays
      !********************************************************************************************************************************
      ios = init_array_int(ig_cpu_neighbor_info,26*ig_nhexa_outer,3,"ig_cpu_neighbor_info") !DAVID: deallocation in subroutine create_mpi_buffer()

      ios = init_array_int(ig_hexa_gnode_glonum,ig_nhexa,ig_hexa_nnode,"ig_hexa_gnode_glonum")

      ios = init_array_int(ig_hexa_gll_glonum,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"ig_hexa_gll_glonum")

      if (ig_nquad_parax > 0 .and. lg_boundary_absorption) then

         ios = init_array_int(ig_quadp_gnode_glonum,ig_nquad_parax,ig_quad_nnode,"ig_quadp_gnode_glonum")

         ios = init_array_int(ig_quadp_gll_glonum,ig_nquad_parax,IG_NGLL,IG_NGLL,"ig_quadp_gll_glonum")

         ios = init_array_int(ig_quadp_neighbor_hexa,ig_nquad_parax,"ig_quadp_neighbor_hexa")

         ios = init_array_int(ig_quadp_neighbor_hexaface,ig_nquad_parax,"ig_quadp_neighbor_hexaface")

      endif
 
      if (ig_nquad_fsurf > 0) then

         ios = init_array_int(ig_quadf_gnode_glonum,ig_nquad_fsurf,ig_quad_nnode,"ig_quadf_gnode_glonum")
         ios = init_array_int(ig_quadf_gll_glonum,ig_nquad_fsurf,IG_NGLL,IG_NGLL,"ig_quadf_gll_glonum")

      endif

      ios = init_array_int(ig_hexa_material_number,ig_nhexa,"ig_hexa_material_number")


      !
      !
      !********************************************************************************************************************************
      !initialize mesh and create GLL nodes numbering
      !********************************************************************************************************************************
      if (ig_myrank == 0) write(IG_LST_UNIT,'(" ",/,a)') "Generating global gll nodes numbering..."

      iparax = 0
      ifsurf = 0
      do ihexa = 1,ig_nhexa
         !
         !read material number of hexa 'ihexa'
         irec = irec + 1
         call read_bin_int(myunit_dat,irec,ig_hexa_material_number(ihexa))
         !
         !read geometric nodes of hexa 'ihexa'
         do inode = 1,ig_hexa_nnode
            irec = irec + 1
            call read_bin_int(myunit_dat,irec,ig_hexa_gnode_glonum(inode,ihexa))
         enddo
         !
         !generate gll node for current hexa
         call init_gll_number(ihexa,ig_ngll_total)
         !
         !loop on 6 faces to propagate gll node of hexa 'ihexa' to its neighbours
         do iface = 1,NFACE

            irec = irec + 1
            call read_bin_int(myunit_dat,irec,ielt_type) !either hexa or quad_parax or quad_fsurf

            select case (ielt_type)
               case(1) !another hexa is connected

                  irec = irec + 1
                  call read_bin_int(myunit_dat,irec,ielt_number)
                  if (ielt_number == 0) then
                     write(info,'(a)') "Error in subroutine init_mesh: ielt_number == 0 for face"
                     call error_stop(info)
                  endif

                  irec = irec + 1
                  call read_bin_int(myunit_dat,irec,ielt_face)
                  if (ielt_face == 0) then
                     write(info,'(a)') "Error in subroutine init_mesh: ielt_face == 0 for face"
                     call error_stop(info)
                  endif

                  irec = irec + 1
                  call read_bin_int(myunit_dat,irec,ielt_coty)

                  irec = irec + 1
                  call read_bin_int(myunit_dat,irec,ielt_cpu)

                  if (ielt_cpu == ig_myrank) then

                     if (ielt_number > ihexa) then
                        call propagate_gll_nodes_face(ihexa,iface,ielt_number,ielt_face,ielt_coty)
                     endif

                  else

                     ig_nneighbor_all_kind = ig_nneighbor_all_kind + 1

                     if (ig_nneighbor_all_kind > 26*ig_nhexa_outer) then
                        write(info,'(a)') "Error in subroutine init_mesh: ig_nneighbor_all_kind too large for face"
                        call error_stop(info)
                     endif

                     ig_cpu_neighbor_info(1, ig_nneighbor_all_kind) = ihexa
                     ig_cpu_neighbor_info(2, ig_nneighbor_all_kind) = iface
                     ig_cpu_neighbor_info(3, ig_nneighbor_all_kind) = ielt_cpu

                  endif

               case(2) !quad_parax is connected
                  if (ig_nquad_parax > 0 .and. lg_boundary_absorption) then

                     iparax                        = iparax + 1
                     ig_quadp_neighbor_hexa(iparax)      = ihexa
                     ig_quadp_neighbor_hexaface(iparax) = iface

                     if (ig_quad_nnode == 4) then
                        select case (iface)
                           case(1)
                              ig_quadp_gnode_glonum(1,iparax) = ig_hexa_gnode_glonum(1,ihexa)
                              ig_quadp_gnode_glonum(2,iparax) = ig_hexa_gnode_glonum(2,ihexa)
                              ig_quadp_gnode_glonum(3,iparax) = ig_hexa_gnode_glonum(3,ihexa)
                              ig_quadp_gnode_glonum(4,iparax) = ig_hexa_gnode_glonum(4,ihexa)
                           case(2)
                              ig_quadp_gnode_glonum(1,iparax) = ig_hexa_gnode_glonum(2,ihexa)
                              ig_quadp_gnode_glonum(2,iparax) = ig_hexa_gnode_glonum(1,ihexa)
                              ig_quadp_gnode_glonum(3,iparax) = ig_hexa_gnode_glonum(5,ihexa)
                              ig_quadp_gnode_glonum(4,iparax) = ig_hexa_gnode_glonum(6,ihexa)
                           case(3)
                              ig_quadp_gnode_glonum(1,iparax) = ig_hexa_gnode_glonum(3,ihexa)
                              ig_quadp_gnode_glonum(2,iparax) = ig_hexa_gnode_glonum(2,ihexa)
                              ig_quadp_gnode_glonum(3,iparax) = ig_hexa_gnode_glonum(6,ihexa)
                              ig_quadp_gnode_glonum(4,iparax) = ig_hexa_gnode_glonum(7,ihexa)
                           case(4)
                              ig_quadp_gnode_glonum(1,iparax) = ig_hexa_gnode_glonum(4,ihexa)
                              ig_quadp_gnode_glonum(2,iparax) = ig_hexa_gnode_glonum(3,ihexa)
                              ig_quadp_gnode_glonum(3,iparax) = ig_hexa_gnode_glonum(7,ihexa)
                              ig_quadp_gnode_glonum(4,iparax) = ig_hexa_gnode_glonum(8,ihexa)
                           case(5)
                              ig_quadp_gnode_glonum(1,iparax) = ig_hexa_gnode_glonum(1,ihexa)
                              ig_quadp_gnode_glonum(2,iparax) = ig_hexa_gnode_glonum(4,ihexa)
                              ig_quadp_gnode_glonum(3,iparax) = ig_hexa_gnode_glonum(8,ihexa)
                              ig_quadp_gnode_glonum(4,iparax) = ig_hexa_gnode_glonum(5,ihexa)
                           case(6)
                              ig_quadp_gnode_glonum(1,iparax) = ig_hexa_gnode_glonum(8,ihexa)
                              ig_quadp_gnode_glonum(2,iparax) = ig_hexa_gnode_glonum(7,ihexa)
                              ig_quadp_gnode_glonum(3,iparax) = ig_hexa_gnode_glonum(6,ihexa)
                              ig_quadp_gnode_glonum(4,iparax) = ig_hexa_gnode_glonum(5,ihexa)
                        end select
                     else
                        ! TODO : idem for 9 nodes quads
                        write(info,'(a)') "Error in init_mesh(): not yet ready for 27 nodes elements"
                        call error_stop(info)
                     endif
                     call propagate_gll_nodes_quad(ihexa,iface,iparax,ig_quadp_gll_glonum,ig_nquad_parax)
                  endif

               case(3) !quad_fsurf is connected

                  ifsurf = ifsurf + 1

                  if (ig_quad_nnode == 4) then
                     select case (iface)
                        case(1)
                           ig_quadf_gnode_glonum(1,ifsurf) = ig_hexa_gnode_glonum(1,ihexa) !see sub. init_element
                           ig_quadf_gnode_glonum(2,ifsurf) = ig_hexa_gnode_glonum(2,ihexa)
                           ig_quadf_gnode_glonum(3,ifsurf) = ig_hexa_gnode_glonum(3,ihexa)
                           ig_quadf_gnode_glonum(4,ifsurf) = ig_hexa_gnode_glonum(4,ihexa)
                        case(2)
                           ig_quadf_gnode_glonum(1,ifsurf) = ig_hexa_gnode_glonum(1,ihexa)
                           ig_quadf_gnode_glonum(2,ifsurf) = ig_hexa_gnode_glonum(5,ihexa)
                           ig_quadf_gnode_glonum(3,ifsurf) = ig_hexa_gnode_glonum(6,ihexa)
                           ig_quadf_gnode_glonum(4,ifsurf) = ig_hexa_gnode_glonum(2,ihexa)
                        case(3)
                           ig_quadf_gnode_glonum(1,ifsurf) = ig_hexa_gnode_glonum(2,ihexa)
                           ig_quadf_gnode_glonum(2,ifsurf) = ig_hexa_gnode_glonum(6,ihexa)
                           ig_quadf_gnode_glonum(3,ifsurf) = ig_hexa_gnode_glonum(7,ihexa)
                           ig_quadf_gnode_glonum(4,ifsurf) = ig_hexa_gnode_glonum(3,ihexa)
                        case(4)
                           ig_quadf_gnode_glonum(1,ifsurf) = ig_hexa_gnode_glonum(4,ihexa)
                           ig_quadf_gnode_glonum(2,ifsurf) = ig_hexa_gnode_glonum(3,ihexa)
                           ig_quadf_gnode_glonum(3,ifsurf) = ig_hexa_gnode_glonum(7,ihexa)
                           ig_quadf_gnode_glonum(4,ifsurf) = ig_hexa_gnode_glonum(8,ihexa)
                        case(5)
                           ig_quadf_gnode_glonum(1,ifsurf) = ig_hexa_gnode_glonum(1,ihexa)
                           ig_quadf_gnode_glonum(2,ifsurf) = ig_hexa_gnode_glonum(4,ihexa)
                           ig_quadf_gnode_glonum(3,ifsurf) = ig_hexa_gnode_glonum(8,ihexa)
                           ig_quadf_gnode_glonum(4,ifsurf) = ig_hexa_gnode_glonum(5,ihexa)
                        case(6)
                           ig_quadf_gnode_glonum(1,ifsurf) = ig_hexa_gnode_glonum(5,ihexa)
                           ig_quadf_gnode_glonum(2,ifsurf) = ig_hexa_gnode_glonum(8,ihexa)
                           ig_quadf_gnode_glonum(3,ifsurf) = ig_hexa_gnode_glonum(7,ihexa)
                           ig_quadf_gnode_glonum(4,ifsurf) = ig_hexa_gnode_glonum(6,ihexa)
                     end select
                  else
                     ! TODO : idem for 9 nodes quads
                     write(info,'(a)') "Error in init_mesh(): not yet ready for 27 nodes elements"
                     call error_stop(info)
                  endif
                  call propagate_gll_nodes_quad(ihexa,iface,ifsurf,ig_quadf_gll_glonum,ig_nquad_fsurf)

               case default !no element is connected

            end select
         enddo
 
         do iedge = 1,NEDGE

            irec = irec + 1
            call read_bin_int(myunit_dat,irec,nbneigh) !nb neighbors of iedge
            
            do ineigh = 1, nbneigh

               irec = irec + 1
               call read_bin_int(myunit_dat,irec,ielt_type) !hexa or none

               select case (ielt_type)
                  case(1) !another hexa is connected

                     irec = irec + 1
                     call read_bin_int(myunit_dat,irec,ielt_number)

                     irec = irec + 1
                     call read_bin_int(myunit_dat,irec,ielt_edge)

                     irec = irec + 1
                     call read_bin_int(myunit_dat,irec,ielt_coty)

                     irec = irec + 1
                     call read_bin_int(myunit_dat,irec,ielt_cpu)

                     if (ielt_cpu == ig_myrank) then

                        if (ielt_number > ihexa) then
                           call propagate_gll_nodes_edge(ihexa,iedge,ielt_number,ielt_edge,ielt_coty)
                        endif

                     else

                        ig_nneighbor_all_kind = ig_nneighbor_all_kind + 1

                        if (ig_nneighbor_all_kind > 26 * ig_nhexa_outer) then
                           write(info,'(a)') "Error in subroutine init_mesh: ig_nneighbor_all_kind too large for edge"
                           call error_stop(info)
                        endif

                        ig_cpu_neighbor_info(1, ig_nneighbor_all_kind) = ihexa
                        ig_cpu_neighbor_info(2, ig_nneighbor_all_kind) = NFACE + iedge
                        ig_cpu_neighbor_info(3, ig_nneighbor_all_kind) = ielt_cpu

                     endif

                  case default !no element is connected

               end select

            enddo
         enddo

         do inode = 1,NNODE

            irec = irec + 1
            call read_bin_int(myunit_dat,irec,nbneigh) ! nb neighbors of iedge
            
            do ineigh = 1, nbneigh

               irec = irec + 1
               call read_bin_int(myunit_dat,irec,ielt_type) !hexa or none

               select case (ielt_type)

                  case(1) !another hexa is connected

                     irec = irec + 1
                     call read_bin_int(myunit_dat,irec,ielt_number)

                     irec = irec + 1
                     call read_bin_int(myunit_dat,irec,ielt_corner)

                     irec = irec + 1
                     call read_bin_int(myunit_dat,irec,ielt_cpu)

                     if (ielt_cpu == ig_myrank) then

                        if (ielt_number > ihexa) then
                           call propagate_gll_nodes_corner(ihexa,inode,ielt_number,ielt_corner)
                        endif

                     else

                        ig_nneighbor_all_kind = ig_nneighbor_all_kind + 1

                        if (ig_nneighbor_all_kind > 26*ig_nhexa_outer) then
                           write(info,'(a)') "Error in subroutine init_mesh: ig_nneighbor_all_kind too large for node"
                           call error_stop(info)
                        endif

                        ig_cpu_neighbor_info(1, ig_nneighbor_all_kind) = ihexa
                        ig_cpu_neighbor_info(2, ig_nneighbor_all_kind) = NFACE + NEDGE + inode
                        ig_cpu_neighbor_info(3, ig_nneighbor_all_kind) = ielt_cpu

                     endif

                  case default !no element is connected

               end select

            enddo
         enddo
      enddo
      close(myunit_dat)
      if (ig_myrank == 0) write(IG_LST_UNIT,'(a)') "Done"


      !
      !
      !***********************************************************
      !initialize displacement, velocity and acceleration arrays
      !***********************************************************
      ios = init_array_real(rg_gll_displacement,ig_ngll_total,IG_NDOF,"rg_gll_displacement")

      ios = init_array_real(rg_gll_velocity    ,ig_ngll_total,IG_NDOF,"rg_gll_velocity")

      ios = init_array_real(rg_gll_acceleration,ig_ngll_total,IG_NDOF,"rg_gll_acceleration")

      ios = init_array_real(rg_gll_acctmp      ,ig_ngll_total,IG_NDOF,"rg_gll_acctmp")


      !
      !
      !******************************
      !write gll info in *.lst file
      !******************************
      write(info,'(a)') "gll nodes"
      call info_all_cpu(ig_ngll_total,info)


      !
      !
      !*************************************************
      !debug mode : write array ig_hexa_gll_glonum
      !*************************************************
      if (LG_OUTPUT_DEBUG_FILE) then
         open(unit=get_newunit(myunit_debug),file="debug."//trim(cg_prefix)//".global.gll."//trim(cg_myrank))
         do ihexa = 1,ig_nhexa
            write(unit=myunit_debug,fmt='(a,i10)') "hexa ",ihexa
            do igll = 1,IG_NGLL
               write(unit=myunit_debug,fmt='(a,i10)') "igll ",igll
               do jgll = 1,IG_NGLL
                  write(unit=myunit_debug,fmt='(10I10)') (ig_hexa_gll_glonum(kgll,jgll,igll,ihexa),kgll=1,IG_NGLL)
               enddo
            enddo
         enddo
         close(myunit_debug)
      endif

      !
      !
      !*************************************************
      !debug mode : write free surface geom nodes
      !*************************************************
      if (LG_OUTPUT_DEBUG_FILE) then
         open(unit=get_newunit(myunit_debug),file="debug."//trim(cg_prefix)//".freesurface.geom."//trim(cg_myrank))
         do ifsurf = 1,ig_nquad_fsurf
            do inode = 1,ig_quad_nnode
               write(unit=myunit_debug,fmt='(2(E14.7,1X))') rg_gnode_x(ig_quadf_gnode_glonum(inode,ifsurf))&
                                                           ,rg_gnode_y(ig_quadf_gnode_glonum(inode,ifsurf))
            enddo
         enddo
         close(myunit_debug)
      endif

      return
!***********************************************************************************************************************************************************************************
   end subroutine init_mesh
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine increments GLL numbering of hexahedron elements in cpu myrank (see variable mod_global_variables::ig_ngll_total) 
!!and creates the local to global GLL indirection array mod_global_variables::ig_hexa_gll_glonum for hexahedron elements.
!>@param ihexa    : hexahedron element number in cpu myrank
!>@param ngll_cpu : local name of global variable mod_global_variables::ig_ngll_total
!***********************************************************************************************************************************************************************************
   subroutine init_gll_number(ihexa,ngll_cpu)
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only : &
                                       ig_hexa_gll_glonum&
                                      ,IG_NGLL
       
      implicit none

      integer, intent(in)                        :: ihexa
      integer, intent(inout)                     :: ngll_cpu
                                                 
      integer                                    :: k
      integer                                    :: l
      integer                                    :: m

      do k = 1,IG_NGLL        !along zeta
         do l = 1,IG_NGLL     !along eta
            do m = 1,IG_NGLL  !along xi

               if (ig_hexa_gll_glonum(m,l,k,ihexa) == 0) then !GLL has not been numbered yet

                  ngll_cpu                           = ngll_cpu + 1
                  ig_hexa_gll_glonum(m,l,k,ihexa) = ngll_cpu

               endif

            enddo
         enddo
      enddo

      return
!***********************************************************************************************************************************************************************************
   end subroutine init_gll_number
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine creates the local to global GLL indirection array for quadrangle elements by propagating existing GLL numbering of hexahedron elements.
!>@param ihexa              : hexahedron element number to which iquad is connected
!>@param iface              : face number of element ihexa to which iquad is connected
!>@param iquad              : quadrangle element connected to hexahedron element ihexa
!>@param global_gll_of_quad : local to global GLL indirection array for quadrangle elements
!>@param number_of_quad     : number of quadrangle elements in array global_gll_of_quad
!***********************************************************************************************************************************************************************************
    subroutine propagate_gll_nodes_quad(ihexa,iface,iquad,global_gll_of_quad,number_of_quad)
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only : IG_NGLL&
                                      ,error_stop&
                                      ,ig_hexa_gll_glonum

      implicit none

      integer, intent(in)                                               :: ihexa
      integer, intent(in)                                               :: iface
      integer, intent(in)                                               :: iquad
      integer, intent(in)                                               :: number_of_quad
      integer, intent(inout), dimension(IG_NGLL,IG_NGLL,number_of_quad) :: global_gll_of_quad
      
      integer                                                           :: i,j,k
 
      character(len=255)                                                :: info
      
      select case(iface)
         case(1)
            do j = 1,IG_NGLL
               do i = 1,IG_NGLL
                  global_gll_of_quad(i,j,iquad) = ig_hexa_gll_glonum(i,j,1,ihexa)
               enddo
            enddo
            return

         case(2)
            do k = 1,IG_NGLL
               do i = 1,IG_NGLL
                  global_gll_of_quad(k,i,iquad) = ig_hexa_gll_glonum(i,1,k,ihexa)
               enddo
            enddo
            return

         case(3)
            do k = 1,IG_NGLL
               do j = 1,IG_NGLL
                  global_gll_of_quad(k,j,iquad) = ig_hexa_gll_glonum(IG_NGLL,j,k,ihexa)
               enddo
            enddo
            return

         case(4)
            do k = 1,IG_NGLL
               do i = 1,IG_NGLL
                  global_gll_of_quad(i,k,iquad) = ig_hexa_gll_glonum(i,IG_NGLL,k,ihexa)
               enddo
            enddo
            return

         case(5)
            do k = 1,IG_NGLL
               do j = 1,IG_NGLL
                  global_gll_of_quad(j,k,iquad) = ig_hexa_gll_glonum(1,j,k,ihexa)
               enddo
            enddo
            return

         case(6)
            do j = 1,IG_NGLL
               do i = 1,IG_NGLL
                  global_gll_of_quad(j,i,iquad) = ig_hexa_gll_glonum(i,j,IG_NGLL,ihexa)
               enddo
            enddo
            return

         case default
            write(info,'(a)') "error in subroutine propagate_gll_nodes_quad. Invalid face number"
            call error_stop(info)

      end select
      
!***********************************************************************************************************************************************************************************
   end subroutine propagate_gll_nodes_quad
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine propagates existing GLL numbering of hexahedron elements to the face of its neighbors hexahedron elements.
!>@param ihexa_old : hexahedron element already numbered to which ihexa_new is connected
!>@param iface_old : face number of element ihexa_old to which ihexa_new is connected
!>@param ihexa_new : hexahedron element to be numbered connected to ihexa_old
!>@param iface_new : face number of element ihexa_new connected to iface_old
!>@param icoty_new : connexion type to redirect correctly GLL numbering between connected faces of hexahedron elements
!***********************************************************************************************************************************************************************************
    subroutine propagate_gll_nodes_face(ihexa_old,iface_old,ihexa_new,iface_new,icoty_new)
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only : IG_NGLL&
                                      ,ig_hexa_gll_glonum

      implicit none

      integer, intent(in)    :: ihexa_new
      integer, intent(in)    :: ihexa_old
      integer, intent(in)    :: iface_new
      integer, intent(in)    :: iface_old
      integer, intent(in)    :: icoty_new
     
      integer                :: ext_vec
      integer                :: ext_sgn
      integer                :: int_vec
      integer                :: int_sgn
      integer                :: fixed_val
      integer                :: i_ext
      integer                :: i_int
      integer                :: i_ext_rev
      integer                :: i_int_rev
      integer                :: i_ext_target
      integer                :: i_int_target
      integer                :: igll_source
      integer                :: mid_gll

     
!     icoty_new : 
!     0 . 0 .SX.VX.VX.SI.VI.VI
!     7   6  5  4  3  2  1  0   th bit 
!     128 64 32 16 8  4  2  1  
      
      ext_vec = ishft(iand(24, icoty_new), -3)
      ext_sgn = ishft(iand(32, icoty_new), -5)
      
      int_vec = iand(3, icoty_new)
      int_sgn = ishft(iand(4, icoty_new), -2)
 
      mid_gll   = ceiling(IG_NGLL/2.0)
      fixed_val = 0
      select case(iface_new)
         case(1)
            if (ig_hexa_gll_glonum(mid_gll,mid_gll,1,ihexa_new) /= 0) return
            fixed_val = 1
         case(2)
            if (ig_hexa_gll_glonum(mid_gll,1,mid_gll,ihexa_new) /= 0) return
            fixed_val = 1
         case(3)
            if (ig_hexa_gll_glonum(IG_NGLL,mid_gll,mid_gll,ihexa_new) /= 0) return
            fixed_val = IG_NGLL
         case(4)
            if (ig_hexa_gll_glonum(mid_gll,IG_NGLL,mid_gll,ihexa_new) /= 0) return
            fixed_val = IG_NGLL
         case(5)
            if (ig_hexa_gll_glonum(1,mid_gll,mid_gll,ihexa_new) /= 0) return
            fixed_val = 1
         case(6)
            if (ig_hexa_gll_glonum(mid_gll,mid_gll,IG_NGLL,ihexa_new) /= 0) return
            fixed_val = IG_NGLL
      end select
   
      do i_ext = 1,IG_NGLL

         i_ext_rev = (IG_NGLL+1)-i_ext

         do i_int = 1,IG_NGLL
            i_int_rev = (IG_NGLL+1)-i_int
            select case(iface_old)
               case(1)
                  igll_source = ig_hexa_gll_glonum(i_ext, i_int, 1, ihexa_old)
               case(2)
                  igll_source = ig_hexa_gll_glonum(i_int, 1, i_ext, ihexa_old)
               case(3)
                  igll_source = ig_hexa_gll_glonum(IG_NGLL, i_int_rev, i_ext_rev, ihexa_old)
               case(4)
                  igll_source = ig_hexa_gll_glonum(i_ext_rev, IG_NGLL, i_int_rev, ihexa_old)
               case(5)
                  igll_source = ig_hexa_gll_glonum(1, i_ext, i_int, ihexa_old)
               case(6)
                  igll_source = ig_hexa_gll_glonum(i_int_rev, i_ext_rev, IG_NGLL, ihexa_old)
            end select
            
            if (ext_sgn == 1) then
               i_ext_target = i_ext
            else
               i_ext_target = i_ext_rev
            endif
            if (int_sgn == 1) then
               i_int_target = i_int
            else
               i_int_target = i_int_rev
            endif
            
            select case(ext_vec)
               case(1)
                  select case(int_vec)
                     case(2)
                        if (ig_hexa_gll_glonum(fixed_val, i_int_target, i_ext_target, ihexa_new) /= igll_source&
                         .and. ig_hexa_gll_glonum(fixed_val, i_int_target, i_ext_target, ihexa_new) /= 0) write(*,*) 'propagate_gll_nodes_face() ', ihexa_old,iface_old,ihexa_new,iface_new
                        ig_hexa_gll_glonum(fixed_val, i_int_target, i_ext_target, ihexa_new) = igll_source
                     case(3)
                        if (ig_hexa_gll_glonum(i_int_target, fixed_val, i_ext_target, ihexa_new) /= igll_source&
                          .and. ig_hexa_gll_glonum(i_int_target, fixed_val, i_ext_target, ihexa_new) /= 0) write(*,*) 'propagate_gll_nodes_face() ', ihexa_old,iface_old,ihexa_new,iface_new
                        ig_hexa_gll_glonum(i_int_target, fixed_val, i_ext_target, ihexa_new) = igll_source
                  end select
               case(2)
                  select case(int_vec)
                     case(1)
                        if (ig_hexa_gll_glonum(fixed_val, i_ext_target, i_int_target, ihexa_new) /= igll_source&
                          .and. ig_hexa_gll_glonum(fixed_val, i_ext_target, i_int_target, ihexa_new) /= 0) write(*,*) 'propagate_gll_nodes_face() ', ihexa_old,iface_old,ihexa_new,iface_new
                        ig_hexa_gll_glonum(fixed_val, i_ext_target, i_int_target, ihexa_new) = igll_source
                     case(3)
                        if (ig_hexa_gll_glonum(i_int_target, i_ext_target, fixed_val, ihexa_new) /= igll_source&
                          .and. ig_hexa_gll_glonum(i_int_target, i_ext_target, fixed_val, ihexa_new) /= 0) write(*,*) 'propagate_gll_nodes_face() ', ihexa_old,iface_old,ihexa_new,iface_new
                        ig_hexa_gll_glonum(i_int_target, i_ext_target, fixed_val, ihexa_new) = igll_source
                  end select
               case(3)
                  select case(int_vec)
                     case(1)
                        if (ig_hexa_gll_glonum(i_ext_target, fixed_val, i_int_target, ihexa_new) /= igll_source&
                          .and. ig_hexa_gll_glonum(i_ext_target, fixed_val, i_int_target, ihexa_new) /= 0) write(*,*) 'propagate_gll_nodes_face() ', ihexa_old,iface_old,ihexa_new,iface_new
                        ig_hexa_gll_glonum(i_ext_target, fixed_val, i_int_target, ihexa_new) = igll_source
                     case(2)
                        if (ig_hexa_gll_glonum(i_ext_target, i_int_target, fixed_val, ihexa_new) /= igll_source&
                          .and. ig_hexa_gll_glonum(i_ext_target, i_int_target, fixed_val, ihexa_new) /= 0) write(*,*) 'propagate_gll_nodes_face() ', ihexa_old,iface_old,ihexa_new,iface_new
                        ig_hexa_gll_glonum(i_ext_target, i_int_target, fixed_val, ihexa_new) = igll_source
                  end select
            end select
         enddo
      enddo
      
      return
!***********************************************************************************************************************************************************************************
   end subroutine propagate_gll_nodes_face
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine propagates existing GLL numbering of hexahedron elements to the edge of its neighbors hexahedron elements.
!>@param ihexa_old : hexahedron element already numbered to which ihexa_new is connected
!>@param iedge_old : edge number of element ihexa_old to which ihexa_new is connected
!>@param ihexa_new : hexahedron element to be numbered connected to ihexa_old
!>@param iedge_new : edge number of element ihexa_new connected to iedge_old
!>@param icoty_new : connexion type to redirect correctly GLL numbering between connected edges of hexahedron elements
!***********************************************************************************************************************************************************************************
   subroutine propagate_gll_nodes_edge(ihexa_old,iedge_old,ihexa_new,iedge_new,icoty_new)
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only : IG_NGLL&
                                      ,ig_myrank&
                                      ,ig_hexa_gll_glonum&
                                      ,LG_OUTPUT_DEBUG_FILE

      implicit none

      integer, intent(in)    :: ihexa_new
      integer, intent(in)    :: ihexa_old
      integer, intent(in)    :: iedge_new
      integer, intent(in)    :: iedge_old
      integer, intent(in)    :: icoty_new
      
      integer, pointer       :: pedge_old(:)
      integer, pointer       :: pedge_new(:)
      integer                :: i, j

      select case(iedge_old)
         case(1)
            pedge_old => ig_hexa_gll_glonum(1:IG_NGLL,1,1,ihexa_old)
         case(2)
            pedge_old => ig_hexa_gll_glonum(IG_NGLL,1:IG_NGLL,1,ihexa_old)
         case(3)
            pedge_old => ig_hexa_gll_glonum(1:IG_NGLL,IG_NGLL,1,ihexa_old)
         case(4)
            pedge_old => ig_hexa_gll_glonum(1,1:IG_NGLL,1,ihexa_old)
         case(5)
            pedge_old => ig_hexa_gll_glonum(1:IG_NGLL,1,IG_NGLL,ihexa_old)
         case(6)
            pedge_old => ig_hexa_gll_glonum(IG_NGLL,1:IG_NGLL,IG_NGLL,ihexa_old)
         case(7)
            pedge_old => ig_hexa_gll_glonum(1:IG_NGLL,IG_NGLL,IG_NGLL,ihexa_old)
         case(8)
            pedge_old => ig_hexa_gll_glonum(1,1:IG_NGLL,IG_NGLL,ihexa_old)
         case(9)
            pedge_old => ig_hexa_gll_glonum(1,1,1:IG_NGLL,ihexa_old)
         case(10)
            pedge_old => ig_hexa_gll_glonum(IG_NGLL,1,1:IG_NGLL,ihexa_old)
         case(11)
            pedge_old => ig_hexa_gll_glonum(IG_NGLL,IG_NGLL,1:IG_NGLL,ihexa_old)
         case(12)
            pedge_old => ig_hexa_gll_glonum(1,IG_NGLL,1:IG_NGLL,ihexa_old)
      end select

      select case(iedge_new)
         case(1)
            pedge_new => ig_hexa_gll_glonum(1:IG_NGLL,1,1,ihexa_new)
         case(2)
            pedge_new => ig_hexa_gll_glonum(IG_NGLL,1:IG_NGLL,1,ihexa_new)
         case(3)
            pedge_new => ig_hexa_gll_glonum(1:IG_NGLL,IG_NGLL,1,ihexa_new)
         case(4)
            pedge_new => ig_hexa_gll_glonum(1,1:IG_NGLL,1,ihexa_new)
         case(5)
            pedge_new => ig_hexa_gll_glonum(1:IG_NGLL,1,IG_NGLL,ihexa_new)
         case(6)
            pedge_new => ig_hexa_gll_glonum(IG_NGLL,1:IG_NGLL,IG_NGLL,ihexa_new)
         case(7)
            pedge_new => ig_hexa_gll_glonum(1:IG_NGLL,IG_NGLL,IG_NGLL,ihexa_new)
         case(8)
            pedge_new => ig_hexa_gll_glonum(1,1:IG_NGLL,IG_NGLL,ihexa_new)
         case(9)
            pedge_new => ig_hexa_gll_glonum(1,1,1:IG_NGLL,ihexa_new)
         case(10)
            pedge_new => ig_hexa_gll_glonum(IG_NGLL,1,1:IG_NGLL,ihexa_new)
         case(11)
            pedge_new => ig_hexa_gll_glonum(IG_NGLL,IG_NGLL,1:IG_NGLL,ihexa_new)
         case(12)
            pedge_new => ig_hexa_gll_glonum(1,IG_NGLL,1:IG_NGLL,ihexa_new)
      end select
      
      do i = 1,IG_NGLL
         if (icoty_new == 1) then
            j = i
         else 
            j = IG_NGLL+1-i
         endif
         if ( (pedge_new(j) /= pedge_old(i) .and. pedge_new(j) /= 0) .or. pedge_old(i) == 0 ) write(*,*) 'propagate_gll_nodes_edge',ihexa_old,iedge_old,ihexa_new,iedge_new
         pedge_new(j) = pedge_old(i)
      enddo
      
      return
!***********************************************************************************************************************************************************************************
   end subroutine propagate_gll_nodes_edge
!***********************************************************************************************************************************************************************************

!
! 
!>@brief This subroutine propagates existing GLL numbering of hexahedron elements to corner (i.e., vertex) of its neighbors hexahedron elements.
!>@param ihexa_old   : hexahedron element already numbered to which ihexa_new is connected
!>@param icorner_old : corner number of element ihexa_old to which ihexa_new is connected
!>@param ihexa_new   : hexahedron element to be numbered connected to ihexa_old
!>@param icorner_new : corner number of element ihexa_new connected to icorner_old
!***********************************************************************************************************************************************************************************
   subroutine propagate_gll_nodes_corner(ihexa_old,icorner_old,ihexa_new,icorner_new)
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only : IG_NGLL&
                                      ,ig_myrank&
                                      ,ig_hexa_gll_glonum&
                                      ,LG_OUTPUT_DEBUG_FILE

      implicit none

      integer, intent(in)    :: ihexa_new
      integer, intent(in)    :: ihexa_old
      integer, intent(in)    :: icorner_new
      integer, intent(in)    :: icorner_old
      
      integer, pointer       :: pcorner_old
      integer, pointer       :: pcorner_new

      select case(icorner_old)
         case(1)
            pcorner_old => ig_hexa_gll_glonum(1,1,1,ihexa_old)
         case(2)
            pcorner_old => ig_hexa_gll_glonum(IG_NGLL,1,1,ihexa_old)
         case(3)
            pcorner_old => ig_hexa_gll_glonum(IG_NGLL,IG_NGLL,1,ihexa_old)
         case(4)
            pcorner_old => ig_hexa_gll_glonum(1,IG_NGLL,1,ihexa_old)
         case(5)
            pcorner_old => ig_hexa_gll_glonum(1,1,IG_NGLL,ihexa_old)
         case(6)
            pcorner_old => ig_hexa_gll_glonum(IG_NGLL,1,IG_NGLL,ihexa_old)
         case(7)
            pcorner_old => ig_hexa_gll_glonum(IG_NGLL,IG_NGLL,IG_NGLL,ihexa_old)
         case(8)
            pcorner_old => ig_hexa_gll_glonum(1,IG_NGLL,IG_NGLL,ihexa_old)
      end select

      select case(icorner_new)
         case(1)
            pcorner_new => ig_hexa_gll_glonum(1,1,1,ihexa_new)
         case(2)
            pcorner_new => ig_hexa_gll_glonum(IG_NGLL,1,1,ihexa_new)
         case(3)
            pcorner_new => ig_hexa_gll_glonum(IG_NGLL,IG_NGLL,1,ihexa_new)
         case(4)
            pcorner_new => ig_hexa_gll_glonum(1,IG_NGLL,1,ihexa_new)
         case(5)
            pcorner_new => ig_hexa_gll_glonum(1,1,IG_NGLL,ihexa_new)
         case(6)
            pcorner_new => ig_hexa_gll_glonum(IG_NGLL,1,IG_NGLL,ihexa_new)
         case(7)
            pcorner_new => ig_hexa_gll_glonum(IG_NGLL,IG_NGLL,IG_NGLL,ihexa_new)
         case(8)
            pcorner_new => ig_hexa_gll_glonum(1,IG_NGLL,IG_NGLL,ihexa_new)
      end select
      if ( (pcorner_new /= pcorner_old .and. pcorner_new /= 0) .or. pcorner_old == 0 ) write(*,*) 'propagate_gll_nodes_corner',ihexa_old,icorner_old,ihexa_new,icorner_new
      pcorner_new = pcorner_old
      
      return
!***********************************************************************************************************************************************************************************
   end subroutine propagate_gll_nodes_corner
!***********************************************************************************************************************************************************************************

!
!
!>@brief subroutine to set up convention of hexahedron and quadrangle elements
!>@param ilnnhe : local name of global variable mod_global_variables::ig_hexa_nnode
!>@param ilnnqu : local name of global variable mod_global_variables::ig_quad_nnode
!***********************************************************************************************************************************************************************************
    subroutine init_element(ilnnhe,ilnnqu)
!***********************************************************************************************************************************************************************************
 
      use mpi

      use mod_global_variables, only : ig_hexa_gnode_xiloc&
                                      ,ig_hexa_gnode_etloc&
                                      ,ig_hexa_gnode_zeloc&
                                      ,ig_line_nnode&
                                      ,ig_hexa_node2gll&
                                      ,rg_gnode_abscissa&
                                      ,ig_quad_gnode_xiloc&
                                      ,ig_quad_gnode_etloc&
                                      ,rg_gnode_abscissa_dist&
                                      ,IG_NGLL&
                                      ,error_stop
      
      implicit none
      
      integer, intent(in) :: ilnnhe
      integer, intent(in) :: ilnnqu

      integer             :: i
      integer             :: j

      character(len=255)  :: info
      !
      !
      !**********************************************************************************
      !fill local position of geometric node of hexa8 or hexa27 on xi, eta and zeta axis
      !**********************************************************************************
      allocate(ig_hexa_gnode_xiloc(ilnnhe),ig_hexa_gnode_etloc(ilnnhe),ig_hexa_gnode_zeloc(ilnnhe))

      if (ilnnhe == 8) then

         ig_line_nnode = 2
!
!------->fill rg_gnode_abscissa: local coordinate of geometric nodes
         allocate(rg_gnode_abscissa(ig_line_nnode))
         rg_gnode_abscissa(1) = +1.0
         rg_gnode_abscissa(2) = -1.0
!
!------->local position of geometric nodes
         ig_hexa_gnode_xiloc(1) = 1
         ig_hexa_gnode_etloc(1) = 1
         ig_hexa_gnode_zeloc(1) = 1
         ig_hexa_gnode_xiloc(2) = 2
         ig_hexa_gnode_etloc(2) = 1
         ig_hexa_gnode_zeloc(2) = 1
         ig_hexa_gnode_xiloc(3) = 2
         ig_hexa_gnode_etloc(3) = 2
         ig_hexa_gnode_zeloc(3) = 1
         ig_hexa_gnode_xiloc(4) = 1
         ig_hexa_gnode_etloc(4) = 2
         ig_hexa_gnode_zeloc(4) = 1
         ig_hexa_gnode_xiloc(5) = 1
         ig_hexa_gnode_etloc(5) = 1
         ig_hexa_gnode_zeloc(5) = 2
         ig_hexa_gnode_xiloc(6) = 2
         ig_hexa_gnode_etloc(6) = 1
         ig_hexa_gnode_zeloc(6) = 2
         ig_hexa_gnode_xiloc(7) = 2
         ig_hexa_gnode_etloc(7) = 2
         ig_hexa_gnode_zeloc(7) = 2
         ig_hexa_gnode_xiloc(8) = 1
         ig_hexa_gnode_etloc(8) = 2
         ig_hexa_gnode_zeloc(8) = 2

      elseif (ilnnhe == 27) then

         ig_line_nnode = 3
!
!------->fill rg_gnode_abscissa: local coordinate of geometric nodes
         allocate(rg_gnode_abscissa(ig_line_nnode))
         rg_gnode_abscissa(1) = +1.0
         rg_gnode_abscissa(2) =  0.0
         rg_gnode_abscissa(3) = -1.0
!
!------->local position of geometric nodes
         ig_hexa_gnode_xiloc( 1) = 1
         ig_hexa_gnode_etloc( 1) = 1
         ig_hexa_gnode_zeloc( 1) = 1
         ig_hexa_gnode_xiloc( 2) = 3
         ig_hexa_gnode_etloc( 2) = 1
         ig_hexa_gnode_zeloc( 2) = 1
         ig_hexa_gnode_xiloc( 3) = 3
         ig_hexa_gnode_etloc( 3) = 3
         ig_hexa_gnode_zeloc( 3) = 1
         ig_hexa_gnode_xiloc( 4) = 1
         ig_hexa_gnode_etloc( 4) = 3
         ig_hexa_gnode_zeloc( 4) = 1
         ig_hexa_gnode_xiloc( 5) = 1
         ig_hexa_gnode_etloc( 5) = 1
         ig_hexa_gnode_zeloc( 5) = 3
         ig_hexa_gnode_xiloc( 6) = 3
         ig_hexa_gnode_etloc( 6) = 1
         ig_hexa_gnode_zeloc( 6) = 3
         ig_hexa_gnode_xiloc( 7) = 3
         ig_hexa_gnode_etloc( 7) = 3
         ig_hexa_gnode_zeloc( 7) = 3
         ig_hexa_gnode_xiloc( 8) = 1
         ig_hexa_gnode_etloc( 8) = 3
         ig_hexa_gnode_zeloc( 8) = 3
         ig_hexa_gnode_xiloc( 9) = 2
         ig_hexa_gnode_etloc( 9) = 1
         ig_hexa_gnode_zeloc( 9) = 1
         ig_hexa_gnode_xiloc(10) = 3
         ig_hexa_gnode_etloc(10) = 2
         ig_hexa_gnode_zeloc(10) = 1
         ig_hexa_gnode_xiloc(11) = 2
         ig_hexa_gnode_etloc(11) = 3
         ig_hexa_gnode_zeloc(11) = 1
         ig_hexa_gnode_xiloc(12) = 1
         ig_hexa_gnode_etloc(12) = 2
         ig_hexa_gnode_zeloc(12) = 1
         ig_hexa_gnode_xiloc(13) = 2
         ig_hexa_gnode_etloc(13) = 1
         ig_hexa_gnode_zeloc(13) = 3
         ig_hexa_gnode_xiloc(14) = 3
         ig_hexa_gnode_etloc(14) = 2
         ig_hexa_gnode_zeloc(14) = 3
         ig_hexa_gnode_xiloc(15) = 2
         ig_hexa_gnode_etloc(15) = 3
         ig_hexa_gnode_zeloc(15) = 3
         ig_hexa_gnode_xiloc(16) = 1
         ig_hexa_gnode_etloc(16) = 2
         ig_hexa_gnode_zeloc(16) = 3
         ig_hexa_gnode_xiloc(17) = 1
         ig_hexa_gnode_etloc(17) = 1
         ig_hexa_gnode_zeloc(17) = 2
         ig_hexa_gnode_xiloc(18) = 3
         ig_hexa_gnode_etloc(18) = 1
         ig_hexa_gnode_zeloc(18) = 2
         ig_hexa_gnode_xiloc(19) = 3
         ig_hexa_gnode_etloc(19) = 3
         ig_hexa_gnode_zeloc(19) = 2
         ig_hexa_gnode_xiloc(20) = 1
         ig_hexa_gnode_etloc(20) = 3
         ig_hexa_gnode_zeloc(20) = 2
         ig_hexa_gnode_xiloc(21) = 2
         ig_hexa_gnode_etloc(21) = 2 
         ig_hexa_gnode_zeloc(21) = 1 
         ig_hexa_gnode_xiloc(22) = 2 
         ig_hexa_gnode_etloc(22) = 1 
         ig_hexa_gnode_zeloc(22) = 2 
         ig_hexa_gnode_xiloc(23) = 3 
         ig_hexa_gnode_etloc(23) = 2 
         ig_hexa_gnode_zeloc(23) = 2 
         ig_hexa_gnode_xiloc(24) = 2 
         ig_hexa_gnode_etloc(24) = 3 
         ig_hexa_gnode_zeloc(24) = 2 
         ig_hexa_gnode_xiloc(25) = 1 
         ig_hexa_gnode_etloc(25) = 2 
         ig_hexa_gnode_zeloc(25) = 2 
         ig_hexa_gnode_xiloc(26) = 2 
         ig_hexa_gnode_etloc(26) = 2 
         ig_hexa_gnode_zeloc(26) = 3 
         ig_hexa_gnode_xiloc(27) = 2 
         ig_hexa_gnode_etloc(27) = 2 
         ig_hexa_gnode_zeloc(27) = 2 
      else
         write(info,'(a)') "error in init_element: invalid number of geometrical nodes for hexa"
         call error_stop(info)
      endif
      allocate(rg_gnode_abscissa_dist(ig_line_nnode,ig_line_nnode))
      do i = 1,ig_line_nnode
         do j = 1,ig_line_nnode
            rg_gnode_abscissa_dist(j,i) = 0.0
            if (i /= j) rg_gnode_abscissa_dist(j,i) = 1.0/(rg_gnode_abscissa(i) - rg_gnode_abscissa(j))
         enddo
      enddo
      !!******************************************************************
      !fill local position of geometric node of quad4 or quad9 on xi, eta
      !!******************************************************************
      !rg_gnode_abscissa are those for hexa
      allocate(ig_quad_gnode_xiloc(ilnnqu),ig_quad_gnode_etloc(ilnnqu))
      if (ilnnqu == 4) then
         ig_quad_gnode_xiloc(1) = 1
         ig_quad_gnode_etloc(1) = 1
         ig_quad_gnode_xiloc(2) = 2
         ig_quad_gnode_etloc(2) = 1
         ig_quad_gnode_xiloc(3) = 2
         ig_quad_gnode_etloc(3) = 2
         ig_quad_gnode_xiloc(4) = 1
         ig_quad_gnode_etloc(4) = 2
      elseif (ilnnqu == 9) then
         ig_quad_gnode_xiloc(1) = 1
         ig_quad_gnode_etloc(1) = 1
         ig_quad_gnode_xiloc(2) = 3
         ig_quad_gnode_etloc(2) = 1
         ig_quad_gnode_xiloc(3) = 3
         ig_quad_gnode_etloc(3) = 3
         ig_quad_gnode_xiloc(4) = 1
         ig_quad_gnode_etloc(4) = 3
         ig_quad_gnode_xiloc(5) = 2
         ig_quad_gnode_etloc(5) = 1
         ig_quad_gnode_xiloc(6) = 3
         ig_quad_gnode_etloc(6) = 2
         ig_quad_gnode_xiloc(7) = 2
         ig_quad_gnode_etloc(7) = 3
         ig_quad_gnode_xiloc(8) = 1
         ig_quad_gnode_etloc(8) = 2
         ig_quad_gnode_xiloc(9) = 2
         ig_quad_gnode_etloc(9) = 2
      else
         write(info,'(a)') "error in init_element: invalid number of geometrical nodes for quad"
         call error_stop(info)
      endif
      !
      !
      !***********************************************************************
      !indirection from geom nodes to GLL nodes (see docs/Convention.html)
      !***********************************************************************
      !                 xi       eta      zeta
      !                 m   -     l   -    k
      !
      !geom node
      !    1            1         1        1
      !    2          IG_NGLL     1        1
      !    3          IG_NGLL  IG_NGLL     1
      !    4            1      IG_NGLL     1
      !    5            1         1     IG_NGLL
      !    6          IG_NGLL     1     IG_NGLL
      !    7          IG_NGLL  IG_NGLL  IG_NGLL
      !    8            1      IG_NGLL  IG_NGLL
      !

      ig_hexa_node2gll(1,1) = 1
      ig_hexa_node2gll(2,1) = 1
      ig_hexa_node2gll(3,1) = 1

      ig_hexa_node2gll(1,2) = 1
      ig_hexa_node2gll(2,2) = 1
      ig_hexa_node2gll(3,2) = IG_NGLL

      ig_hexa_node2gll(1,3) = 1
      ig_hexa_node2gll(2,3) = IG_NGLL
      ig_hexa_node2gll(3,3) = IG_NGLL

      ig_hexa_node2gll(1,4) = 1
      ig_hexa_node2gll(2,4) = IG_NGLL
      ig_hexa_node2gll(3,4) = 1

      ig_hexa_node2gll(1,5) = IG_NGLL
      ig_hexa_node2gll(2,5) = 1
      ig_hexa_node2gll(3,5) = 1

      ig_hexa_node2gll(1,6) = IG_NGLL
      ig_hexa_node2gll(2,6) = 1
      ig_hexa_node2gll(3,6) = IG_NGLL

      ig_hexa_node2gll(1,7) = IG_NGLL
      ig_hexa_node2gll(2,7) = IG_NGLL
      ig_hexa_node2gll(3,7) = IG_NGLL

      ig_hexa_node2gll(1,8) = IG_NGLL
      ig_hexa_node2gll(2,8) = IG_NGLL
      ig_hexa_node2gll(3,8) = 1

      return
!***********************************************************************************************************************************************************************************
   end subroutine init_element
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine reads integer in binary files *.dat using a buffer of size BUFFER_READ_SIZE.
!>@param u  : FORTRAN unit file to be read
!>@param ir : record number to read in buffer array
!>@param i  : output integer
!***********************************************************************************************************************************************************************************
   subroutine read_bin_int(u,ir,i)
!***********************************************************************************************************************************************************************************

      implicit none

      integer, intent(in   )  :: u
      integer, intent(inout)  :: ir
      integer, intent(  out)  :: i

      integer                 :: ios

      if (ir == 1) then
         read(u,iostat=ios) buffer
         i = buffer(ir)
      elseif (ir < BUFFER_READ_SIZE) then
         i = buffer(ir)
      else
         i = buffer(ir)
         ir = 0
      endif

      return 

!***********************************************************************************************************************************************************************************
   end subroutine read_bin_int
!***********************************************************************************************************************************************************************************:w


end module mod_init_mesh
