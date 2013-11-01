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
!!This file contains a module to compute information about receivers.

!>@brief
!!This module contains subroutines to compute information about receivers and to write receivers' time history into binary files.
module mod_receiver

   use mpi

   implicit none

   private

   public :: init_hexa_receiver
   public :: init_quad_receiver
   public :: search_closest_hexa_gll
   public :: search_closest_quad_gll
   public :: compute_info_hexa_receiver
   public :: compute_info_quad_receiver
   public :: write_receiver_output

   contains

!
!
!>@brief
!!This subroutine reads all receivers in file *.vor; determines to which cpu they belong; computes their local coordinates @f$\xi,\eta,\zeta@f$ in hexahedron elements
!!and open receivers' output files in which displacements, velocities and accelerations are written.
!>@return mod_global_variables::tg_receiver_hexa
!***********************************************************************************************************************************************************************************
   subroutine init_hexa_receiver()
!***********************************************************************************************************************************************************************************
      
      use mpi
      
      use mod_global_variables, only :&
                                       cg_prefix&
                                      ,tg_receiver_hexa&
                                      ,ig_nreceiver_hexa&
                                      ,IG_LST_UNIT&
                                      ,ig_myrank&
                                      ,ig_ncpu&
                                      ,type_receiver_hexa&
                                      ,ig_receiver_saving_incr&
                                      ,LG_OUTPUT_DEBUG_FILE&
                                      ,cg_myrank&
                                      ,ig_hexa_receiver_unit&
                                      ,get_newunit&
                                      ,ig_ndt
      
      implicit none
      
      real   , dimension(ig_ncpu)                   :: rldum1
      real                                                   :: maxdmin
      real                                                   :: rldmin
                                                             
      integer, dimension(1)                                  :: min_loc
      integer, dimension(ig_ncpu)                   :: ilrcpu
      integer                                                :: ios
      integer                                                :: irec
      integer                                                :: i
      integer                                                :: j
      integer                                                :: ilnrec
      integer, dimension(MPI_STATUS_SIZE)                    :: statut
                                                             
      character(len=1)                                       :: ctmp
      character(len=6)                                       :: crec

      type(type_receiver_hexa), allocatable, dimension(:) :: tlrinf

!
!---->initialize total number of receivers in cpu myrank
      ig_nreceiver_hexa = 0
!
!
!---->if file *.vor exists, then we read receiver x,y,z location
      open(unit=100,file=trim(cg_prefix)//".vor",status='old',iostat=ios)

      if( (ios == 0) .and. (ig_receiver_saving_incr <= ig_ndt) ) then

         if (ig_myrank == 0) then
            write(IG_LST_UNIT,'(a)') " "
            write(IG_LST_UNIT,'(a)') "information about receivers in file "//trim(cg_prefix)//".vor"
         endif

         do while (.true.)
            read(unit=100,fmt=*,iostat=ios) ctmp
            if (ios /= 0) exit
            ig_nreceiver_hexa = ig_nreceiver_hexa + 1
         enddo
         rewind(100)

         allocate(tlrinf(ig_nreceiver_hexa))
         ilnrec = 0

         do irec = 1,ig_nreceiver_hexa

!
!---------->initialize tlrinf
            tlrinf(irec)%x          = 0.0
            tlrinf(irec)%y          = 0.0
            tlrinf(irec)%z          = 0.0
            tlrinf(irec)%xi         = 0.0
            tlrinf(irec)%eta        = 0.0
            tlrinf(irec)%zeta       = 0.0
            tlrinf(irec)%dmin       = 0.0
            tlrinf(irec)%lag(:,:,:) = 0.0   
            tlrinf(irec)%gll(:,:,:) = 0     
            tlrinf(irec)%cpu        = 0
            tlrinf(irec)%iel        = 0
            tlrinf(irec)%kgll       = 0
            tlrinf(irec)%lgll       = 0
            tlrinf(irec)%mgll       = 0
            tlrinf(irec)%rglo       = 0

            read(unit=100,fmt=*) tlrinf(irec)%x,tlrinf(irec)%y,tlrinf(irec)%z
!
!---------->find the closest element and gll point to the receiver irec
            call search_closest_hexa_gll(tlrinf(irec)%x   ,tlrinf(irec)%y   ,tlrinf(irec)%z   ,tlrinf(irec)%dmin&
                                                     ,tlrinf(irec)%kgll,tlrinf(irec)%lgll,tlrinf(irec)%mgll,tlrinf(irec)%iel )
!
!---------->check to which cpu the receiver belongs (i.e., cpu which has the smallest dmin)
            call mpi_allgather(tlrinf(irec)%dmin,1,mpi_real,rldum1,1,mpi_real,mpi_comm_world,ios)
            min_loc = minloc(rldum1(1:ig_ncpu))

            if (min_loc(1) == ig_myrank+1) then
               tlrinf(irec)%cpu  = ig_myrank + 1
               ilnrec            = ilnrec + 1
               if (ig_myrank /= 0) call mpi_send(tlrinf(irec)%iel,1,mpi_integer,0,100,mpi_comm_world,ios)
            else
               tlrinf(irec)%cpu = min_loc(1)
            endif
!
!---------->write some info about receivers in *.lst
            if (ig_myrank == 0) then
               if (min_loc(1) /= ig_myrank+1) call mpi_recv(tlrinf(irec)%iel,1,mpi_integer,min_loc(1)-1,100,mpi_comm_world,statut,ios)
               write(IG_LST_UNIT,'(3(a,i8))') "receiver ",irec," computed by cpu ",min_loc(1)-1," hexa ",tlrinf(irec)%iel
            endif

         enddo
         close(100)
!
!------->transfer array tlrinf which knows all the receivers to array tg_receiver_hexa that needs only to know the receiver of cpu myrank
         if (ilnrec > 0) then
            allocate(tg_receiver_hexa(ilnrec))
            j = 0
            do irec = 1,ig_nreceiver_hexa
               if (tlrinf(irec)%cpu == ig_myrank+1) then
                  j = j + 1
                  tg_receiver_hexa(j)      = tlrinf(irec)
                  tg_receiver_hexa(j)%rglo = irec
               endif
            enddo
            ig_nreceiver_hexa = ilnrec
         else
            ig_nreceiver_hexa = 0
         endif
         deallocate(tlrinf)
!
!------->cpu 0 gather the number of receiver (ig_nreceiver_hexa) of the other cpus
         call mpi_gather(ig_nreceiver_hexa,1,mpi_integer,ilrcpu,1,mpi_integer,0,mpi_comm_world,ios)
!
!------->compute information for receivers that belongs to cpu myrank and send 'dmin' to cpu 0
         do irec = 1,ig_nreceiver_hexa

            call compute_info_hexa_receiver(tg_receiver_hexa(irec))

            if (ig_myrank /= 0) then
               call mpi_send(tg_receiver_hexa(irec)%dmin,1,mpi_real,0,101,mpi_comm_world,ios)
            endif

         enddo

         if (LG_OUTPUT_DEBUG_FILE) then
            open(unit=100,file="debug.cpu."//trim(cg_myrank)//".volume.receivers.info")
            do irec = 1,ig_nreceiver_hexa
               write(unit=100,fmt='(/,a,i10    )') "Receiver     ",tg_receiver_hexa(irec)%rglo
               write(unit=100,fmt='( a,i10     )') "    in core  ",tg_receiver_hexa(irec)%cpu-1
               write(unit=100,fmt='( a,i10     )') "    in hexa  ",tg_receiver_hexa(irec)%iel
               write(unit=100,fmt='(        3a )') "       X       ","       Y       ","       Z       "
               write(unit=100,fmt='(3(e14.7,1x))') tg_receiver_hexa(irec)%x,tg_receiver_hexa(irec)%y,tg_receiver_hexa(irec)%z
               write(unit=100,fmt='(        3a )') "      XI       ","      ETA      ","     ZETA      "
               write(unit=100,fmt='(3(e14.7,1x))') tg_receiver_hexa(irec)%xi,tg_receiver_hexa(irec)%eta,tg_receiver_hexa(irec)%zeta
            enddo
            close(100)
         endif
!
!------->cpu 0 gathers dmin of all cpus and write the maximum dmin in *.lst
         if (ig_myrank == 0) then
            maxdmin = 0.0
            do i = 1,ig_ncpu
               if (i == 1) then
                  do irec = 1,ig_nreceiver_hexa
                     maxdmin = max(maxdmin,tg_receiver_hexa(irec)%dmin)
                  enddo
               else
                  do irec = 1,ilrcpu(i)
                     call mpi_recv(rldmin,1,mpi_real,i-1,101,mpi_comm_world,statut,ios)
                     maxdmin = max(maxdmin,rldmin)
                  enddo 
               endif
            enddo
            write(IG_LST_UNIT,'(a,e14.7)') "maximum localisation error of all receivers in volume = ",maxdmin
            call flush(IG_LST_UNIT)
         endif
!
!------->open receivers file once and for all
         if ( ig_nreceiver_hexa > 0 ) then

            allocate(ig_hexa_receiver_unit(ig_nreceiver_hexa))

            do irec = 1,ig_nreceiver_hexa

               write(crec,'(i6.6)') tg_receiver_hexa(irec)%rglo

               open(unit       = get_newunit(ig_hexa_receiver_unit(irec))&
                   ,file       = trim(cg_prefix)//".vor."//trim(adjustl(crec))//".gpl"&
                   ,status     = 'replace'&
                   ,action     = 'write'&
                   ,access     = 'sequential'&
                   ,form       = 'unformatted'&
                   ,buffered   = 'yes'&
                   ,recordtype = 'stream')

            enddo

         endif

      else
         if (ig_myrank == 0) then
            write(IG_LST_UNIT,'(/a)') "no volume receiver computed"
            call flush(IG_LST_UNIT)
         endif
      endif
      
      return
!***********************************************************************************************************************************************************************************
   end subroutine init_hexa_receiver
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine reads all receivers in file *.fsr; determines to which cpu they belong; computes their local coordinates @f$\xi,\eta@f$ in quadrangle elements
!!and open receivers' output files in which displacements, velocities and accelerations are written.
!>@return mod_global_variables::tg_receiver_quad
!***********************************************************************************************************************************************************************************
   subroutine init_quad_receiver()
!***********************************************************************************************************************************************************************************
      
      use mpi
      
      use mod_global_variables, only :&
                                       cg_prefix&
                                      ,tg_receiver_quad&
                                      ,ig_quadf_gll_glonum&
                                      ,ig_nreceiver_quad&
                                      ,IG_LST_UNIT&
                                      ,ig_myrank&
                                      ,ig_ncpu&
                                      ,type_receiver_quad&
                                      ,ig_receiver_saving_incr&
                                      ,LG_OUTPUT_DEBUG_FILE&
                                      ,cg_myrank&
                                      ,ig_quad_receiver_unit&
                                      ,get_newunit&
                                      ,ig_ndt
      
      implicit none
      
      real   , dimension(ig_ncpu)                   :: rldum1
      real                                                   :: maxdmin
      real                                                   :: rldmin
                                                             
      integer, dimension(1)                                  :: min_loc
      integer, dimension(ig_ncpu)                   :: ilrcpu
      integer                                                :: ios
      integer                                                :: irec
      integer                                                :: i
      integer                                                :: j
      integer                                                :: ilnrec
      integer, dimension(MPI_STATUS_SIZE)                    :: statut
                                                             
      character(len=1)                                       :: ctmp
      character(len=6)                                       :: crec

      type(type_receiver_quad), allocatable, dimension(:) :: tlrinf

!
!---->initialize total number of receivers in cpu myrank
      ig_nreceiver_quad = 0
!
!
!---->if file *receiver exists, then we read receiver x,y,z location
      open(unit=100,file=trim(cg_prefix)//".fsr",status='old',iostat=ios)

      if( (ios == 0) .and. (ig_receiver_saving_incr <= ig_ndt) ) then

         if (ig_myrank == 0) then
            write(IG_LST_UNIT,'(a)') " "
            write(IG_LST_UNIT,'(a)') "information about receivers in file "//trim(cg_prefix)//".fsr"
         endif

         do while (.true.)
            read(unit=100,fmt=*,iostat=ios) ctmp
            if (ios /= 0) exit
            ig_nreceiver_quad = ig_nreceiver_quad + 1
         enddo
         rewind(100)

         allocate(tlrinf(ig_nreceiver_quad))
         ilnrec = 0

         do irec = 1,ig_nreceiver_quad

!
!---------->initialize tlrinf
            tlrinf(irec)%x        = 0.0
            tlrinf(irec)%y        = 0.0
            tlrinf(irec)%z        = 0.0
            tlrinf(irec)%xi       = 0.0
            tlrinf(irec)%eta      = 0.0
            tlrinf(irec)%dmin     = 0.0
            tlrinf(irec)%lag(:,:) = 0.0   
            tlrinf(irec)%gll(:,:) = 0     
            tlrinf(irec)%cpu      = 0
            tlrinf(irec)%iel      = 0
            tlrinf(irec)%lgll     = 0
            tlrinf(irec)%mgll     = 0
            tlrinf(irec)%rglo     = 0

            read(unit=100,fmt=*) tlrinf(irec)%x,tlrinf(irec)%y
!
!---------->find the closest element and gll point to the receiver irec
            call search_closest_quad_gll(tlrinf(irec)%x,tlrinf(irec)%y,ig_quadf_gll_glonum,tlrinf(irec)%dmin,tlrinf(irec)%lgll,tlrinf(irec)%mgll,tlrinf(irec)%iel)
!
!---------->check to which cpu the receiver belongs (i.e., cpu which has the smallest dmin)
            call mpi_allgather(tlrinf(irec)%dmin,1,mpi_real,rldum1,1,mpi_real,mpi_comm_world,ios)
            min_loc = minloc(rldum1(1:ig_ncpu))

            if (min_loc(1) == ig_myrank+1) then
               tlrinf(irec)%cpu  = ig_myrank + 1
               ilnrec            = ilnrec + 1
               if (ig_myrank /= 0) call mpi_send(tlrinf(irec)%iel,1,mpi_integer,0,110,mpi_comm_world,ios)
            else
               tlrinf(irec)%cpu = min_loc(1)
            endif
!
!---------->write some info about receivers in *.lst
            if (ig_myrank == 0) then
               if (min_loc(1) /= ig_myrank+1) call mpi_recv(tlrinf(irec)%iel,1,mpi_integer,min_loc(1)-1,110,mpi_comm_world,statut,ios)
               write(IG_LST_UNIT,'(3(a,i8))') "receiver ",irec," computed by cpu ",min_loc(1)-1," quad ",tlrinf(irec)%iel
            endif

         enddo
         close(100)
!
!------->transfer array tlrinf which knows all the receivers to array tg_receiver_quad that needs only to know the receiver of cpu myrank
         if (ilnrec > 0) then
            allocate(tg_receiver_quad(ilnrec))
            j = 0
            do irec = 1,ig_nreceiver_quad
               if (tlrinf(irec)%cpu == ig_myrank+1) then
                  j = j + 1
                  tg_receiver_quad(j)      = tlrinf(irec)
                  tg_receiver_quad(j)%rglo = irec
               endif
            enddo
            ig_nreceiver_quad = ilnrec
         else
            ig_nreceiver_quad = 0
         endif
         deallocate(tlrinf)
!
!------->cpu 0 gather the number of receiver (ig_nreceiver_quad) of the other cpus
         call mpi_gather(ig_nreceiver_quad,1,mpi_integer,ilrcpu,1,mpi_integer,0,mpi_comm_world,ios)
!
!------->compute information for receivers that belongs to cpu myrank and send 'dmin' to cpu 0
         do irec = 1,ig_nreceiver_quad

            call compute_info_quad_receiver(tg_receiver_quad(irec))

            if (ig_myrank /= 0) then
               call mpi_send(tg_receiver_quad(irec)%dmin,1,mpi_real,0,111,mpi_comm_world,ios)
            endif

         enddo

         if (LG_OUTPUT_DEBUG_FILE) then
            open(unit=100,file="debug.cpu."//trim(cg_myrank)//".freesurface.receivers.info")
            do irec = 1,ig_nreceiver_quad
               write(unit=100,fmt='(/,a,i10    )') "Receiver     ",tg_receiver_quad(irec)%rglo
               write(unit=100,fmt='( a,i10     )') "    in core  ",tg_receiver_quad(irec)%cpu-1
               write(unit=100,fmt='( a,i10     )') "    in quad  ",tg_receiver_quad(irec)%iel
               write(unit=100,fmt='(        2a )') "       X       ","       Y       "
               write(unit=100,fmt='(2(e14.7,1x))') tg_receiver_quad(irec)%x,tg_receiver_quad(irec)%y
               write(unit=100,fmt='(        2a )') "      XI       ","      ETA      "
               write(unit=100,fmt='(2(e14.7,1x))') tg_receiver_quad(irec)%xi,tg_receiver_quad(irec)%eta
            enddo
            close(100)
         endif
!
!------->cpu 0 gathers dmin of all cpus and write the maximum dmin in *.lst
         if (ig_myrank == 0) then
            maxdmin = 0.0
            do i = 1,ig_ncpu
               if (i == 1) then
                  do irec = 1,ig_nreceiver_quad
                     maxdmin = max(maxdmin,tg_receiver_quad(irec)%dmin)
                  enddo
               else
                  do irec = 1,ilrcpu(i)
                     call mpi_recv(rldmin,1,mpi_real,i-1,111,mpi_comm_world,statut,ios)
                     maxdmin = max(maxdmin,rldmin)
                  enddo 
               endif
            enddo
            write(IG_LST_UNIT,'(a,e14.7)') "maximum localisation error of all receivers in free surface = ",maxdmin
            call flush(IG_LST_UNIT)
         endif
!
!------->open receivers file once and for all
         if ( ig_nreceiver_quad > 0 ) then

            allocate(ig_quad_receiver_unit(ig_nreceiver_quad))

            do irec = 1,ig_nreceiver_quad

               write(crec,'(i6.6)') tg_receiver_quad(irec)%rglo

               open(unit       = get_newunit(ig_quad_receiver_unit(irec))&
                   ,file       = trim(cg_prefix)//".fsr."//trim(adjustl(crec))//".gpl"&
                   ,status     = 'replace'&
                   ,action     = 'write'&
                   ,access     = 'sequential'&
                   ,form       = 'unformatted'&
                   ,buffered   = 'yes'&
                   ,recordtype = 'stream')

            enddo

         endif

      else
         if (ig_myrank == 0) then
            write(IG_LST_UNIT,'(/a)') "no free surface receiver computed"
            call flush(IG_LST_UNIT)
         endif
      endif
      
      return
!***********************************************************************************************************************************************************************************
   end subroutine init_quad_receiver
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine searches among all hexahedron elements in cpu myrank the closest GLL node to a point of coordinates @f$x,y,z@f$
!>@param x    : @f$x@f$-coordinates of a point
!>@param y    : @f$y@f$-coordinates of a point
!>@param z    : @f$z@f$-coordinates of a point
!>@param dmin : minimal Euclidean distance from point to all GLL nodes in cpu myrank
!>@param kgll : closest GLL node index along @f$\zeta@f$-coordinate 
!>@param lgll : closest GLL node index along @f$\eta @f$-coordinate
!>@param mgll : closest GLL node index along @f$\xi  @f$-coordinate
!>@param iel  : number of hexahedron element in cpu myrank whose GLL node has the minimal dmin
!***********************************************************************************************************************************************************************************
   subroutine search_closest_hexa_gll(x,y,z,dmin,kgll,lgll,mgll,iel)
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only :&
                                      IG_NGLL&
                                     ,rg_gll_coordinate&
                                     ,ig_nhexa&
                                     ,ig_hexa_gll_glonum
      
      implicit none
      
      real   , intent(in)  :: x
      real   , intent(in)  :: y
      real   , intent(in)  :: z
      real   , intent(out) :: dmin
      integer, intent(out) :: kgll
      integer, intent(out) :: lgll
      integer, intent(out) :: mgll
      integer, intent(out) :: iel
      
      real                 :: d
      integer              :: ihexa
      integer              :: k
      integer              :: l
      integer              :: m
      integer              :: igll

      dmin = huge(dmin)

      do ihexa = 1,ig_nhexa
         do k = 1,IG_NGLL
            do l = 1,IG_NGLL
               do m = 1,IG_NGLL

                  igll = ig_hexa_gll_glonum(m,l,k,ihexa)

                  d = sqrt( (rg_gll_coordinate(1,igll) - x)**2 &
                           +(rg_gll_coordinate(2,igll) - y)**2 &
                           +(rg_gll_coordinate(3,igll) - z)**2 )
                  if (d < dmin) then
                     dmin = d
                     kgll = k
                     lgll = l
                     mgll = m
                     iel  = ihexa
                  endif

               enddo
            enddo
         enddo
      enddo
      
      return
!***********************************************************************************************************************************************************************************
   end subroutine search_closest_hexa_gll
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine searches among all free surface quadrangle elements in cpu myrank the closest GLL node to a point of coordinates @f$x,y@f$
!>@param x          : @f$x@f$-coordinates of a point
!>@param y          : @f$y@f$-coordinates of a point
!>@param global_gll : indirection from local GLL nodes to global GLL nodes for quadrangle elements
!>@param dmin       : minimal Euclidean distance from point to all GLL nodes in cpu myrank
!>@param lgll       : closest GLL node index along @f$\eta @f$-coordinate
!>@param mgll       : closest GLL node index along @f$\xi  @f$-coordinate
!>@param iel        : number of free surface quadrangle element in cpu myrank whose GLL node has the minimal dmin
!***********************************************************************************************************************************************************************************
   subroutine search_closest_quad_gll(x,y,global_gll,dmin,lgll,mgll,iel)
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only :&
                                      IG_NGLL&
                                     ,rg_gll_coordinate
      
      implicit none
      
      real   , intent(in)                   :: x
      real   , intent(in)                   :: y
      integer, intent(in), dimension(:,:,:) :: global_gll

      real   , intent(out)                  :: dmin
      integer, intent(out)                  :: lgll
      integer, intent(out)                  :: mgll
      integer, intent(out)                  :: iel
                                            
      real                                  :: d
                                            
      integer                               :: iquad
      integer                               :: l
      integer                               :: m
      integer                               :: nquad
      integer                               :: igll

      dmin  = huge(dmin)
      nquad = size(global_gll,3)

      do iquad = 1,nquad
         do l = 1,IG_NGLL
            do m = 1,IG_NGLL

               igll  = global_gll(m,l,iquad)
               d     = sqrt( (rg_gll_coordinate(1,igll) - x)**2 + (rg_gll_coordinate(2,igll) - y)**2 )

               if (d < dmin) then
                  dmin = d
                  lgll = l
                  mgll = m
                  iel  = iquad
               endif

            enddo
         enddo
      enddo

      return
!***********************************************************************************************************************************************************************************
   end subroutine search_closest_quad_gll
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine computes local coordinates @f$\xi,\eta,\zeta@f$ of a receiver in the hexahedron element it belongs.
!!It also stores global GLL nodes numbers for future interpolations as well as Lagrange polynonial values at the receiver.
!>@param tlxinf : data structure for a receiver located inside a hexahedron element
!***********************************************************************************************************************************************************************************
   subroutine compute_info_hexa_receiver(tlxinf)
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only :&
                                      rg_gll_abscissa&
                                     ,IG_NGLL&
                                     ,ig_hexa_gll_glonum&
                                     ,type_receiver_hexa

      use mod_jacobian        , only : compute_hexa_jacobian 

      use mod_lagrange        , only : lagrap
      
      use mod_coordinate      , only : compute_hexa_point_coord
      
      implicit none
      
      real                                       :: eps
      real                                       :: dmin
      real                                       :: xisol
      real                                       :: etsol
      real                                       :: zesol
      real                                       :: newx
      real                                       :: newy
      real                                       :: newz
      real                                       :: dxidx
      real                                       :: dxidy
      real                                       :: dxidz
      real                                       :: detdx
      real                                       :: detdy
      real                                       :: detdz
      real                                       :: dzedx
      real                                       :: dzedy
      real                                       :: dzedz
      real                                       :: dx
      real                                       :: dy
      real                                       :: dz
      real                                       :: dxi
      real                                       :: det
      real                                       :: dze
      real                                       :: deter
                                                 
      integer                                    :: ihexa 
      integer                                    :: igll
      integer                                    :: k
      integer                                    :: l
      integer                                    :: m
      integer                                    :: iter
      integer, parameter                         :: ITER_MAX=100
    
      type(type_receiver_hexa), intent(inout) :: tlxinf

!
!---->initialize element number and coordinates (needed if dmin < eps)
      ihexa= tlxinf%iel
      newx = tlxinf%x
      newy = tlxinf%y
      newz = tlxinf%z
      eps  = 1.0e-2
      iter = 0
!
!---->solve the nonlinear system to find local coordinates
      xisol = rg_gll_abscissa(tlxinf%mgll)
      etsol = rg_gll_abscissa(tlxinf%lgll)
      zesol = rg_gll_abscissa(tlxinf%kgll)
      dmin  = tlxinf%dmin

      do while(dmin > eps)

         call compute_hexa_jacobian(ihexa,xisol,etsol,zesol,dxidx,dxidy,dxidz,detdx,detdy,detdz,dzedx,dzedy,dzedz,deter)

         call compute_hexa_point_coord(ihexa,xisol,etsol,zesol,newx,newy,newz)

         dx    = tlxinf%x - newx
         dy    = tlxinf%y - newy
         dz    = tlxinf%z - newz
         dxi   = dxidx*dx + dxidy*dy + dxidz*dz
         det   = detdx*dx + detdy*dy + detdz*dz
         dze   = dzedx*dx + dzedy*dy + dzedz*dz
         xisol = xisol + dxi
         etsol = etsol + det
         zesol = zesol + dze
         dmin  = sqrt( (newx-tlxinf%x)**2 + (newy-tlxinf%y)**2 + (newz-tlxinf%z)**2 )

         iter  = iter + 1

         if (mod(iter,ITER_MAX) == 0) then
            eps = eps*2.0
         endif

      enddo

      tlxinf%xi   = xisol
      tlxinf%eta  = etsol
      tlxinf%zeta = zesol
      tlxinf%x    = newx
      tlxinf%y    = newy
      tlxinf%z    = newz
      tlxinf%dmin = dmin
!
!---->store global GLL nodes number
      do k = 1,IG_NGLL
         do l = 1,IG_NGLL
            do m = 1,IG_NGLL
               igll              = ig_hexa_gll_glonum(m,l,k,ihexa)
               tlxinf%gll(m,l,k) = igll
            enddo
         enddo
      enddo
!     
!---->compute and store Lagrange polynomial values of a receiver at xisol, etsol and zesol
      do k = 1,IG_NGLL
         do l = 1,IG_NGLL
            do m = 1,IG_NGLL
               tlxinf%lag(m,l,k) = lagrap(m,tlxinf%xi,IG_NGLL)*lagrap(l,tlxinf%eta,IG_NGLL)*lagrap(k,tlxinf%zeta,IG_NGLL)
            enddo
         enddo
      enddo

      return
!***********************************************************************************************************************************************************************************
   end subroutine compute_info_hexa_receiver
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine computes local coordinates @f$\xi,\eta@f$ of a receiver in the quadrangle element it belongs.
!!It also stores global GLL nodes numbers for future interpolations as well as Lagrange polynonial values at the receiver.
!>@param tlxinf : data structure for a receiver located inside a quadrangle element
!***********************************************************************************************************************************************************************************
   subroutine compute_info_quad_receiver(tlxinf)
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only : rg_gll_abscissa&
                                      ,IG_NGLL&
                                      ,ig_quadf_gll_glonum&
                                      ,type_receiver_quad

      use mod_jacobian        , only : compute_quad_jacobian

      use mod_lagrange        , only : lagrap
      
      use mod_coordinate      , only :&
                                      compute_quad_point_coord&
                                     ,compute_quad_point_coord_z
      
      implicit none
      
      real                                       :: eps
      real                                       :: dmin
      real                                       :: xisol
      real                                       :: etsol
      real                                       :: newx
      real                                       :: newy
      real                                       :: newz
      real                                       :: dxidx
      real                                       :: dxidy
      real                                       :: detdx
      real                                       :: detdy
      real                                       :: deter
      real                                       :: dx
      real                                       :: dy
      real                                       :: dxi
      real                                       :: det
                                                 
      integer                                    :: k
      integer                                    :: l
      integer                                    :: iquad
      integer                                    :: igll
      integer                                    :: iter
      integer, parameter                         :: ITER_MAX=100
    
      type(type_receiver_quad), intent(inout) :: tlxinf

!
!---->initialize coordinates, epsilon and iteration
      newx  = tlxinf%x
      newy  = tlxinf%y
      dmin  = tlxinf%dmin
      iquad = tlxinf%iel
      xisol = rg_gll_abscissa(tlxinf%mgll)
      etsol = rg_gll_abscissa(tlxinf%lgll)
      eps   = 1.0e-2
      iter  = 0
!
!---->solve the nonlinear system to find local coordinates
      do while(dmin > eps)

         call compute_quad_jacobian(iquad,xisol,etsol,dxidx,dxidy,detdx,detdy,deter)

         call compute_quad_point_coord(iquad,xisol,etsol,newx,newy,newz)

         dx    = tlxinf%x - newx
         dy    = tlxinf%y - newy
         dxi   = dxidx*dx + dxidy*dy
         det   = detdx*dx + detdy*dy
         xisol = xisol + dxi
         etsol = etsol + det
         dmin  = sqrt( (newx-tlxinf%x)**2 + (newy-tlxinf%y)**2 )

         iter  = iter + 1
         if (mod(iter,ITER_MAX) == 0) then
            eps = eps*2.0
         endif

      enddo

      if (xisol < -1.0) xisol = -1.0
      if (xisol > +1.0) xisol = +1.0
      tlxinf%xi   = xisol

      if (etsol < -1.0) etsol = -1.0
      if (etsol > +1.0) etsol = +1.0
      tlxinf%eta  = etsol

      tlxinf%x    = newx
      tlxinf%y    = newy
      tlxinf%z    = compute_quad_point_coord_z(iquad,xisol,etsol)
      tlxinf%dmin = dmin
!
!---->store global GLL nodes number
      do k = 1,IG_NGLL
         do l = 1,IG_NGLL
            igll            = ig_quadf_gll_glonum(l,k,iquad)
            tlxinf%gll(l,k) = igll
         enddo
      enddo
!     
!---->compute and store Lagrange polynomial values of a receiver at xisol and etsol
      do k = 1,IG_NGLL
         do l = 1,IG_NGLL
            tlxinf%lag(l,k) = lagrap(l,tlxinf%xi,IG_NGLL)*lagrap(k,tlxinf%eta,IG_NGLL)
         enddo
      enddo

      return
!***********************************************************************************************************************************************************************************
   end subroutine compute_info_quad_receiver
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine writes @f$x,y,z@f$-displacements, velocities and accelerations at receivers
!***********************************************************************************************************************************************************************************
   subroutine write_receiver_output()
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only :&
                                      ig_nreceiver_hexa&
                                     ,ig_nreceiver_quad&
                                     ,rg_simu_current_time&
                                     ,tg_receiver_hexa&
                                     ,tg_receiver_quad&
                                     ,rg_gll_displacement&
                                     ,rg_gll_velocity&
                                     ,rg_gll_acceleration&
                                     ,IG_NGLL&
                                     ,IG_NDOF&
                                     ,ig_hexa_receiver_unit&
                                     ,ig_quad_receiver_unit

      use mod_gll_value

      use mod_lagrange, only :&
                              hexa_lagrange_interp&
                             ,quad_lagrange_interp

      implicit none 

      real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL) :: gll_dis_hexa
      real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL) :: gll_vel_hexa
      real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL) :: gll_acc_hexa

      real, dimension(IG_NDOF,IG_NGLL,IG_NGLL)         :: gll_dis_quad
      real, dimension(IG_NDOF,IG_NGLL,IG_NGLL)         :: gll_vel_quad
      real, dimension(IG_NDOF,IG_NGLL,IG_NGLL)         :: gll_acc_quad

      real                                             :: dis_x
      real                                             :: dis_y
      real                                             :: dis_z
      real                                             :: vel_x
      real                                             :: vel_y
      real                                             :: vel_z
      real                                             :: acc_x
      real                                             :: acc_y
      real                                             :: acc_z

      integer                                          :: irec

!
!
!********************************************************
!---->receivers inside hexahedron elements
!********************************************************
      do irec = 1,ig_nreceiver_hexa

         call get_hexa_gll_value(rg_gll_displacement,tg_receiver_hexa(irec)%gll,gll_dis_hexa)
         call get_hexa_gll_value(rg_gll_velocity    ,tg_receiver_hexa(irec)%gll,gll_vel_hexa)
         call get_hexa_gll_value(rg_gll_acceleration,tg_receiver_hexa(irec)%gll,gll_acc_hexa)

!
!------->interpolate displacement, velocity and acceleration for receivers inside cpu myrank
         call hexa_lagrange_interp(gll_dis_hexa,tg_receiver_hexa(irec)%lag,dis_x,dis_y,dis_z)
         call hexa_lagrange_interp(gll_vel_hexa,tg_receiver_hexa(irec)%lag,vel_x,vel_y,vel_z)
         call hexa_lagrange_interp(gll_acc_hexa,tg_receiver_hexa(irec)%lag,acc_x,acc_y,acc_z)
!
!------->write binary values (files have been open inside the subroutine init_hexa_receiver)
         write(unit=ig_hexa_receiver_unit(irec)) rg_simu_current_time,dis_x,dis_y,dis_z,vel_x,vel_y,vel_z,acc_x,acc_y,acc_z
         call flush(ig_hexa_receiver_unit(irec))

      enddo

!
!
!********************************************************
!---->receivers inside quadrangle elements
!********************************************************
      do irec = 1,ig_nreceiver_quad

         call get_quad_gll_value(rg_gll_displacement,tg_receiver_quad(irec)%gll,gll_dis_quad)
         call get_quad_gll_value(rg_gll_velocity    ,tg_receiver_quad(irec)%gll,gll_vel_quad)
         call get_quad_gll_value(rg_gll_acceleration,tg_receiver_quad(irec)%gll,gll_acc_quad)

!
!------->interpolate displacement, velocity and acceleration for receivers inside cpu myrank
         call quad_lagrange_interp(gll_dis_quad,tg_receiver_quad(irec)%lag,dis_x,dis_y,dis_z)
         call quad_lagrange_interp(gll_vel_quad,tg_receiver_quad(irec)%lag,vel_x,vel_y,vel_z)
         call quad_lagrange_interp(gll_acc_quad,tg_receiver_quad(irec)%lag,acc_x,acc_y,acc_z)
!
!------->write binary values (files have been open inside the subroutine init_quad_receiver)
         write(unit=ig_quad_receiver_unit(irec)) rg_simu_current_time,dis_x,dis_y,dis_z,vel_x,vel_y,vel_z,acc_x,acc_y,acc_z
         call flush(ig_quad_receiver_unit(irec))

      enddo

      return

!***********************************************************************************************************************************************************************************
   end subroutine write_receiver_output
!***********************************************************************************************************************************************************************************

end module mod_receiver
