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
!!This file contains a module to compute information about sources.

!>@brief
!!This module contains subroutines to compute information about sources.
module mod_source

   use mpi

   implicit none

   public  :: compute_double_couple_source
   public  :: init_double_couple_source
   public  :: init_single_force_source
   private :: compute_source_local_coordinate

   contains

!
!
!>@brief
!!This subroutine pre-computes all double couple point sources defined by type mod_global_variables::type_double_couple_source as equivalent single force point sources (@f$ F^{ext}_{n+1} @f$)
!!distributed over hexahedron elements (see mod_global_variables::rg_dcsource_gll_force).
!>@return mod_global_variables::rg_dcsource_gll_force
!>@reference : Komatitsch and Tromp, 1999 \cite Komatitsch1999a; Komatitsch and Tromp, 2002 \cite Komatitsch2002b.
!***********************************************************************************************************************************************************************************
   subroutine compute_double_couple_source()
!***********************************************************************************************************************************************************************************
      
      use mpi
   
      use mod_global_variables, only : ig_ndcsource&
                                      ,IG_NGLL&
                                      ,IG_NDOF&
                                      ,tg_dcsource&
                                      ,rg_dcsource_gll_force&
                                      ,rg_hexa_gll_dxidx&
                                      ,rg_hexa_gll_dxidy&
                                      ,rg_hexa_gll_dxidz&
                                      ,rg_hexa_gll_detdx&
                                      ,rg_hexa_gll_detdy&
                                      ,rg_hexa_gll_detdz&
                                      ,rg_hexa_gll_dzedx&
                                      ,rg_hexa_gll_dzedy&
                                      ,rg_hexa_gll_dzedz&
                                      ,ig_myrank&
                                      ,error_stop&
                                      ,IG_LST_UNIT

      use mod_init_memory

      use mod_lagrange, only : lagrad,lagrap

      implicit none

      real    :: ftmp(3,3,IG_NGLL,IG_NGLL,IG_NGLL)
              
      integer :: iso
      integer :: iel
      integer :: k
      integer :: l
      integer :: m
      integer :: n
      integer :: o
      integer :: p
      integer :: ios



      !**************************************************************************************************************!
      !References: Komatitsch and Tromp, 1999 \cite Komatitsch1999a; Komatitsch and Tromp, 2002 \cite Komatitsch2002b!
      !**************************************************************************************************************!


      !
      !
      !***********************************************************
      !implementation of the double couple point source 
      !***********************************************************
      if (ig_myrank == 0) then
         write(IG_LST_UNIT,'(" ",/,a)') "computing double-couple sources if any..."
         call flush(IG_LST_UNIT)
      endif
   
      if (ig_ndcsource > 0) then

         ios = init_array_real(rg_dcsource_gll_force,ig_ndcsource,IG_NGLL,IG_NGLL,IG_NGLL,IG_NDOF,"rg_dcsource_gll_force")

      endif

      !
      !inject double couple points sources at any location (i.e., not necessarily at a gll node)
      do iso = 1,ig_ndcsource
   
         iel = tg_dcsource(iso)%iel
      ! 
      !->compute first part at gll node first, then an interpolation (cf. below) is done at the xi, eta, zeta of the source
         do k = 1,IG_NGLL
            do l = 1,IG_NGLL
               do m = 1,IG_NGLL
                  ftmp(1,1,m,l,k) = tg_dcsource(iso)%mxx*rg_hexa_gll_dxidx(m,l,k,iel) &
                                  + tg_dcsource(iso)%mxy*rg_hexa_gll_dxidy(m,l,k,iel) &
                                  + tg_dcsource(iso)%mxz*rg_hexa_gll_dxidz(m,l,k,iel)
                  ftmp(2,1,m,l,k) = tg_dcsource(iso)%mxx*rg_hexa_gll_detdx(m,l,k,iel) &
                                  + tg_dcsource(iso)%mxy*rg_hexa_gll_detdy(m,l,k,iel) &
                                  + tg_dcsource(iso)%mxz*rg_hexa_gll_detdz(m,l,k,iel)
                  ftmp(3,1,m,l,k) = tg_dcsource(iso)%mxx*rg_hexa_gll_dzedx(m,l,k,iel) &
                                  + tg_dcsource(iso)%mxy*rg_hexa_gll_dzedy(m,l,k,iel) &
                                  + tg_dcsource(iso)%mxz*rg_hexa_gll_dzedz(m,l,k,iel)
       
                  ftmp(1,2,m,l,k) = tg_dcsource(iso)%mxy*rg_hexa_gll_dxidx(m,l,k,iel) &
                                  + tg_dcsource(iso)%myy*rg_hexa_gll_dxidy(m,l,k,iel) &
                                  + tg_dcsource(iso)%myz*rg_hexa_gll_dxidz(m,l,k,iel)
                  ftmp(2,2,m,l,k) = tg_dcsource(iso)%mxy*rg_hexa_gll_detdx(m,l,k,iel) &
                                  + tg_dcsource(iso)%myy*rg_hexa_gll_detdy(m,l,k,iel) &
                                  + tg_dcsource(iso)%myz*rg_hexa_gll_detdz(m,l,k,iel)
                  ftmp(3,2,m,l,k) = tg_dcsource(iso)%mxy*rg_hexa_gll_dzedx(m,l,k,iel) &
                                  + tg_dcsource(iso)%myy*rg_hexa_gll_dzedy(m,l,k,iel) &
                                  + tg_dcsource(iso)%myz*rg_hexa_gll_dzedz(m,l,k,iel)
       
                  ftmp(1,3,m,l,k) = tg_dcsource(iso)%mxz*rg_hexa_gll_dxidx(m,l,k,iel) &
                                  + tg_dcsource(iso)%myz*rg_hexa_gll_dxidy(m,l,k,iel) &
                                  + tg_dcsource(iso)%mzz*rg_hexa_gll_dxidz(m,l,k,iel)
                  ftmp(2,3,m,l,k) = tg_dcsource(iso)%mxz*rg_hexa_gll_detdx(m,l,k,iel) &
                                  + tg_dcsource(iso)%myz*rg_hexa_gll_detdy(m,l,k,iel) &
                                  + tg_dcsource(iso)%mzz*rg_hexa_gll_detdz(m,l,k,iel)
                  ftmp(3,3,m,l,k) = tg_dcsource(iso)%mxz*rg_hexa_gll_dzedx(m,l,k,iel) &
                                  + tg_dcsource(iso)%myz*rg_hexa_gll_dzedy(m,l,k,iel) &
                                  + tg_dcsource(iso)%mzz*rg_hexa_gll_dzedz(m,l,k,iel)
               enddo
            enddo
         enddo
      ! 
      !->interpolation at the xi, eta, zeta of the source + multiplication by the test function
         do k = 1,IG_NGLL
            do l = 1,IG_NGLL
               do m = 1,IG_NGLL
                  rg_dcsource_gll_force(:,m,l,k,iso) = 0.0
                  do n = 1,IG_NGLL
                     do o = 1,IG_NGLL
                        do p = 1,IG_NGLL
                           rg_dcsource_gll_force(1,m,l,k,iso) = rg_dcsource_gll_force(1,m,l,k,iso) &
                           + lagrap(p,tg_dcsource(iso)%xi,IG_NGLL) &
                           * lagrap(o,tg_dcsource(iso)%et,IG_NGLL) &
                           * lagrap(n,tg_dcsource(iso)%ze,IG_NGLL) &
                           *(ftmp(1,1,p,o,n) * lagrad(m,tg_dcsource(iso)%xi,IG_NGLL) &
                                             * lagrap(l,tg_dcsource(iso)%et,IG_NGLL) &
                                             * lagrap(k,tg_dcsource(iso)%ze,IG_NGLL) &
                           + ftmp(2,1,p,o,n) * lagrap(m,tg_dcsource(iso)%xi,IG_NGLL) &
                                             * lagrad(l,tg_dcsource(iso)%et,IG_NGLL) &
                                             * lagrap(k,tg_dcsource(iso)%ze,IG_NGLL) &
                           + ftmp(3,1,p,o,n) * lagrap(m,tg_dcsource(iso)%xi,IG_NGLL) &
                                             * lagrap(l,tg_dcsource(iso)%et,IG_NGLL) &
                                             * lagrad(k,tg_dcsource(iso)%ze,IG_NGLL))
        
                           rg_dcsource_gll_force(2,m,l,k,iso) = rg_dcsource_gll_force(2,m,l,k,iso) &
                           + lagrap(p,tg_dcsource(iso)%xi,IG_NGLL) &
                           * lagrap(o,tg_dcsource(iso)%et,IG_NGLL) &
                           * lagrap(n,tg_dcsource(iso)%ze,IG_NGLL) &
                           *(ftmp(1,2,p,o,n) * lagrad(m,tg_dcsource(iso)%xi,IG_NGLL) &
                                             * lagrap(l,tg_dcsource(iso)%et,IG_NGLL) &
                                             * lagrap(k,tg_dcsource(iso)%ze,IG_NGLL) &
                           + ftmp(2,2,p,o,n) * lagrap(m,tg_dcsource(iso)%xi,IG_NGLL) &
                                             * lagrad(l,tg_dcsource(iso)%et,IG_NGLL) &
                                             * lagrap(k,tg_dcsource(iso)%ze,IG_NGLL) &
                           + ftmp(3,2,p,o,n) * lagrap(m,tg_dcsource(iso)%xi,IG_NGLL) &
                                             * lagrap(l,tg_dcsource(iso)%et,IG_NGLL) &
                                             * lagrad(k,tg_dcsource(iso)%ze,IG_NGLL)) 
       
                           rg_dcsource_gll_force(3,m,l,k,iso) = rg_dcsource_gll_force(3,m,l,k,iso) &
                           + lagrap(p,tg_dcsource(iso)%xi,IG_NGLL) &
                           * lagrap(o,tg_dcsource(iso)%et,IG_NGLL) &
                           * lagrap(n,tg_dcsource(iso)%ze,IG_NGLL) &
                           *(ftmp(1,3,p,o,n) * lagrad(m,tg_dcsource(iso)%xi,IG_NGLL) &
                                             * lagrap(l,tg_dcsource(iso)%et,IG_NGLL) &
                                             * lagrap(k,tg_dcsource(iso)%ze,IG_NGLL) &
                           + ftmp(2,3,p,o,n) * lagrap(m,tg_dcsource(iso)%xi,IG_NGLL) &
                                             * lagrad(l,tg_dcsource(iso)%et,IG_NGLL) &
                                             * lagrap(k,tg_dcsource(iso)%ze,IG_NGLL) &
                           + ftmp(3,3,p,o,n) * lagrap(m,tg_dcsource(iso)%xi,IG_NGLL) &
                                             * lagrap(l,tg_dcsource(iso)%et,IG_NGLL) &
                                             * lagrad(k,tg_dcsource(iso)%ze,IG_NGLL)) 
                        enddo
                     enddo
                  enddo
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
   end subroutine compute_double_couple_source
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine reads all double couple point sources in file *.dcs; sets double couple point sources moment tensor; determines to which cpu belong double couple point sources and
!!computes local coordinates @f$\xi,\eta,\zeta@f$ of double couple point sources inside the hexahedron element it belongs.
!>@return mod_global_variables::tg_dcsource
!>@return mod_global_variables::ig_ndcsource
!>@return mod_global_variables::rg_dcsource_user_func
!***********************************************************************************************************************************************************************************
   subroutine init_double_couple_source()
!***********************************************************************************************************************************************************************************
      
      use mpi
      
      use mod_global_variables, only :&
                                      IG_LST_UNIT&
                                     ,tg_dcsource&
                                     ,ig_ndcsource&
                                     ,ig_ncpu&
                                     ,ig_myrank&
                                     ,type_double_couple_source&
                                     ,cg_prefix&
                                     ,rg_dt&
                                     ,ig_ndt&
                                     ,RG_PI&
                                     ,rg_dcsource_user_func
      
      use mod_source_function , only : interp_linear
   
      use mod_receiver        , only : search_closest_hexa_gll 
      
      implicit none
      
      real, allocatable, dimension(:)               :: dum1
      real, allocatable, dimension(:)               :: dum2
      real, allocatable, dimension(:)               :: dum3
      real             , dimension(ig_ncpu)         :: dummy
      real                                          :: st
      real                                          :: di
      real                                          :: ra
      real                                          :: m0
      real                                          :: rldmin
      real                                          :: maxdmin
                                            
      integer          , dimension(MPI_STATUS_SIZE) :: statut
      integer                                       :: min_loc(1)
      integer          , dimension(ig_ncpu)         :: ilrcpu
      integer                                       :: i
      integer                                       :: j
      integer                                       :: iso
      integer                                       :: ios
      integer                                       :: ios1
      integer                                       :: ilndco
      integer                                       :: itempo
                                                    
      character(len=1)                              :: ctempo
                                                    
      type(type_double_couple_source), allocatable  :: tldoco(:)
   
      !
      !
      !********************************************************************
      !find element containing double couple point source in cpu 'myrank'
      !********************************************************************

      ilndco = 0

      open(unit=100,file=trim(cg_prefix)//".dcs",status='old',iostat=ios)
   
      if (ios == 0) then

         if (ig_myrank == 0) write(IG_LST_UNIT,'(a)') " "
   
         read(unit=100,fmt='(i10)') ig_ndcsource

         allocate(tldoco(ig_ndcsource))
   
         do iso = 1,ig_ndcsource
   
            !initialize tldoco
            tldoco(iso)%x           = 0.0
            tldoco(iso)%y           = 0.0
            tldoco(iso)%z           = 0.0
            tldoco(iso)%mxx         = 0.0       
            tldoco(iso)%myy         = 0.0       
            tldoco(iso)%mzz         = 0.0
            tldoco(iso)%mxy         = 0.0
            tldoco(iso)%mxz         = 0.0
            tldoco(iso)%myz         = 0.0
            tldoco(iso)%shift_time  = 0.0
            tldoco(iso)%rise_time   = 0.0
            tldoco(iso)%xi          = 0.0
            tldoco(iso)%et          = 0.0
            tldoco(iso)%ze          = 0.0
            tldoco(iso)%dmin        = 0.0
            tldoco(iso)%str         = 0.0
            tldoco(iso)%dip         = 0.0
            tldoco(iso)%rak         = 0.0
            tldoco(iso)%mw          = 0.0
            tldoco(iso)%icur        = 0
            tldoco(iso)%cpu         = 0
            tldoco(iso)%iel         = 0
            tldoco(iso)%kgll        = 0
            tldoco(iso)%lgll        = 0
            tldoco(iso)%mgll        = 0
   
   
            read(unit=100,fmt=*) tldoco(iso)%x,tldoco(iso)%y,tldoco(iso)%z,tldoco(iso)%mw,tldoco(iso)%str,tldoco(iso)%dip,tldoco(iso)%rak,tldoco(iso)%shift_time,tldoco(iso)%rise_time,tldoco(iso)%icur
   
            st              = tldoco(iso)%str*RG_PI/180.0
            di              = tldoco(iso)%dip*RG_PI/180.0
            ra              = tldoco(iso)%rak*RG_PI/180.0
            m0              = 10**(1.5*tldoco(iso)%mw + 9.1)
            tldoco(iso)%mxx = +m0*(sin(di)*cos(ra)*sin(2.0*st) -     sin(2.0*di)*sin(ra)*(cos(st))**2) !+myy for aki
            tldoco(iso)%mxy = +m0*(sin(di)*cos(ra)*cos(2.0*st) + 0.5*sin(2.0*di)*sin(ra)* sin(2.0*st)) !+mxy 
            tldoco(iso)%mxz = +m0*(cos(di)*cos(ra)*sin(    st) -     cos(2.0*di)*sin(ra)* cos(st))     !-myz 
            tldoco(iso)%myy = -m0*(sin(di)*cos(ra)*sin(2.0*st) +     sin(2.0*di)*sin(ra)*(sin(st))**2) !+mxx
            tldoco(iso)%myz = +m0*(cos(di)*cos(ra)*cos(    st) +     cos(2.0*di)*sin(ra)* sin(st))     !-mxz
            tldoco(iso)%mzz = +m0*sin(2.0*di)*sin(ra)                                                  !+mzz
   
            if ( (iso == 1) .and. ((tldoco(iso)%icur == 9) .or. (tldoco(iso)%icur == 10)) ) then
               itempo = 0
               open(unit=99,file=trim(cg_prefix)//".dcf")
               do while (.true.)
                  read(unit= 99,fmt=*,iostat=ios1) ctempo
                  if (ios1 /= 0) exit
                  itempo = itempo + 1
               enddo
               rewind(99)
               allocate(rg_dcsource_user_func(ig_ndt+1),dum3(ig_ndt+1),dum1(itempo),dum2(itempo))
               do i = 1,itempo
                  read(unit=99,fmt=*) dum1(i),dum2(i)
               enddo
               close(99)
               do i = 1,ig_ndt+1
                  dum3(i) = (i-1)*rg_dt
               enddo
               call interp_linear(1,itempo,dum1,dum2,ig_ndt+1,dum3,rg_dcsource_user_func)
               if (ig_myrank == 0) then
                  open(unit=90,file=trim(cg_prefix)//".dcf.interp")
                  do i = 1,ig_ndt+1
                     write(unit=90,fmt='(2(E14.7,1X))') dum3(i),rg_dcsource_user_func(i)
                  enddo
                  close(90)
               endif
            endif
         !
         !->find the closest element and gll point to the source iso
            call search_closest_hexa_gll(tldoco(iso)%x,tldoco(iso)%y,tldoco(iso)%z &
                                                     ,tldoco(iso)%dmin,tldoco(iso)%kgll,tldoco(iso)%lgll,tldoco(iso)%mgll,tldoco(iso)%iel)
         !
         !->check which cpu has the smallest dmin. this cpu will compute this source and the other won't
            call mpi_allgather(tldoco(iso)%dmin,1,mpi_real,dummy,1,mpi_real,mpi_comm_world,ios)
            min_loc = minloc(dummy(1:ig_ncpu))
            if ( (min_loc(1)-1) == ig_myrank) then
               tldoco(iso)%cpu  = ig_myrank
               ilndco = ilndco + 1
            else
               tldoco(iso)%cpu  = -1    
            endif
            if (ig_myrank == 0) then
               write(IG_LST_UNIT,'(a,i8,a,i0)') "double couple source ",iso," computed by cpu ",min_loc(1)-1
            endif
         enddo
         close(100)
         !
         !->transfer array tldoco which knows all the sources to array tg_dcsource that needs only to know the sources of cpu 'ig_myrank'
         if (ilndco > 0) then
            allocate(tg_dcsource(ilndco))
            j = 0
            do iso = 1,ig_ndcsource
               if (tldoco(iso)%cpu == ig_myrank) then
                  j = j + 1
                  tg_dcsource(j) = tldoco(iso)
               endif
            enddo
            ig_ndcsource = ilndco
         else
            ig_ndcsource = 0
         endif
         deallocate(tldoco)
         !
         !->cpu 0 gather the number of source (ig_ndcsource) of the other cpus
         call mpi_gather(ig_ndcsource,1,mpi_integer,ilrcpu,1,mpi_integer,0,mpi_comm_world,ios)
         !
         do iso = 1,ig_ndcsource 
            call compute_source_local_coordinate(tg_dcsource(iso))
            if (ig_myrank /= 0) call mpi_send(tg_dcsource(iso)%dmin,1,mpi_real,0,100,mpi_comm_world,ios)
         enddo
         !
         !->cpu 0 gathers dmin of all cpus and write the maximum dmin in *.lst
         if (ig_myrank == 0) then
            maxdmin = 0.0
            do i = 1,ig_ncpu
               if (i == 1) then
                  do iso = 1,ig_ndcsource
                     maxdmin = max(maxdmin,tg_dcsource(iso)%dmin)
                  enddo
               else
                  do iso = 1,ilrcpu(i)
                     call mpi_recv(rldmin,1,mpi_real,i-1,100,mpi_comm_world,statut,ios)
                     maxdmin = max(maxdmin,rldmin)
                  enddo 
               endif
            enddo
            write(IG_LST_UNIT,'(a,e14.7)') "maximum localisation error of all double couple point sources = ",maxdmin
            call flush(IG_LST_UNIT)
         endif
      else
         if (ig_myrank == 0) write(IG_LST_UNIT,'(/,a)') "no double couple point source computed"
         ig_ndcsource = 0
      endif

      return
!***********************************************************************************************************************************************************************************
   end subroutine init_double_couple_source
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine reads all single force point sources in file *.sfs; determines to which cpu belong single force point sources and
!!computes local coordinates @f$\xi,\eta,\zeta@f$ of single force point sources inside the hexahedron element it belongs.
!>@return mod_global_variables::tg_sfsource
!>@return mod_global_variables::ig_nsfsource
!>@return mod_global_variables::rg_sfsource_user_func
!***********************************************************************************************************************************************************************************
   subroutine init_single_force_source()
!***********************************************************************************************************************************************************************************
   
      use mpi
      
      use mod_global_variables, only :&
                                      IG_LST_UNIT&
                                     ,tg_sfsource&
                                     ,ig_nsfsource&
                                     ,ig_ncpu&
                                     ,ig_myrank&
                                     ,type_single_force_source&
                                     ,cg_prefix&
                                     ,rg_dt&
                                     ,ig_ndt&
                                     ,rg_sfsource_user_func&
                                     ,ig_hexa_gll_glonum
      
      use mod_source_function , only : interp_linear
    
      use mod_receiver        , only : search_closest_hexa_gll 
      
      implicit none
      
      real, allocatable, dimension(:)     :: dum1,dum2,dum3
      real, dimension(ig_ncpu)            :: dummy
      real                                :: maxdmin
      real                                :: rldmin
      
      integer, dimension(mpi_status_size) :: statut
      integer, dimension(1)               :: min_loc
      integer, dimension(ig_ncpu)         :: ilrcpu
      integer                             :: i
      integer                             :: j
      integer                             :: ipf
      integer                             :: ios
      integer                             :: ios1
      integer                             :: ilnpfo
      integer                             :: itempo
      
      character(len=1)                    :: ctempo
      
      type(type_single_force_source), allocatable, dimension(:) :: tlpofo

      !
      !
      !*****************************************************************
      !find element containing single force point source in cpu myrank
      !*****************************************************************
      
      ilnpfo = 0

      open(unit=100,file=trim(cg_prefix)//".sfs",status='old',iostat=ios)

      if (ios == 0) then

         if (ig_myrank == 0) write(IG_LST_UNIT,'(a)') " "

         read(unit=100,fmt='(i10)') ig_nsfsource

         allocate(tlpofo(ig_nsfsource))

         do ipf = 1,ig_nsfsource
         
           !initialize tlpofo
            tlpofo(ipf)%x          = 0.0
            tlpofo(ipf)%y          = 0.0
            tlpofo(ipf)%z          = 0.0
            tlpofo(ipf)%fac        = 0.0
            tlpofo(ipf)%rise_time  = 0.0
            tlpofo(ipf)%shift_time = 0.0
            tlpofo(ipf)%var1       = 0.0
            tlpofo(ipf)%dmin       = 0.0
            tlpofo(ipf)%icur       = 0.0
            tlpofo(ipf)%idir       = 0
            tlpofo(ipf)%cpu        = 0
            tlpofo(ipf)%iel        = 0
            tlpofo(ipf)%iequ       = 0
            tlpofo(ipf)%kgll       = 0
            tlpofo(ipf)%lgll       = 0
            tlpofo(ipf)%mgll       = 0
         
            read(unit=100,fmt=*) tlpofo(ipf)%x,tlpofo(ipf)%y,tlpofo(ipf)%z,tlpofo(ipf)%fac,tlpofo(ipf)%idir,tlpofo(ipf)%shift_time,tlpofo(ipf)%rise_time,tlpofo(ipf)%icur

            if ( (ipf == 1) .and. ((tlpofo(ipf)%icur == 9) .or. (tlpofo(ipf)%icur == 10)) ) then
               itempo = 0
               open(unit=99,file=trim(cg_prefix)//".pff")
               do while (.true.)
                  read(unit= 99,fmt=*,iostat=ios1) ctempo
                  if (ios1 /= 0) exit
                  itempo = itempo + 1
               enddo
               rewind(99)
               allocate(rg_sfsource_user_func(ig_ndt+1),dum3(ig_ndt+1),dum1(itempo),dum2(itempo))
               do i = 1,itempo
                  read(unit=99,fmt=*) dum1(i),dum2(i)
               enddo
               close(99)
               do i = 1,ig_ndt+1
                  dum3(i) = (i-1)*rg_dt
               enddo
               call interp_linear(1,itempo,dum1,dum2,ig_ndt+1,dum3,rg_sfsource_user_func)
               if (ig_myrank == 0) then
                  open(unit=90,file=trim(cg_prefix)//".pff.interp")
                  do i = 1,ig_ndt+1
                     write(unit=90,fmt='(2(E14.7,1X))') dum3(i),rg_sfsource_user_func(i)
                  enddo
                  close(90)
               endif
            endif
         !
         !---->find the closest element and gll point to the source ipf
            call search_closest_hexa_gll(tlpofo(ipf)%x,tlpofo(ipf)%y,tlpofo(ipf)%z &
                                              ,tlpofo(ipf)%dmin,tlpofo(ipf)%kgll,tlpofo(ipf)%lgll,tlpofo(ipf)%mgll,tlpofo(ipf)%iel)
         !
         !->check which cpu has the smallest dmin. this cpu will compute this source and the other won't
            call mpi_allgather(tlpofo(ipf)%dmin,1,mpi_real,dummy,1,mpi_real,mpi_comm_world,ios)
            min_loc = minloc(dummy(1:ig_ncpu))
            if ( (min_loc(1)-1) == ig_myrank) then
               tlpofo(ipf)%cpu  = ig_myrank
               ilnpfo = ilnpfo + 1
            else
               tlpofo(ipf)%cpu  = -1    
            endif
            if (ig_myrank == 0) then
               write(IG_LST_UNIT,'(a,i8,a,i0)') "single force point source ",ipf," computed by cpu ",min_loc(1)-1
            endif
         enddo
         close(100)
         !
         !->transfer array tlpofo which knows all the sources to array tg_sfsource that needs only to know the sources of cpu 'ig_myrank'
         if (ilnpfo > 0) then
            allocate(tg_sfsource(ilnpfo))
            j = 0
            do ipf = 1,ig_nsfsource
               if (tlpofo(ipf)%cpu == ig_myrank) then
                  j = j + 1
                  tg_sfsource(j) = tlpofo(ipf)
               endif
            enddo
            ig_nsfsource = ilnpfo
         else
            ig_nsfsource = 0
         endif
         deallocate(tlpofo)
         !
         !->cpu 0 gather the number of single force point source (ig_nsfsource) of the other cpus
         call mpi_gather(ig_nsfsource,1,mpi_integer,ilrcpu,1,mpi_integer,0,mpi_comm_world,ios)
         !
         do ipf = 1,ig_nsfsource 
         !  call pfoinf(tg_sfsource(ipf))
            tg_sfsource(ipf)%iequ = ig_hexa_gll_glonum(tg_sfsource(ipf)%mgll,tg_sfsource(ipf)%lgll,tg_sfsource(ipf)%kgll,tg_sfsource(ipf)%iel)
            if (ig_myrank /= 0) call mpi_send(tg_sfsource(ipf)%dmin,1,mpi_real,0,100,mpi_comm_world,ios)
         enddo
         !
         !->cpu 0 gathers dmin of all cpus and write the maximum dmin in *.lst
         if (ig_myrank == 0) then
            maxdmin = 0.0
            do i = 1,ig_ncpu
               if (i == 1) then
                  do ipf = 1,ig_nsfsource
                     maxdmin = max(maxdmin,tg_sfsource(ipf)%dmin)
                  enddo
               else
                  do ipf = 1,ilrcpu(i)
                     call mpi_recv(rldmin,1,mpi_real,i-1,100,mpi_comm_world,statut,ios)
                     maxdmin = max(maxdmin,rldmin)
                  enddo 
               endif
            enddo
            write(IG_LST_UNIT,'(a,e14.7)') "maximum localisation error of all single force point sources = ",maxdmin
            call flush(IG_LST_UNIT)
         endif
      else
         if (ig_myrank == 0) write(IG_LST_UNIT,'(/,a)') "no single force point source computed"
         ig_nsfsource = 0
      endif
      
      return
!***********************************************************************************************************************************************************************************
   end subroutine init_single_force_source
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine solves a nonlinear system by a gradient method to compute local coordinates @f$\xi,\eta,\zeta@f$ of a point source inside the hexahedron element it belongs.
!!The starting point of the gradient method is the closest GLL node.
!!Local coordinates are stored in mod_global_variables::tg_dcsource.
!>@param tlsrc : data type that contains information about point sources
!***********************************************************************************************************************************************************************************
   subroutine compute_source_local_coordinate(tlsrc)
!***********************************************************************************************************************************************************************************
      
      use mpi
   
      use mod_global_variables, only :&
                                      rg_gll_abscissa&
                                     ,ig_myrank&
                                     ,type_double_couple_source

      use mod_jacobian        , only : compute_hexa_jacobian
   
      use mod_coordinate      , only : compute_hexa_point_coord
   
      
      implicit none
   
      type(type_double_couple_source), intent(inout) :: tlsrc
      
      real :: eps
      real :: dmin
      real :: xisol
      real :: etsol
      real :: zesol
      real :: newx
      real :: newy
      real :: newz
      real :: dxidx
      real :: dxidy
      real :: dxidz
      real :: detdx
      real :: detdy
      real :: detdz
      real :: dzedx
      real :: dzedy
      real :: dzedz
      real :: dx
      real :: dy
      real :: dz
      real :: dxi
      real :: det
      real :: dze
      real :: deter

      integer :: ihexa
      integer :: iter

      integer, parameter :: ITER_MAX=100
                              
      !
      !->initialize value
      xisol = rg_gll_abscissa(tlsrc%mgll)
      etsol = rg_gll_abscissa(tlsrc%lgll)
      zesol = rg_gll_abscissa(tlsrc%kgll)
      newx  = tlsrc%x
      newy  = tlsrc%y
      newz  = tlsrc%z
      dmin  = tlsrc%dmin
      ihexa = tlsrc%iel  
      eps  = 1.0e-2
      iter = 0
   
      !      
      !->solve the nonlinear system
      do while(dmin > EPS)
   
         call compute_hexa_jacobian(ihexa,xisol,etsol,zesol,dxidx,dxidy,dxidz,detdx,detdy,detdz,dzedx,dzedy,dzedz,deter)
   
         call compute_hexa_point_coord(ihexa,xisol,etsol,zesol,newx,newy,newz)
   
         dx    = tlsrc%x - newx
         dy    = tlsrc%y - newy
         dz    = tlsrc%z - newz
         dxi   = dxidx*dx + dxidy*dy + dxidz*dz
         det   = detdx*dx + detdy*dy + detdz*dz
         dze   = dzedx*dx + dzedy*dy + dzedz*dz
         xisol = xisol + dxi
         etsol = etsol + det
         zesol = zesol + dze
         dmin  = sqrt( (newx-tlsrc%x)**2 + (newy-tlsrc%y)**2 + (newz-tlsrc%z)**2 )
   
         iter  = iter + 1
   
         if (mod(iter,ITER_MAX) == 0) then
            eps = eps*2.0
         endif
   
      enddo
   
      tlsrc%xi   = xisol
      tlsrc%et   = etsol
      tlsrc%ze   = zesol
      tlsrc%dmin = dmin 
      
      return
   
!***********************************************************************************************************************************************************************************
   end subroutine compute_source_local_coordinate
!***********************************************************************************************************************************************************************************


end module mod_source
