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
!!This file contains a module to initialize the medium where wave propagation takes place.

!>@brief
!!This module contains subroutines to initialize medium for hexahedron and quadrangle elements. It also contains subroutines to write medium in binary VTK xml format readable by ParaView.
module mod_init_medium
   
   use mpi

   implicit none

   private

   public  :: init_hexa_medium
   public  :: init_quadp_medium
   private :: init_medium_from_mat_file
   private :: write_medium_vtk_cell_xml
   private :: write_medium_vtk_node_xml
   private :: write_medium_collection
   
   contains


!
!
!>@brief
!!This subroutine fills medium properties (elastic or viscoelastic) at GLL nodes of hexahedron elements based on material number found in file *.inp generated by CUBIT mesh generator.
!!If files *.cpu.*.mat are present in the simulation's directory, then material number of hexahedron elements are read from files *.cpu*.mat
!>@return mod_global_variables::rg_hexa_gll_rho    
!>@return mod_global_variables::rg_hexa_gll_rhovs2 
!>@return mod_global_variables::rg_hexa_gll_rhovp2
!>@return mod_global_variables::rg_hexa_gll_wkqs
!>@return mod_global_variables::rg_hexa_gll_wkqp
!>@return mod_global_variables::rg_mem_var_exp
!>@return mod_global_variables::tg_elastic_material
!>@return mod_global_variables::tg_viscoelastic_material
!***********************************************************************************************************************************************************************************
      subroutine init_hexa_medium()
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only :&
                                      ig_nmaterial&
                                     ,tg_elastic_material&
                                     ,ig_myrank&
                                     ,cg_prefix&
                                     ,cg_myrank&
                                     ,IG_NGLL&
                                     ,ig_material_type&
                                     ,ig_hexa_material_number&
                                     ,LG_VISCO&
                                     ,LG_OUTPUT_MEDIUM_VTK&
                                     ,ig_nhexa&
                                     ,rg_hexa_gll_rho&
                                     ,rg_hexa_gll_rhovs2&
                                     ,rg_hexa_gll_rhovp2&
                                     ,rg_hexa_gll_wkqs&
                                     ,rg_hexa_gll_wkqp&
                                     ,rg_hexa_gll_ksixx&
                                     ,rg_hexa_gll_ksiyy&
                                     ,rg_hexa_gll_ksizz&
                                     ,rg_hexa_gll_ksixy&
                                     ,rg_hexa_gll_ksixz&
                                     ,rg_hexa_gll_ksiyz&
                                     ,rg_mem_var_exp&
                                     ,tg_viscoelastic_material&
                                     ,RG_RELAX_COEFF&
                                     ,IG_NRELAX&
                                     ,RG_PI&
                                     ,rg_dt&
                                     ,error_stop&
                                     ,IG_LST_UNIT

      use mod_init_memory

      implicit none
   
      complex            :: ctmpqs
      complex            :: ctmpqp
   
      real               :: dtmpqs
      real               :: dtmpqp
   
      integer            :: ios
      integer            :: iel
      integer            :: imat
      integer            :: imem_var
      integer            :: k
      integer            :: l
      integer            :: m
      integer            :: n
   
      logical            :: is_file_exist
   
      character(len=255) :: info


      !
      !
      !************************************************************************************************
      !initialize memory
      !************************************************************************************************

      ios = init_array_real(rg_hexa_gll_rho,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_rho") !deallocated in subroutine init_mass_matrix.f90

      ios = init_array_real(rg_hexa_gll_rhovs2,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_rhovs2")

      ios = init_array_real(rg_hexa_gll_rhovp2,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,"rg_hexa_gll_rhovp2")


      if (LG_VISCO) then

         ios = init_array_real(rg_hexa_gll_ksixx,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,IG_NRELAX,"rg_hexa_gll_ksixx")

         ios = init_array_real(rg_hexa_gll_ksiyy,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,IG_NRELAX,"rg_hexa_gll_ksiyy")

         ios = init_array_real(rg_hexa_gll_ksizz,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,IG_NRELAX,"rg_hexa_gll_ksizz")

         ios = init_array_real(rg_hexa_gll_ksixy,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,IG_NRELAX,"rg_hexa_gll_ksixy")

         ios = init_array_real(rg_hexa_gll_ksixz,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,IG_NRELAX,"rg_hexa_gll_ksixz")

         ios = init_array_real(rg_hexa_gll_ksiyz,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,IG_NRELAX,"rg_hexa_gll_ksiyz")

         ios = init_array_real(rg_hexa_gll_wkqs,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,IG_NRELAX,"rg_hexa_gll_wkqs")

         ios = init_array_real(rg_hexa_gll_wkqp,ig_nhexa,IG_NGLL,IG_NGLL,IG_NGLL,IG_NRELAX,"rg_hexa_gll_wkqp")
        
!
!------->precompute exponential used for computing internal forces
         do imem_var = 1,IG_NRELAX
            rg_mem_var_exp(imem_var) = exp(-rg_dt/RG_RELAX_COEFF(imem_var,1))
         enddo
   
      endif 


      !
      !
      !***********************************************************************************************************
      !->if a geomodeler file *.vox or *.mat exists, then the medium (ig_hexa_material_number) is read from such files
      !***********************************************************************************************************
      is_file_exist = .false.

      inquire(file=trim(cg_prefix)//".cpu."//trim(cg_myrank)//".mat", exist=is_file_exist)
      if (is_file_exist) then
         if (ig_myrank == 0) then
            write(IG_LST_UNIT,'(a)') " --> medium generated from geomodeler's material (file *.mat)"
         endif
         call init_medium_from_mat_file(ig_hexa_material_number)
      endif
 
      if ( .not.(is_file_exist) .and. (ig_myrank == 0) ) then
         write(IG_LST_UNIT,'(a)') " --> medium generated from cubit block's information"
      endif


      !
      !
      !************************************************************************************************
      !->check if the medium (ig_hexa_material_number) is entirely defined
      !************************************************************************************************
      do iel = 1,ig_nhexa

         imat = ig_hexa_material_number(iel)

         if ( (imat <= 0) .or. (imat > ig_nmaterial) ) then
            write(info,'(a)') "error in subroutine init_hexa_medium: undefined material for hexa"
            call error_stop(info)
         endif

      enddo


      !
      !
      !************************************************************************************************
      !Compute elastic and viscoelastic parameter pour media_type = 0
      !************************************************************************************************
      do imat = 1,ig_nmaterial
!
!
!------->elastic
         tg_elastic_material(imat)%lam2 = tg_elastic_material(imat)%dens*tg_elastic_material(imat)%svel**2
         tg_elastic_material(imat)%lam1 = tg_elastic_material(imat)%dens*tg_elastic_material(imat)%pvel**2 - 2.0*tg_elastic_material(imat)%lam2
         tg_elastic_material(imat)%pois = (tg_elastic_material(imat)%pvel**2 - 2.0*tg_elastic_material(imat)%svel**2)/(2.0*(tg_elastic_material(imat)%pvel**2-tg_elastic_material(imat)%svel**2))
         tg_elastic_material(imat)%youn = 2.0*tg_elastic_material(imat)%lam2*(1.0+tg_elastic_material(imat)%pois)
         tg_elastic_material(imat)%bulk = tg_elastic_material(imat)%lam1+2.0/3.0*tg_elastic_material(imat)%lam2
         tg_elastic_material(imat)%rvel = (0.862+1.14*tg_elastic_material(imat)%pois)/(1.0+tg_elastic_material(imat)%pois)*tg_elastic_material(imat)%svel
         tg_elastic_material(imat)%pwmo = tg_elastic_material(imat)%lam1 + 2.0*tg_elastic_material(imat)%lam2
   
         if (ig_myrank == 0) then
            if (ig_material_type(imat) == 1) write(IG_LST_UNIT,'(" ",/,a,i0,a)') "material ",imat," is elastic"
            if (ig_material_type(imat) == 2) write(IG_LST_UNIT,'(" ",/,a,i0,a)') "material ",imat," is viscoelactic"
            write(IG_LST_UNIT,'(a,f15.8  ,a)') "s-wave velocity           = ",tg_elastic_material(imat)%svel," m/s"
            write(IG_LST_UNIT,'(a,f15.8  ,a)') "p-wave velocity           = ",tg_elastic_material(imat)%pvel," m/s"
            write(IG_LST_UNIT,'(a,f15.8  ,a)') "poisson ratio             = ",tg_elastic_material(imat)%pois," "
            write(IG_LST_UNIT,'(a,f15.8  ,a)') "density                   = ",tg_elastic_material(imat)%dens," kg/m3"
            write(IG_LST_UNIT,'(a,e15.7e3,a)') "lambda (first lame coef.) = ",tg_elastic_material(imat)%lam1," Pa"
            write(IG_LST_UNIT,'(a,e15.7e3,a)') "s-wave modulus            = ",tg_elastic_material(imat)%lam2," Pa"
            write(IG_LST_UNIT,'(a,e15.7e3,a)') "p-wave modulus            = ",tg_elastic_material(imat)%pwmo," Pa"
            write(IG_LST_UNIT,'(a,e15.7e3,a)') "bulk-modulus              = ",tg_elastic_material(imat)%bulk," Pa"
         endif
!
!
!------->viscoelastic
         if (LG_VISCO) then
   
               allocate(tg_viscoelastic_material(imat)%wkqs(IG_NRELAX),tg_viscoelastic_material(imat)%wkqp(IG_NRELAX))
   
               dtmpqs = (3.071+1.433*tg_viscoelastic_material(imat)%qs**(-1.158)*log(tg_viscoelastic_material(imat)%qs/5.0))/(1.0+0.415*tg_viscoelastic_material(imat)%qs)
               dtmpqp = (3.071+1.433*tg_viscoelastic_material(imat)%qp**(-1.158)*log(tg_viscoelastic_material(imat)%qp/5.0))/(1.0+0.415*tg_viscoelastic_material(imat)%qp)
               do k = 1,IG_NRELAX
                  tg_viscoelastic_material(imat)%wkqs(k) = dtmpqs*(dtmpqs*RG_RELAX_COEFF(k,2)+RG_RELAX_COEFF(k,3))
                  tg_viscoelastic_material(imat)%wkqp(k) = dtmpqp*(dtmpqp*RG_RELAX_COEFF(k,2)+RG_RELAX_COEFF(k,3))
               enddo
   
               ctmpqs = cmplx(0.0,0.0)
               ctmpqp = cmplx(0.0,0.0)
               do k = 1,IG_NRELAX
                  ctmpqs = ctmpqs + cmplx(tg_viscoelastic_material(imat)%wkqs(k),0.0)/(1.0 + cmplx(0.0,1.0)*2.0*RG_PI*tg_viscoelastic_material(imat)%freq*RG_RELAX_COEFF(k,1))
                  ctmpqp = ctmpqp + cmplx(tg_viscoelastic_material(imat)%wkqp(k),0.0)/(1.0 + cmplx(0.0,1.0)*2.0*RG_PI*tg_viscoelastic_material(imat)%freq*RG_RELAX_COEFF(k,1))
               enddo
   
               tg_viscoelastic_material(imat)%uswm = tg_elastic_material(imat)%lam2/abs(1.0-ctmpqs)
               tg_viscoelastic_material(imat)%upwm = tg_elastic_material(imat)%pwmo/abs(1.0-ctmpqp)
               if (ig_myrank == 0) then
                  write(IG_LST_UNIT,'(a,f15.8  ,a)') "frequency for viscosity   = ",tg_viscoelastic_material(imat)%freq," hz"
                  write(IG_LST_UNIT,'(a,f15.8  ,a)') "s-wave quality factor     = ",tg_viscoelastic_material(imat)%qs  ," "
                  write(IG_LST_UNIT,'(a,f15.8  ,a)') "p-wave quality factor     = ",tg_viscoelastic_material(imat)%qp  ," "
                  write(IG_LST_UNIT,'(a,e15.7e3,a)') "unrelaxed s-wave modulus  = ",tg_viscoelastic_material(imat)%uswm," pa"
                  write(IG_LST_UNIT,'(a,e15.7e3,a)') "unrelaxed p-wave modulus  = ",tg_viscoelastic_material(imat)%upwm," pa"
               endif
   
         endif !endif of material type viscoelastic
   
      enddo !enddo of loop on material


      !
      !
      !************************************************************************************************
      !Distribute material properties to gll nodes
      !************************************************************************************************
!
!
!------->elastic
      if (.not.LG_VISCO) then
         do iel = 1,ig_nhexa
            do k = 1,IG_NGLL        !zeta
               do l = 1,IG_NGLL     !eta
                  do m = 1,IG_NGLL  !xi
         
                     imat                          = ig_hexa_material_number(iel)
                     rg_hexa_gll_rho   (m,l,k,iel) = tg_elastic_material(imat)%dens
                     rg_hexa_gll_rhovs2(m,l,k,iel) = tg_elastic_material(imat)%dens*tg_elastic_material(imat)%svel**2
                     rg_hexa_gll_rhovp2(m,l,k,iel) = tg_elastic_material(imat)%dens*tg_elastic_material(imat)%pvel**2
    
                  enddo
               enddo
            enddo
         enddo
!
!
!------->viscoelastic
      else
         do iel = 1,ig_nhexa
            do k = 1,IG_NGLL        !zeta
               do l = 1,IG_NGLL     !eta
                  do m = 1,IG_NGLL  !xi
         
                     imat                          = ig_hexa_material_number(iel)
                     rg_hexa_gll_rho   (m,l,k,iel) = tg_elastic_material(imat)%dens
                     rg_hexa_gll_rhovs2(m,l,k,iel) = tg_viscoelastic_material(imat)%uswm
                     rg_hexa_gll_rhovp2(m,l,k,iel) = tg_viscoelastic_material(imat)%upwm
   
                     do n = 1,IG_NRELAX
                        rg_hexa_gll_wkqs(n,m,l,k,iel) = tg_viscoelastic_material(imat)%wkqs(n)
                        rg_hexa_gll_wkqp(n,m,l,k,iel) = tg_viscoelastic_material(imat)%wkqp(n)
                     enddo
         
                  enddo
               enddo
            enddo
         enddo
      endif

      !
      !
      !************************************************************************************************
      !write material in VTK xml format
      !************************************************************************************************

      if (LG_OUTPUT_MEDIUM_VTK) call write_medium_vtk_cell_xml()
   
      return
!***********************************************************************************************************************************************************************************
      end subroutine init_hexa_medium
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine fills medium properties at GLL nodes of paraxial quadrangle elements (i.e., absorbing boundary) based on material number of hexahedron elements it belongs.
!>@return mod_global_variables::rg_quadp_gll_rhovs 
!>@return mod_global_variables::rg_quadp_gll_rhovp 
!***********************************************************************************************************************************************************************************
      subroutine init_quadp_medium()
!***********************************************************************************************************************************************************************************
      use mpi
      use mod_global_variables, only :&
                                      IG_NGLL&
                                     ,ig_myrank&
                                     ,ig_nquad_parax&
                                     ,ig_quadp_neighbor_hexa&
                                     ,ig_quadp_neighbor_hexaface&
                                     ,rg_quadp_gll_rhovs&
                                     ,rg_quadp_gll_rhovp&
                                     ,rg_hexa_gll_rho&
                                     ,rg_hexa_gll_rhovs2&
                                     ,rg_hexa_gll_rhovp2&
                                     ,error_stop&
                                     ,lg_boundary_absorption

      use mod_init_memory

      implicit none
   
      integer            :: ios
      integer            :: i
      integer            :: j
      integer            :: k
      integer            :: iquad
      integer            :: ihexa
      integer            :: iface

      character(len=255) :: info

      !
      !See numbering convention in subroutine propagate_gll_nodes_quad of module_mesh_gll_num.f90
      !

      if (ig_nquad_parax > 0 .and. lg_boundary_absorption) then

         ios = init_array_real(rg_quadp_gll_rhovs,ig_nquad_parax,IG_NGLL,IG_NGLL,"rg_quadp_gll_rhovs")

         ios = init_array_real(rg_quadp_gll_rhovp,ig_nquad_parax,IG_NGLL,IG_NGLL,"rg_quadp_gll_rhovp")

         do iquad = 1,ig_nquad_parax
    
            ihexa = ig_quadp_neighbor_hexa    (iquad)
            iface = ig_quadp_neighbor_hexaface(iquad)
    
            select case(iface)
               case(1)
                  do j = 1,IG_NGLL
                     do i = 1,IG_NGLL
                        rg_quadp_gll_rhovs(i,j,iquad) = rg_hexa_gll_rho(i,j,1,ihexa)*sqrt(rg_hexa_gll_rhovs2(i,j,1,ihexa)/rg_hexa_gll_rho(i,j,1,ihexa))
                        rg_quadp_gll_rhovp(i,j,iquad) = rg_hexa_gll_rho(i,j,1,ihexa)*sqrt(rg_hexa_gll_rhovp2(i,j,1,ihexa)/rg_hexa_gll_rho(i,j,1,ihexa))
                     enddo
                  enddo

               case(2)
                  do k = 1,IG_NGLL
                     do i = 1,IG_NGLL
                        rg_quadp_gll_rhovs(k,i,iquad) = rg_hexa_gll_rho(i,1,k,ihexa)*sqrt(rg_hexa_gll_rhovs2(i,1,k,ihexa)/rg_hexa_gll_rho(i,1,k,ihexa))
                        rg_quadp_gll_rhovp(k,i,iquad) = rg_hexa_gll_rho(i,1,k,ihexa)*sqrt(rg_hexa_gll_rhovp2(i,1,k,ihexa)/rg_hexa_gll_rho(i,1,k,ihexa))
                     enddo
                  enddo

               case(3)
                  do k = 1,IG_NGLL
                     do j = 1,IG_NGLL
                        rg_quadp_gll_rhovs(k,j,iquad) = rg_hexa_gll_rho(IG_NGLL,j,k,ihexa)*sqrt(rg_hexa_gll_rhovs2(IG_NGLL,j,k,ihexa)/rg_hexa_gll_rho(IG_NGLL,j,k,ihexa))
                        rg_quadp_gll_rhovp(k,j,iquad) = rg_hexa_gll_rho(IG_NGLL,j,k,ihexa)*sqrt(rg_hexa_gll_rhovp2(IG_NGLL,j,k,ihexa)/rg_hexa_gll_rho(IG_NGLL,j,k,ihexa))
                     enddo
                  enddo

               case(4)
                  do k = 1,IG_NGLL
                     do i = 1,IG_NGLL
                        rg_quadp_gll_rhovs(i,k,iquad) = rg_hexa_gll_rho(i,IG_NGLL,k,ihexa)*sqrt(rg_hexa_gll_rhovs2(i,IG_NGLL,k,ihexa)/rg_hexa_gll_rho(i,IG_NGLL,k,ihexa))
                        rg_quadp_gll_rhovp(i,k,iquad) = rg_hexa_gll_rho(i,IG_NGLL,k,ihexa)*sqrt(rg_hexa_gll_rhovp2(i,IG_NGLL,k,ihexa)/rg_hexa_gll_rho(i,IG_NGLL,k,ihexa))
                     enddo
                  enddo

               case(5)
                  do k = 1,IG_NGLL
                     do j = 1,IG_NGLL
                        rg_quadp_gll_rhovs(j,k,iquad) = rg_hexa_gll_rho(1,j,k,ihexa)*sqrt(rg_hexa_gll_rhovs2(1,j,k,ihexa)/rg_hexa_gll_rho(1,j,k,ihexa))
                        rg_quadp_gll_rhovp(j,k,iquad) = rg_hexa_gll_rho(1,j,k,ihexa)*sqrt(rg_hexa_gll_rhovp2(1,j,k,ihexa)/rg_hexa_gll_rho(1,j,k,ihexa))
                     enddo
                  enddo

               case(6)
                  do j = 1,IG_NGLL
                     do i = 1,IG_NGLL
                        rg_quadp_gll_rhovs(j,i,iquad) = rg_hexa_gll_rho(i,j,IG_NGLL,ihexa)*sqrt(rg_hexa_gll_rhovs2(i,j,IG_NGLL,ihexa)/rg_hexa_gll_rho(i,j,IG_NGLL,ihexa))
                        rg_quadp_gll_rhovp(j,i,iquad) = rg_hexa_gll_rho(i,j,IG_NGLL,ihexa)*sqrt(rg_hexa_gll_rhovp2(i,j,IG_NGLL,ihexa)/rg_hexa_gll_rho(i,j,IG_NGLL,ihexa))
                     enddo
                  enddo

               case default
                  write(info,'(a)') "error in subroutine init_quadp_medium: invalid face number"
                  call error_stop(info)

            end select
    
         enddo
   
      endif
      
      if (allocated(ig_quadp_neighbor_hexa    )) deallocate(ig_quadp_neighbor_hexa)
      if (allocated(ig_quadp_neighbor_hexaface)) deallocate(ig_quadp_neighbor_hexaface)
      
      return
!***********************************************************************************************************************************************************************************
      end subroutine init_quadp_medium
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine fills medium properties at GLL nodes of hexahedron elements based on material number in files *.cpu*.mat
!>@param il_mat_num_of_hexa : array filled with material number of hexahedron elements in cpu myrank.
!***********************************************************************************************************************************************************************************
      subroutine init_medium_from_mat_file(il_mat_num_of_hexa)
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only :&
                                      ig_nhexa&
                                     ,cg_myrank&
                                     ,cg_prefix
   
      implicit none
   
      integer, intent(out)         :: il_mat_num_of_hexa(:)
      integer                      :: ihexa
   !
   !->fill the array 'il_mat_num_of_hexa' using the file *.mat specific for myrank
      open(unit=90,file=trim(cg_prefix)//".cpu."//trim(cg_myrank)//".mat")
      do ihexa = 1,ig_nhexa
         read(unit=90,fmt=*) il_mat_num_of_hexa(ihexa)
      enddo
      close(90)
   
      return
!***********************************************************************************************************************************************************************************
      end subroutine init_medium_from_mat_file
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine writes S and P-wave velocities of the physical domain of cpu myrank in UnstructuredGrid binary vtk_cell xml file readable by ParaView.
!>@return files *.vtu
!***********************************************************************************************************************************************************************************
      subroutine write_medium_vtk_cell_xml()
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only :&
                                      ig_hexa_material_number&
                                     ,ig_nhexa&
                                     ,ig_hexa_nnode&
                                     ,ig_hexa_gnode_glonum&
                                     ,ig_mesh_nnode&
                                     ,rg_gnode_x&
                                     ,rg_gnode_y&
                                     ,rg_gnode_z&
                                     ,tg_elastic_material&
                                     ,cg_prefix&
                                     ,ig_myrank&
                                     ,cg_myrank

      use mod_vtk_io

      implicit none

      real           , dimension(ig_nhexa)               :: var

      integer(kind=1), dimension(ig_nhexa)               :: cell_type
      integer        , dimension(ig_nhexa)               :: offset
      integer        , dimension(ig_hexa_nnode*ig_nhexa) :: connectivity
      integer                                            :: ihexa
      integer                                            :: inode
      integer                                            :: imat
      integer                                            :: icon
      integer                                            :: ios
                                                         
      character(len=255)                                 :: fname

!
!---->file name to store the medium of each cpu
      fname = trim(cg_prefix)//'.medium.cpu.'//trim(cg_myrank)//".vtu"

!
!---->fill the array cell_type with 12 --> hexahedron 8 nodes
      do ihexa = 1,ig_nhexa
         cell_type(ihexa) = 12
      enddo

!
!---->fill the array offset : sum of vertices over all hexahedra
      do ihexa = 1,ig_nhexa
         offset(ihexa) = ihexa*ig_hexa_nnode
      enddo

!
!---->fill the array connectivity
      icon = 0
      do ihexa = 1,ig_nhexa
         do inode = 1,ig_hexa_nnode

            icon               = icon + 1
            connectivity(icon) = ig_hexa_gnode_glonum(inode,ihexa) - 1

         enddo
      enddo

      ios = VTK_INI_XML(                                   &
                        output_format = 'BINARY'           &
                       ,filename      = fname              &
                       ,mesh_topology = 'UnstructuredGrid' )
                                                           
      ios = VTK_GEO_XML(                   &
                        NN = ig_mesh_nnode &
                       ,NC = ig_nhexa      &
                       ,X  = rg_gnode_x    &
                       ,Y  = rg_gnode_y    &
                       ,Z  = rg_gnode_z    )
                                                           
      ios = VTK_CON_XML(                         &
                        NC        = ig_nhexa     &
                       ,connect   = connectivity &
                       ,offset    = offset       &
                       ,cell_type = cell_type    )
    
      ios = VTK_DAT_XML(                          &
                        var_location     = 'CELL' &
                       ,var_block_action = 'OPEN' )

!
!---->fill the array var with P-wave velocity value
      do ihexa = 1,ig_nhexa

         imat         = ig_hexa_material_number(ihexa)
         var(ihexa)   = tg_elastic_material(imat)%pvel

      enddo

      ios = VTK_VAR_XML(                   &
                        NC_NN   = ig_nhexa &
                       ,varname = 'Vp'     &
                       ,var     = var      &
                                           )
!
!---->fill the array var with S-wave velocity value
      do ihexa = 1,ig_nhexa
         imat         = ig_hexa_material_number(ihexa)
         var(ihexa)   = tg_elastic_material(imat)%svel
      enddo

      ios = VTK_VAR_XML(                   &
                        NC_NN   = ig_nhexa &
                       ,varname = 'Vs'     &
                       ,var     = var      &
                                           )

      ios = VTK_DAT_XML(                           &
                        var_location     = 'CELL'  &
                       ,var_block_action = 'CLOSE' &
                                                   )

      ios = VTK_GEO_XML()

      ios = VTK_END_XML()

      if (ig_myrank == 0) call write_medium_collection()

      return
!***********************************************************************************************************************************************************************************
      end subroutine write_medium_vtk_cell_xml
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine writes S and P-wave velocities of the physical domain of cpu myrank in UnstructuredGrid binary vtk_node xml file readable by ParaView.
!>@return files *.vtu
!***********************************************************************************************************************************************************************************
      subroutine write_medium_vtk_node_xml()
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only :&
                                      ig_hexa_material_number&
                                     ,ig_nhexa&
                                     ,ig_hexa_nnode&
                                     ,ig_hexa_gnode_glonum&
                                     ,ig_mesh_nnode&
                                     ,ig_hexa_node2gll&
                                     ,rg_gnode_x&
                                     ,rg_gnode_y&
                                     ,rg_gnode_z&
                                     ,rg_hexa_gll_rho&
                                     ,rg_hexa_gll_rhovs2&
                                     ,cg_prefix&
                                     ,ig_myrank&
                                     ,cg_myrank

      use mod_vtk_io

      implicit none

      real           , dimension(ig_mesh_nnode)          :: var
      real                                               :: myvar

      integer(kind=1), dimension(ig_nhexa)               :: cell_type
      integer        , dimension(ig_nhexa)               :: offset
      integer        , dimension(ig_hexa_nnode*ig_nhexa) :: connectivity
      integer                                            :: ihexa
      integer                                            :: inode
      integer                                            :: global_gnode
      integer                                            :: icon
      integer                                            :: ios
      integer                                            :: k
      integer                                            :: l
      integer                                            :: m
                                                         
      character(len=255)                                 :: fname

!
!---->file name to store the medium of each cpu
      fname = trim(cg_prefix)//'.medium.cpu.'//trim(cg_myrank)//".vtu"

!
!---->fill the array cell_type with 12 --> hexahedron 8 nodes
      do ihexa = 1,ig_nhexa
         cell_type(ihexa) = 12
      enddo

!
!---->fill the array offset : sum of vertices over all hexahedra
      do ihexa = 1,ig_nhexa
         offset(ihexa) = ihexa*ig_hexa_nnode
      enddo

!
!---->fill the array connectivity
      icon = 0
      do ihexa = 1,ig_nhexa
         do inode = 1,ig_hexa_nnode

            icon               = icon + 1
            connectivity(icon) = ig_hexa_gnode_glonum(inode,ihexa) - 1

         enddo
      enddo

!
!---->initialize VTK xml file
      ios = VTK_INI_XML(                                   &
                        output_format = 'BINARY'           &
                       ,filename      = fname              &
                       ,mesh_topology = 'UnstructuredGrid' )
                                                           
!
!---->geometry information for VTK xml file
      ios = VTK_GEO_XML(                   &
                        NN = ig_mesh_nnode &
                       ,NC = ig_nhexa      &
                       ,X  = rg_gnode_x    &
                       ,Y  = rg_gnode_y    &
                       ,Z  = rg_gnode_z    &
                                           )
                                                           
!
!---->connectivity information for VTK xml file
      ios = VTK_CON_XML(                         &
                        NC        = ig_nhexa     &
                       ,connect   = connectivity &
                       ,offset    = offset       &
                       ,cell_type = cell_type    &
                                                 )
    
!
!---->data storage type for VTK xml file : cell or node
      ios = VTK_DAT_XML(                          &
                        var_location     = 'NODE' &
                       ,var_block_action = 'OPEN' &
                                                  )

!
!---->fill the array var with gll_node value

      do ihexa = 1,ig_nhexa

         do inode = 1,ig_hexa_nnode

            k = ig_hexa_node2gll(1,inode)
            l = ig_hexa_node2gll(2,inode)
            m = ig_hexa_node2gll(3,inode)

            global_gnode  = ig_hexa_gnode_glonum(inode,ihexa)
            myvar         = sqrt(rg_hexa_gll_rhovs2(m,l,k,ihexa)/rg_hexa_gll_rho(m,l,k,ihexa))

            var(global_gnode) = myvar

         enddo

      enddo

!
!---->write array var in the xml VTK file
      ios = VTK_VAR_XML(                        &
                        NC_NN   = ig_mesh_nnode &
                       ,varname = 'Vs'          &
                       ,var     = var           &
                                                )
!
!---->finalize xml VTK file
      ios = VTK_DAT_XML(                           &
                        var_location     = 'NODE'  &
                       ,var_block_action = 'CLOSE' &
                                                   )
      ios = VTK_GEO_XML()

      ios = VTK_END_XML()

      if (ig_myrank == 0) call write_medium_collection()

      return
!***********************************************************************************************************************************************************************************
      end subroutine write_medium_vtk_node_xml
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine writes a VTK xml collection file *.pvd which contains all *.vtu files generated by \ref write_medium_vtk_cell_xml or \ref write_medium_vtk_node_xml.
!>@return files *.medium.collection.pvd
!***********************************************************************************************************************************************************************************
      subroutine write_medium_collection()
!***********************************************************************************************************************************************************************************

         use mod_global_variables, only :&
                                         get_newunit&
                                        ,cg_prefix&
                                        ,ig_ncpu


         implicit none 
         
         integer                      :: icpu
         integer                      :: myunit
         
         character(len=255)           :: fname
         character(len=6  )           :: cpart
         
         open(unit=get_newunit(myunit),file=trim(cg_prefix)//".medium.collection.pvd")
         
         write(unit=myunit,fmt='(a)') "<?xml version=""1.0""?>"
         write(unit=myunit,fmt='(a)') "<VTKFile type=""Collection"" version=""0.1"" byte_order=""BigEndian"">"
         write(unit=myunit,fmt='(a)') "  <Collection>"
         
         do icpu = 1,ig_ncpu
         
            write(cpart,'(I6.6)') icpu-1
         
            fname = trim(cg_prefix)//".medium.cpu."//trim(cpart)//".vtu"
         
            write(unit=myunit,fmt='(3a)') "    <DataSet group="""" part=""0"" file=""",trim(fname),"""/>"
         
         enddo
         
         write(unit=myunit,fmt='(a)') "  </Collection>"
         write(unit=myunit,fmt='(a)') "</VTKFile>"
         
         close(myunit)

      return
!***********************************************************************************************************************************************************************************
      end subroutine write_medium_collection
!***********************************************************************************************************************************************************************************

end module mod_init_medium
