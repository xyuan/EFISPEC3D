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
!!This file contains a modulme to compute Newmark explicit time marching scheme, external forces @f$ F^{ext} @f$, internal forces @f$ KU @f$ and boundary traction forces @f$ C\dot{U} @f$ of the system @f$ M\ddot{U} + C\dot{U} + KU = F^{ext} @f$.

!>@brief
!!This module contains subroutines to compute Newmark explicit time marching scheme, external forces @f$ F^{ext} @f$, internal forces @f$ KU @f$ and boundary traction forces @f$ C\dot{U} @f$ of the system @f$ M\ddot{U} + C\dot{U} + KU = F^{ext} @f$.

module mod_solver

   use mpi

   implicit none

   public :: newmark_ini
   public :: newmark_end
   public :: compute_internal_forces_order4
   public :: compute_internal_forces_order5
   public :: compute_internal_forces_order6
   public :: compute_absorption_forces
   public :: compute_external_force

   contains

!
!
!>@brief This subroutine initializes Newmark time marching scheme at step n+1
!>@return displacements at step n+1     : see global variable mod_global_variables::rg_gll_displacement
!>@return external forces flush to zero : see global variable mod_global_variables::rg_gll_acceleration
!***********************************************************************************************************************************************************************************
      subroutine newmark_ini()
!***********************************************************************************************************************************************************************************

         use mpi

         use mod_global_variables, only :&
                                         ig_ngll_total&
                                        ,rg_gll_displacement&
                                        ,rg_gll_velocity&
                                        ,rg_gll_acceleration&
                                        ,rg_dt&
                                        ,rg_dt2&
                                        ,RG_NEWMARK_GAMMA
   
         implicit none
   
         integer :: igll
   
         do igll = 1,ig_ngll_total
            rg_gll_displacement(1,igll) = rg_gll_displacement(1,igll) + rg_dt*rg_gll_velocity(1,igll) + rg_dt2*RG_NEWMARK_GAMMA*rg_gll_acceleration(1,igll) !displacement x
            rg_gll_displacement(2,igll) = rg_gll_displacement(2,igll) + rg_dt*rg_gll_velocity(2,igll) + rg_dt2*RG_NEWMARK_GAMMA*rg_gll_acceleration(2,igll) !dispalcement y
            rg_gll_displacement(3,igll) = rg_gll_displacement(3,igll) + rg_dt*rg_gll_velocity(3,igll) + rg_dt2*RG_NEWMARK_GAMMA*rg_gll_acceleration(3,igll) !displacement z
         enddo
   
         do igll = 1,ig_ngll_total
            rg_gll_acceleration(1,igll) = 0.0 !flush to zero acceleration x
            rg_gll_acceleration(2,igll) = 0.0 !flush to zero acceleration y
            rg_gll_acceleration(3,igll) = 0.0 !flush to zero acceleration z
         enddo
   
         return
!***********************************************************************************************************************************************************************************
      end subroutine newmark_ini
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine finalizes Newmark time marching scheme at step n+1
!***********************************************************************************************************************************************************************************
      subroutine newmark_end()
!***********************************************************************************************************************************************************************************

         use mpi

         use mod_global_variables, only :&
                                         ig_ngll_total&
                                        ,rg_gll_velocity&
                                        ,rg_gll_acceleration&
                                        ,rg_gll_acctmp&
                                        ,rg_dt&
                                        ,RG_NEWMARK_GAMMA&
                                        ,rg_gll_mass_matrix
         implicit none
   
         integer :: igll
   
         do igll = 1,ig_ngll_total
            rg_gll_acceleration(1,igll) = rg_gll_acceleration(1,igll)*rg_gll_mass_matrix(igll) !acceleration x step n+1
            rg_gll_acceleration(2,igll) = rg_gll_acceleration(2,igll)*rg_gll_mass_matrix(igll) !acceleration y step n+1
            rg_gll_acceleration(3,igll) = rg_gll_acceleration(3,igll)*rg_gll_mass_matrix(igll) !acceleration z step n+1
         enddo
    
         do igll = 1,ig_ngll_total
            rg_gll_velocity(1,igll) = rg_gll_velocity(1,igll) + rg_dt*((1.0-RG_NEWMARK_GAMMA)*rg_gll_acctmp(1,igll) + RG_NEWMARK_GAMMA*rg_gll_acceleration(1,igll)) !velocity x step n+1
            rg_gll_velocity(2,igll) = rg_gll_velocity(2,igll) + rg_dt*((1.0-RG_NEWMARK_GAMMA)*rg_gll_acctmp(2,igll) + RG_NEWMARK_GAMMA*rg_gll_acceleration(2,igll)) !velocity y step n+1
            rg_gll_velocity(3,igll) = rg_gll_velocity(3,igll) + rg_dt*((1.0-RG_NEWMARK_GAMMA)*rg_gll_acctmp(3,igll) + RG_NEWMARK_GAMMA*rg_gll_acceleration(3,igll)) !velocity z step n+1
         enddo
    
         do igll = 1,ig_ngll_total
            rg_gll_acctmp(1,igll)   = rg_gll_acceleration(1,igll) !store tmp acceleration at step n for next step n+1
            rg_gll_acctmp(2,igll)   = rg_gll_acceleration(2,igll) !store tmp acceleration at step n for next step n+1
            rg_gll_acctmp(3,igll)   = rg_gll_acceleration(3,igll) !store tmp acceleration at step n for next step n+1
         enddo
   
         return
!***********************************************************************************************************************************************************************************
      end subroutine newmark_end
!***********************************************************************************************************************************************************************************
   
!
!
!>@brief
!!This subroutine computes internal forces @f$ \int _{\Omega}  \boldsymbol{\epsilon}(\mathbf{v}) ^{T} \colon \boldsymbol{\tau} \, d\Omega @f$ for spectral-elements of order 4.
!!Stress-strain relationship can be linear elastic (general isotropic fourth-order Hooke's law for continuous media) or viscoelastic (memory variables method).
!>@param elt_start : first hexahedron element of the loop 
!>@param elt_end   : last  hexahedron element of the loop 
!***********************************************************************************************************************************************************************************
      subroutine compute_internal_forces_order4(elt_start,elt_end)
!***********************************************************************************************************************************************************************************

         use mpi

         use mod_global_variables, only :&
                                         IG_NGLL&
                                        ,IG_NDOF&
                                        ,rg_gll_lagrange_deriv&
                                        ,rg_gll_displacement&
                                        ,rg_gll_acceleration&
                                        ,rg_gll_weight&
                                        ,rg_gll_abscissa&
                                        ,ig_nreceiver_hexa&
                                        ,ig_idt&
                                        ,ig_receiver_saving_incr&
                                        ,ig_nhexa&
                                        ,ig_hexa_gll_glonum&
                                        ,rg_hexa_gll_dxidx&
                                        ,rg_hexa_gll_dxidy&
                                        ,rg_hexa_gll_dxidz&
                                        ,rg_hexa_gll_detdx&
                                        ,rg_hexa_gll_detdy&
                                        ,rg_hexa_gll_detdz&
                                        ,rg_hexa_gll_dzedx&
                                        ,rg_hexa_gll_dzedy&
                                        ,rg_hexa_gll_dzedz&
                                        ,rg_hexa_gll_jacobian_det&
                                        ,ig_myrank&
                                        ,ig_nhexa_outer&
                                        ,rg_hexa_gll_rho&
                                        ,rg_hexa_gll_rhovs2&
                                        ,rg_hexa_gll_rhovp2&
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
                                        ,RG_RELAX_COEFF&
                                        ,rg_mem_var_exp&
                                        ,IG_NRELAX&
                                        ,rg_dt&
                                        ,LG_VISCO&
                                        ,ig_nhexa
         
         implicit none
         
         integer, intent(in) :: elt_start
         integer, intent(in) :: elt_end
         
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpx1
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpx2
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpx3
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpy1
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpy2
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpy3
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpz1
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpz2
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpz3
         real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL) :: rl_displacement_gll
         real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL) :: rl_acceleration_gll

         real :: duxdxi
         real :: duxdet
         real :: duxdze
         real :: duydxi
         real :: duydet
         real :: duydze
         real :: duzdxi
         real :: duzdet
         real :: duzdze
         real :: duxdx
         real :: duydy
         real :: duzdz
         real :: duxdy
         real :: duxdz
         real :: duydx
         real :: duydz
         real :: duzdx
         real :: duzdy
         real :: dxidx
         real :: dxidy
         real :: dxidz
         real :: detdx
         real :: detdy
         real :: detdz
         real :: dzedx
         real :: dzedy
         real :: dzedz
         real :: tauxx
         real :: tauyy
         real :: tauzz
         real :: tauxy
         real :: tauxz
         real :: tauyz
         real :: tauxx_n12
         real :: tauyy_n12
         real :: tauzz_n12
         real :: tauxy_n12
         real :: tauxz_n12
         real :: tauyz_n12
         real :: trace_tau
         real :: tmpx1
         real :: tmpx2
         real :: tmpx3
         real :: tmpx4
         real :: tmpy1
         real :: tmpy2
         real :: tmpy3
         real :: tmpz1
         real :: tmpz2
         real :: tmpz3
         real :: fac1
         real :: fac2
         real :: fac3
         
         integer :: iel
         integer :: k
         integer :: l
         integer :: m
         integer :: igll
         integer :: imem_var
         
         
         do iel = elt_start,elt_end
   
            !
            !------->flush local acceleration to zero
            do k = 1,IG_NGLL
               do l = 1,IG_NGLL
                  do m = 1,IG_NGLL 
   
                     rl_acceleration_gll(1,m,l,k) = 0.0
                     rl_acceleration_gll(2,m,l,k) = 0.0
                     rl_acceleration_gll(3,m,l,k) = 0.0
   
                  enddo
               enddo
            enddo
   
            !
            !------->fill local displacement
            do k = 1,IG_NGLL        !zeta
               do l = 1,IG_NGLL     !eta
                  do m = 1,IG_NGLL  !xi
   
                     igll                         = ig_hexa_gll_glonum(m,l,k,iel)
   
                     rl_displacement_gll(1,m,l,k) = rg_gll_displacement(1,igll)
                     rl_displacement_gll(2,m,l,k) = rg_gll_displacement(2,igll)
                     rl_displacement_gll(3,m,l,k) = rg_gll_displacement(3,igll)
   
                  enddo
               enddo
            enddo
         !
         !
         !******************************************************************************
         !->compute integrale at GLL nodes + assemble forces in global gll grid for hexa
         !******************************************************************************
            do k = 1,IG_NGLL        !zeta
               do l = 1,IG_NGLL     !eta
                  do m = 1,IG_NGLL  !xi
   
         !
         !---------->derivative of displacement with respect to local coordinate xi, eta and zeta at the gll node klm
   
                     duxdxi = rl_displacement_gll(1,1,l,k)*rg_gll_lagrange_deriv(1,m) &
                            + rl_displacement_gll(1,2,l,k)*rg_gll_lagrange_deriv(2,m) &
                            + rl_displacement_gll(1,3,l,k)*rg_gll_lagrange_deriv(3,m) &
                            + rl_displacement_gll(1,4,l,k)*rg_gll_lagrange_deriv(4,m) &
                            + rl_displacement_gll(1,5,l,k)*rg_gll_lagrange_deriv(5,m) 
   
                     duydxi = rl_displacement_gll(2,1,l,k)*rg_gll_lagrange_deriv(1,m) &
                            + rl_displacement_gll(2,2,l,k)*rg_gll_lagrange_deriv(2,m) &
                            + rl_displacement_gll(2,3,l,k)*rg_gll_lagrange_deriv(3,m) &
                            + rl_displacement_gll(2,4,l,k)*rg_gll_lagrange_deriv(4,m) &
                            + rl_displacement_gll(2,5,l,k)*rg_gll_lagrange_deriv(5,m) 
   
                     duzdxi = rl_displacement_gll(3,1,l,k)*rg_gll_lagrange_deriv(1,m) &
                            + rl_displacement_gll(3,2,l,k)*rg_gll_lagrange_deriv(2,m) &
                            + rl_displacement_gll(3,3,l,k)*rg_gll_lagrange_deriv(3,m) &
                            + rl_displacement_gll(3,4,l,k)*rg_gll_lagrange_deriv(4,m) &
                            + rl_displacement_gll(3,5,l,k)*rg_gll_lagrange_deriv(5,m)
   
                     duxdet = rl_displacement_gll(1,m,1,k)*rg_gll_lagrange_deriv(1,l) &
                            + rl_displacement_gll(1,m,2,k)*rg_gll_lagrange_deriv(2,l) &
                            + rl_displacement_gll(1,m,3,k)*rg_gll_lagrange_deriv(3,l) &
                            + rl_displacement_gll(1,m,4,k)*rg_gll_lagrange_deriv(4,l) &
                            + rl_displacement_gll(1,m,5,k)*rg_gll_lagrange_deriv(5,l)
   
                     duydet = rl_displacement_gll(2,m,1,k)*rg_gll_lagrange_deriv(1,l) &
                            + rl_displacement_gll(2,m,2,k)*rg_gll_lagrange_deriv(2,l) &
                            + rl_displacement_gll(2,m,3,k)*rg_gll_lagrange_deriv(3,l) &
                            + rl_displacement_gll(2,m,4,k)*rg_gll_lagrange_deriv(4,l) &
                            + rl_displacement_gll(2,m,5,k)*rg_gll_lagrange_deriv(5,l)
   
                     duzdet = rl_displacement_gll(3,m,1,k)*rg_gll_lagrange_deriv(1,l) &
                            + rl_displacement_gll(3,m,2,k)*rg_gll_lagrange_deriv(2,l) &
                            + rl_displacement_gll(3,m,3,k)*rg_gll_lagrange_deriv(3,l) &
                            + rl_displacement_gll(3,m,4,k)*rg_gll_lagrange_deriv(4,l) &
                            + rl_displacement_gll(3,m,5,k)*rg_gll_lagrange_deriv(5,l)
   
                     duxdze = rl_displacement_gll(1,m,l,1)*rg_gll_lagrange_deriv(1,k) &
                            + rl_displacement_gll(1,m,l,2)*rg_gll_lagrange_deriv(2,k) &
                            + rl_displacement_gll(1,m,l,3)*rg_gll_lagrange_deriv(3,k) &
                            + rl_displacement_gll(1,m,l,4)*rg_gll_lagrange_deriv(4,k) &
                            + rl_displacement_gll(1,m,l,5)*rg_gll_lagrange_deriv(5,k) 
   
                     duydze = rl_displacement_gll(2,m,l,1)*rg_gll_lagrange_deriv(1,k) &
                            + rl_displacement_gll(2,m,l,2)*rg_gll_lagrange_deriv(2,k) &
                            + rl_displacement_gll(2,m,l,3)*rg_gll_lagrange_deriv(3,k) &
                            + rl_displacement_gll(2,m,l,4)*rg_gll_lagrange_deriv(4,k) &
                            + rl_displacement_gll(2,m,l,5)*rg_gll_lagrange_deriv(5,k) 
   
                     duzdze = rl_displacement_gll(3,m,l,1)*rg_gll_lagrange_deriv(1,k) &
                            + rl_displacement_gll(3,m,l,2)*rg_gll_lagrange_deriv(2,k) &
                            + rl_displacement_gll(3,m,l,3)*rg_gll_lagrange_deriv(3,k) &
                            + rl_displacement_gll(3,m,l,4)*rg_gll_lagrange_deriv(4,k) &
                            + rl_displacement_gll(3,m,l,5)*rg_gll_lagrange_deriv(5,k)
   
         !      
         !---------->derivative of displacement at step n+1 with respect to global coordinate x, y and z at the gll node klm
                     dxidx = rg_hexa_gll_dxidx  (m,l,k,iel)
                     dxidy = rg_hexa_gll_dxidy  (m,l,k,iel)
                     dxidz = rg_hexa_gll_dxidz  (m,l,k,iel)
                     detdx = rg_hexa_gll_detdx (m,l,k,iel)
                     detdy = rg_hexa_gll_detdy (m,l,k,iel)
                     detdz = rg_hexa_gll_detdz (m,l,k,iel)
                     dzedx = rg_hexa_gll_dzedx(m,l,k,iel)
                     dzedy = rg_hexa_gll_dzedy(m,l,k,iel)
                     dzedz = rg_hexa_gll_dzedz(m,l,k,iel)
   
                     duxdx = duxdxi*dxidx + duxdet*detdx + duxdze*dzedx
                     duxdy = duxdxi*dxidy + duxdet*detdy + duxdze*dzedy
                     duxdz = duxdxi*dxidz + duxdet*detdz + duxdze*dzedz
                     duydx = duydxi*dxidx + duydet*detdx + duydze*dzedx
                     duydy = duydxi*dxidy + duydet*detdy + duydze*dzedy
                     duydz = duydxi*dxidz + duydet*detdz + duydze*dzedz
                     duzdx = duzdxi*dxidx + duzdet*detdx + duzdze*dzedx
                     duzdy = duzdxi*dxidy + duzdet*detdy + duzdze*dzedy
                     duzdz = duzdxi*dxidz + duzdet*detdz + duzdze*dzedz
         !
         !---------->compute elastic stress (elastic simulation) or unrelaxed elastic stress (viscoelastic simulation)
                     trace_tau = (rg_hexa_gll_rhovp2(m,l,k,iel) - 2.0*rg_hexa_gll_rhovs2(m,l,k,iel))*(duxdx+duydy+duzdz)
                     tauxx     = trace_tau + 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*duxdx
                     tauyy     = trace_tau + 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*duydy
                     tauzz     = trace_tau + 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*duzdz
                     tauxy     =                 rg_hexa_gll_rhovs2(m,l,k,iel)*(duxdy+duydx)
                     tauxz     =                 rg_hexa_gll_rhovs2(m,l,k,iel)*(duxdz+duzdx)
                     tauyz     =                 rg_hexa_gll_rhovs2(m,l,k,iel)*(duydz+duzdy)
         !
         !----------->compute viscoelastic stress
               
                     if(LG_VISCO) then
                   
                         do imem_var = 1,IG_NRELAX
                   
                            tmpx1 = rg_mem_var_exp(imem_var)
                   
                            tmpx2 =     rg_hexa_gll_rhovp2(m,l,k,iel)*rg_hexa_gll_wkqp(imem_var,m,l,k,iel)
                            tmpx3 = 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*rg_hexa_gll_wkqs(imem_var,m,l,k,iel)
                   
                            tmpx4 = (duxdx+duydy+duzdz)*(tmpx2 - tmpx3)
                   
         !         
         !----------------->anelastic stress at step n+1/2 following s. ma and p. liu (2006) using epsnp1 and unrelaxed material modulus
                            tauxx_n12 = tmpx1*rg_hexa_gll_ksixx(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*duxdx + tmpx4)
                            tauyy_n12 = tmpx1*rg_hexa_gll_ksiyy(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*duydy + tmpx4)
                            tauzz_n12 = tmpx1*rg_hexa_gll_ksizz(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*duzdz + tmpx4)
                            tauxy_n12 = tmpx1*rg_hexa_gll_ksixy(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*0.5*(duxdy+duydx))
                            tauxz_n12 = tmpx1*rg_hexa_gll_ksixz(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*0.5*(duxdz+duzdx))
                            tauyz_n12 = tmpx1*rg_hexa_gll_ksiyz(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*0.5*(duydz+duzdy))
         !            
         !----------------->compute final stress at step n+1 according to day and minster (1984) + ma and liu (2006) : tauxx - SUM anelastic stress at step n+1
                            tauxx = tauxx - 0.5*(tauxx_n12 + rg_hexa_gll_ksixx(imem_var,m,l,k,iel))
                            tauyy = tauyy - 0.5*(tauyy_n12 + rg_hexa_gll_ksiyy(imem_var,m,l,k,iel))
                            tauzz = tauzz - 0.5*(tauzz_n12 + rg_hexa_gll_ksizz(imem_var,m,l,k,iel))
                            tauxy = tauxy - 0.5*(tauxy_n12 + rg_hexa_gll_ksixy(imem_var,m,l,k,iel))
                            tauxz = tauxz - 0.5*(tauxz_n12 + rg_hexa_gll_ksixz(imem_var,m,l,k,iel))
                            tauyz = tauyz - 0.5*(tauyz_n12 + rg_hexa_gll_ksiyz(imem_var,m,l,k,iel))
         !            
         !----------------->update memory for stress (step n+1/2)
                            rg_hexa_gll_ksixx(imem_var,m,l,k,iel) = tauxx_n12
                            rg_hexa_gll_ksiyy(imem_var,m,l,k,iel) = tauyy_n12
                            rg_hexa_gll_ksizz(imem_var,m,l,k,iel) = tauzz_n12
                            rg_hexa_gll_ksixy(imem_var,m,l,k,iel) = tauxy_n12
                            rg_hexa_gll_ksixz(imem_var,m,l,k,iel) = tauxz_n12
                            rg_hexa_gll_ksiyz(imem_var,m,l,k,iel) = tauyz_n12
                   
                         enddo
                   
                     endif
         !      
         !---------->store members of integration of the gll node klm
                     intpx1(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxx*dxidx+tauxy*dxidy+tauxz*dxidz)
                     intpx2(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxx*detdx+tauxy*detdy+tauxz*detdz)
                     intpx3(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxx*dzedx+tauxy*dzedy+tauxz*dzedz)
   
                     intpy1(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxy*dxidx+tauyy*dxidy+tauyz*dxidz)
                     intpy2(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxy*detdx+tauyy*detdy+tauyz*detdz)
                     intpy3(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxy*dzedx+tauyy*dzedy+tauyz*dzedz)
   
                     intpz1(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxz*dxidx+tauyz*dxidy+tauzz*dxidz)
                     intpz2(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxz*detdx+tauyz*detdy+tauzz*detdz)
                     intpz3(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxz*dzedx+tauyz*dzedy+tauzz*dzedz)
   
                  enddo !xi
               enddo    !eta
            enddo       !zeta
   
         !
         !->finish integration for hexa (internal forces at step n+1)
            do k = 1,IG_NGLL
               do l = 1,IG_NGLL
                  do m = 1,IG_NGLL
   
                     tmpx1 = intpx1(1,l,k)*rg_gll_lagrange_deriv(m,1)*rg_gll_weight(1) &
                           + intpx1(2,l,k)*rg_gll_lagrange_deriv(m,2)*rg_gll_weight(2) &
                           + intpx1(3,l,k)*rg_gll_lagrange_deriv(m,3)*rg_gll_weight(3) &
                           + intpx1(4,l,k)*rg_gll_lagrange_deriv(m,4)*rg_gll_weight(4) &
                           + intpx1(5,l,k)*rg_gll_lagrange_deriv(m,5)*rg_gll_weight(5)
   
                     tmpy1 = intpy1(1,l,k)*rg_gll_lagrange_deriv(m,1)*rg_gll_weight(1) &
                           + intpy1(2,l,k)*rg_gll_lagrange_deriv(m,2)*rg_gll_weight(2) &
                           + intpy1(3,l,k)*rg_gll_lagrange_deriv(m,3)*rg_gll_weight(3) &
                           + intpy1(4,l,k)*rg_gll_lagrange_deriv(m,4)*rg_gll_weight(4) &
                           + intpy1(5,l,k)*rg_gll_lagrange_deriv(m,5)*rg_gll_weight(5)
   
                     tmpz1 = intpz1(1,l,k)*rg_gll_lagrange_deriv(m,1)*rg_gll_weight(1) &
                           + intpz1(2,l,k)*rg_gll_lagrange_deriv(m,2)*rg_gll_weight(2) &
                           + intpz1(3,l,k)*rg_gll_lagrange_deriv(m,3)*rg_gll_weight(3) &
                           + intpz1(4,l,k)*rg_gll_lagrange_deriv(m,4)*rg_gll_weight(4) &
                           + intpz1(5,l,k)*rg_gll_lagrange_deriv(m,5)*rg_gll_weight(5)
   
                     tmpx2 = intpx2(m,1,k)*rg_gll_lagrange_deriv(l,1)*rg_gll_weight(1) &
                           + intpx2(m,2,k)*rg_gll_lagrange_deriv(l,2)*rg_gll_weight(2) &
                           + intpx2(m,3,k)*rg_gll_lagrange_deriv(l,3)*rg_gll_weight(3) &
                           + intpx2(m,4,k)*rg_gll_lagrange_deriv(l,4)*rg_gll_weight(4) &
                           + intpx2(m,5,k)*rg_gll_lagrange_deriv(l,5)*rg_gll_weight(5)
   
                     tmpy2 = intpy2(m,1,k)*rg_gll_lagrange_deriv(l,1)*rg_gll_weight(1) &
                           + intpy2(m,2,k)*rg_gll_lagrange_deriv(l,2)*rg_gll_weight(2) &
                           + intpy2(m,3,k)*rg_gll_lagrange_deriv(l,3)*rg_gll_weight(3) &
                           + intpy2(m,4,k)*rg_gll_lagrange_deriv(l,4)*rg_gll_weight(4) &
                           + intpy2(m,5,k)*rg_gll_lagrange_deriv(l,5)*rg_gll_weight(5)
   
                     tmpz2 = intpz2(m,1,k)*rg_gll_lagrange_deriv(l,1)*rg_gll_weight(1) &
                           + intpz2(m,2,k)*rg_gll_lagrange_deriv(l,2)*rg_gll_weight(2) &
                           + intpz2(m,3,k)*rg_gll_lagrange_deriv(l,3)*rg_gll_weight(3) &
                           + intpz2(m,4,k)*rg_gll_lagrange_deriv(l,4)*rg_gll_weight(4) &
                           + intpz2(m,5,k)*rg_gll_lagrange_deriv(l,5)*rg_gll_weight(5)
   
                     tmpx3 = intpx3(m,l,1)*rg_gll_lagrange_deriv(k,1)*rg_gll_weight(1) &
                           + intpx3(m,l,2)*rg_gll_lagrange_deriv(k,2)*rg_gll_weight(2) &
                           + intpx3(m,l,3)*rg_gll_lagrange_deriv(k,3)*rg_gll_weight(3) &
                           + intpx3(m,l,4)*rg_gll_lagrange_deriv(k,4)*rg_gll_weight(4) &
                           + intpx3(m,l,5)*rg_gll_lagrange_deriv(k,5)*rg_gll_weight(5)
   
                     tmpy3 = intpy3(m,l,1)*rg_gll_lagrange_deriv(k,1)*rg_gll_weight(1) &
                           + intpy3(m,l,2)*rg_gll_lagrange_deriv(k,2)*rg_gll_weight(2) &
                           + intpy3(m,l,3)*rg_gll_lagrange_deriv(k,3)*rg_gll_weight(3) &
                           + intpy3(m,l,4)*rg_gll_lagrange_deriv(k,4)*rg_gll_weight(4) &
                           + intpy3(m,l,5)*rg_gll_lagrange_deriv(k,5)*rg_gll_weight(5)
   
                     tmpz3 = intpz3(m,l,1)*rg_gll_lagrange_deriv(k,1)*rg_gll_weight(1) &
                           + intpz3(m,l,2)*rg_gll_lagrange_deriv(k,2)*rg_gll_weight(2) &
                           + intpz3(m,l,3)*rg_gll_lagrange_deriv(k,3)*rg_gll_weight(3) &
                           + intpz3(m,l,4)*rg_gll_lagrange_deriv(k,4)*rg_gll_weight(4) &
                           + intpz3(m,l,5)*rg_gll_lagrange_deriv(k,5)*rg_gll_weight(5) 
   
                     fac1 = rg_gll_weight(l)*rg_gll_weight(k)
                     fac2 = rg_gll_weight(m)*rg_gll_weight(k)
                     fac3 = rg_gll_weight(m)*rg_gll_weight(l)
   
                     rl_acceleration_gll(1,m,l,k) = rl_acceleration_gll(1,m,l,k) + (fac1*tmpx1 + fac2*tmpx2 + fac3*tmpx3)
                     rl_acceleration_gll(2,m,l,k) = rl_acceleration_gll(2,m,l,k) + (fac1*tmpy1 + fac2*tmpy2 + fac3*tmpy3)
                     rl_acceleration_gll(3,m,l,k) = rl_acceleration_gll(3,m,l,k) + (fac1*tmpz1 + fac2*tmpz2 + fac3*tmpz3)
   
                  enddo
               enddo
            enddo
   
            do k = 1,IG_NGLL        !zeta
               do l = 1,IG_NGLL     !eta
                  do m = 1,IG_NGLL  !xi
   
                     igll                        = ig_hexa_gll_glonum(m,l,k,iel)
   
                     rg_gll_acceleration(1,igll) = rg_gll_acceleration(1,igll) - rl_acceleration_gll(1,m,l,k)
                     rg_gll_acceleration(2,igll) = rg_gll_acceleration(2,igll) - rl_acceleration_gll(2,m,l,k)
                     rg_gll_acceleration(3,igll) = rg_gll_acceleration(3,igll) - rl_acceleration_gll(3,m,l,k)
   
                  enddo
               enddo
            enddo
   
         enddo !loop on hexahedron elements
   
         return
!***********************************************************************************************************************************************************************************
      end subroutine compute_internal_forces_order4
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine computes internal forces @f$ \int _{\Omega}  \boldsymbol{\epsilon}(\mathbf{v}) ^{T} \colon \boldsymbol{\tau} \, d\Omega @f$ for spectral-elements of order 5.
!!Stress-strain relationship can be linear elastic (general isotropic fourth-order Hooke's law for continuous media) or viscoelastic (memory variables method).
!>@param elt_start : first hexahedron element of the loop 
!>@param elt_end   : last  hexahedron element of the loop 
!***********************************************************************************************************************************************************************************
   subroutine compute_internal_forces_order5(elt_start,elt_end)
!***********************************************************************************************************************************************************************************

         use mpi

         use mod_global_variables, only :&
                                         IG_NGLL&
                                        ,IG_NDOF&
                                        ,rg_gll_lagrange_deriv&
                                        ,rg_gll_displacement&
                                        ,rg_gll_acceleration&
                                        ,rg_gll_weight&
                                        ,rg_gll_abscissa&
                                        ,ig_nreceiver_hexa&
                                        ,ig_idt&
                                        ,ig_receiver_saving_incr&
                                        ,ig_nhexa&
                                        ,ig_hexa_gll_glonum&
                                        ,rg_hexa_gll_dxidx&
                                        ,rg_hexa_gll_dxidy&
                                        ,rg_hexa_gll_dxidz&
                                        ,rg_hexa_gll_detdx&
                                        ,rg_hexa_gll_detdy&
                                        ,rg_hexa_gll_detdz&
                                        ,rg_hexa_gll_dzedx&
                                        ,rg_hexa_gll_dzedy&
                                        ,rg_hexa_gll_dzedz&
                                        ,rg_hexa_gll_jacobian_det&
                                        ,ig_myrank&
                                        ,ig_nhexa_outer&
                                        ,rg_hexa_gll_rho&
                                        ,rg_hexa_gll_rhovs2&
                                        ,rg_hexa_gll_rhovp2&
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
                                        ,RG_RELAX_COEFF&
                                        ,rg_mem_var_exp&
                                        ,IG_NRELAX&
                                        ,rg_dt&
                                        ,LG_VISCO&
                                        ,ig_nhexa
         
         implicit none
         
         integer, intent(in) :: elt_start
         integer, intent(in) :: elt_end
         
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpx1
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpx2
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpx3
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpy1
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpy2
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpy3
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpz1
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpz2
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpz3
         real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL) :: rl_displacement_gll
         real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL) :: rl_acceleration_gll

         real :: duxdxi
         real :: duxdet
         real :: duxdze
         real :: duydxi
         real :: duydet
         real :: duydze
         real :: duzdxi
         real :: duzdet
         real :: duzdze
         real :: duxdx
         real :: duydy
         real :: duzdz
         real :: duxdy
         real :: duxdz
         real :: duydx
         real :: duydz
         real :: duzdx
         real :: duzdy
         real :: dxidx
         real :: dxidy
         real :: dxidz
         real :: detdx
         real :: detdy
         real :: detdz
         real :: dzedx
         real :: dzedy
         real :: dzedz
         real :: tauxx
         real :: tauyy
         real :: tauzz
         real :: tauxy
         real :: tauxz
         real :: tauyz
         real :: tauxx_n12
         real :: tauyy_n12
         real :: tauzz_n12
         real :: tauxy_n12
         real :: tauxz_n12
         real :: tauyz_n12
         real :: trace_tau
         real :: tmpx1
         real :: tmpx2
         real :: tmpx3
         real :: tmpx4
         real :: tmpy1
         real :: tmpy2
         real :: tmpy3
         real :: tmpz1
         real :: tmpz2
         real :: tmpz3
         real :: fac1
         real :: fac2
         real :: fac3
         
         integer :: iel
         integer :: k
         integer :: l
         integer :: m
         integer :: igll
         integer :: imem_var
         
         
         do iel = elt_start,elt_end
   
            !
            !------->flush to zero local acceleration
            do k = 1,IG_NGLL
               do l = 1,IG_NGLL
                  do m = 1,IG_NGLL 
   
                     rl_acceleration_gll(1,m,l,k) = 0.0
                     rl_acceleration_gll(2,m,l,k) = 0.0
                     rl_acceleration_gll(3,m,l,k) = 0.0
   
                  enddo
               enddo
            enddo
   
            !
            !------->fill local displacement
            do k = 1,IG_NGLL        !zeta
               do l = 1,IG_NGLL     !eta
                  do m = 1,IG_NGLL  !xi
   
                     igll                         = ig_hexa_gll_glonum(m,l,k,iel)
   
                     rl_displacement_gll(1,m,l,k) = rg_gll_displacement(1,igll)
                     rl_displacement_gll(2,m,l,k) = rg_gll_displacement(2,igll)
                     rl_displacement_gll(3,m,l,k) = rg_gll_displacement(3,igll)
   
                  enddo
               enddo
            enddo
         !
         !
         !******************************************************************************
         !->compute integrale at gll nodes + assemble force in global gll grid for hexa
         !******************************************************************************
            do k = 1,IG_NGLL        !zeta
               do l = 1,IG_NGLL     !eta
                  do m = 1,IG_NGLL  !xi
         !
         !---------->derivative of displacement with respect to local coordinate xi, eta and zeta at the gll node klm
   
   
                     duxdxi = rl_displacement_gll(1,1,l,k)*rg_gll_lagrange_deriv(1,m) &
                            + rl_displacement_gll(1,2,l,k)*rg_gll_lagrange_deriv(2,m) &
                            + rl_displacement_gll(1,3,l,k)*rg_gll_lagrange_deriv(3,m) &
                            + rl_displacement_gll(1,4,l,k)*rg_gll_lagrange_deriv(4,m) &
                            + rl_displacement_gll(1,5,l,k)*rg_gll_lagrange_deriv(5,m) &
                            + rl_displacement_gll(1,6,l,k)*rg_gll_lagrange_deriv(6,m) 
   
                     duydxi = rl_displacement_gll(2,1,l,k)*rg_gll_lagrange_deriv(1,m) &
                            + rl_displacement_gll(2,2,l,k)*rg_gll_lagrange_deriv(2,m) &
                            + rl_displacement_gll(2,3,l,k)*rg_gll_lagrange_deriv(3,m) &
                            + rl_displacement_gll(2,4,l,k)*rg_gll_lagrange_deriv(4,m) &
                            + rl_displacement_gll(2,5,l,k)*rg_gll_lagrange_deriv(5,m) &
                            + rl_displacement_gll(2,6,l,k)*rg_gll_lagrange_deriv(6,m) 
   
                     duzdxi = rl_displacement_gll(3,1,l,k)*rg_gll_lagrange_deriv(1,m) &
                            + rl_displacement_gll(3,2,l,k)*rg_gll_lagrange_deriv(2,m) &
                            + rl_displacement_gll(3,3,l,k)*rg_gll_lagrange_deriv(3,m) &
                            + rl_displacement_gll(3,4,l,k)*rg_gll_lagrange_deriv(4,m) &
                            + rl_displacement_gll(3,5,l,k)*rg_gll_lagrange_deriv(5,m) &
                            + rl_displacement_gll(3,6,l,k)*rg_gll_lagrange_deriv(6,m)
   
                     duxdet = rl_displacement_gll(1,m,1,k)*rg_gll_lagrange_deriv(1,l) &
                            + rl_displacement_gll(1,m,2,k)*rg_gll_lagrange_deriv(2,l) &
                            + rl_displacement_gll(1,m,3,k)*rg_gll_lagrange_deriv(3,l) &
                            + rl_displacement_gll(1,m,4,k)*rg_gll_lagrange_deriv(4,l) &
                            + rl_displacement_gll(1,m,5,k)*rg_gll_lagrange_deriv(5,l) &
                            + rl_displacement_gll(1,m,6,k)*rg_gll_lagrange_deriv(6,l)
   
                     duydet = rl_displacement_gll(2,m,1,k)*rg_gll_lagrange_deriv(1,l) &
                            + rl_displacement_gll(2,m,2,k)*rg_gll_lagrange_deriv(2,l) &
                            + rl_displacement_gll(2,m,3,k)*rg_gll_lagrange_deriv(3,l) &
                            + rl_displacement_gll(2,m,4,k)*rg_gll_lagrange_deriv(4,l) &
                            + rl_displacement_gll(2,m,5,k)*rg_gll_lagrange_deriv(5,l) &
                            + rl_displacement_gll(2,m,6,k)*rg_gll_lagrange_deriv(6,l)
   
                     duzdet = rl_displacement_gll(3,m,1,k)*rg_gll_lagrange_deriv(1,l) &
                            + rl_displacement_gll(3,m,2,k)*rg_gll_lagrange_deriv(2,l) &
                            + rl_displacement_gll(3,m,3,k)*rg_gll_lagrange_deriv(3,l) &
                            + rl_displacement_gll(3,m,4,k)*rg_gll_lagrange_deriv(4,l) &
                            + rl_displacement_gll(3,m,5,k)*rg_gll_lagrange_deriv(5,l) &
                            + rl_displacement_gll(3,m,6,k)*rg_gll_lagrange_deriv(6,l)
   
                     duxdze = rl_displacement_gll(1,m,l,1)*rg_gll_lagrange_deriv(1,k) &
                            + rl_displacement_gll(1,m,l,2)*rg_gll_lagrange_deriv(2,k) &
                            + rl_displacement_gll(1,m,l,3)*rg_gll_lagrange_deriv(3,k) &
                            + rl_displacement_gll(1,m,l,4)*rg_gll_lagrange_deriv(4,k) &
                            + rl_displacement_gll(1,m,l,5)*rg_gll_lagrange_deriv(5,k) &
                            + rl_displacement_gll(1,m,l,6)*rg_gll_lagrange_deriv(6,k) 
   
                     duydze = rl_displacement_gll(2,m,l,1)*rg_gll_lagrange_deriv(1,k) &
                            + rl_displacement_gll(2,m,l,2)*rg_gll_lagrange_deriv(2,k) &
                            + rl_displacement_gll(2,m,l,3)*rg_gll_lagrange_deriv(3,k) &
                            + rl_displacement_gll(2,m,l,4)*rg_gll_lagrange_deriv(4,k) &
                            + rl_displacement_gll(2,m,l,5)*rg_gll_lagrange_deriv(5,k) &
                            + rl_displacement_gll(2,m,l,6)*rg_gll_lagrange_deriv(6,k) 
   
                     duzdze = rl_displacement_gll(3,m,l,1)*rg_gll_lagrange_deriv(1,k) &
                            + rl_displacement_gll(3,m,l,2)*rg_gll_lagrange_deriv(2,k) &
                            + rl_displacement_gll(3,m,l,3)*rg_gll_lagrange_deriv(3,k) &
                            + rl_displacement_gll(3,m,l,4)*rg_gll_lagrange_deriv(4,k) &
                            + rl_displacement_gll(3,m,l,5)*rg_gll_lagrange_deriv(5,k) &
                            + rl_displacement_gll(3,m,l,6)*rg_gll_lagrange_deriv(6,k)
   
         !      
         !---------->derivative of displacement at step n+1 with respect to global coordinate x, y and z at the gll node klm
                     dxidx = rg_hexa_gll_dxidx  (m,l,k,iel)
                     dxidy = rg_hexa_gll_dxidy  (m,l,k,iel)
                     dxidz = rg_hexa_gll_dxidz  (m,l,k,iel)
                     detdx = rg_hexa_gll_detdx (m,l,k,iel)
                     detdy = rg_hexa_gll_detdy (m,l,k,iel)
                     detdz = rg_hexa_gll_detdz (m,l,k,iel)
                     dzedx = rg_hexa_gll_dzedx(m,l,k,iel)
                     dzedy = rg_hexa_gll_dzedy(m,l,k,iel)
                     dzedz = rg_hexa_gll_dzedz(m,l,k,iel)
   
                     duxdx = duxdxi*dxidx + duxdet*detdx + duxdze*dzedx
                     duxdy = duxdxi*dxidy + duxdet*detdy + duxdze*dzedy
                     duxdz = duxdxi*dxidz + duxdet*detdz + duxdze*dzedz
                     duydx = duydxi*dxidx + duydet*detdx + duydze*dzedx
                     duydy = duydxi*dxidy + duydet*detdy + duydze*dzedy
                     duydz = duydxi*dxidz + duydet*detdz + duydze*dzedz
                     duzdx = duzdxi*dxidx + duzdet*detdx + duzdze*dzedx
                     duzdy = duzdxi*dxidy + duzdet*detdy + duzdze*dzedy
                     duzdz = duzdxi*dxidz + duzdet*detdz + duzdze*dzedz
         !
         !---------->compute elastic stress (elastic simulation) or unrelaxed elastic stress (viscoelastic simulation)
                     trace_tau = (rg_hexa_gll_rhovp2(m,l,k,iel) - 2.0*rg_hexa_gll_rhovs2(m,l,k,iel))*(duxdx+duydy+duzdz)
                     tauxx     = trace_tau + 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*duxdx
                     tauyy     = trace_tau + 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*duydy
                     tauzz     = trace_tau + 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*duzdz
                     tauxy     =                 rg_hexa_gll_rhovs2(m,l,k,iel)*(duxdy+duydx)
                     tauxz     =                 rg_hexa_gll_rhovs2(m,l,k,iel)*(duxdz+duzdx)
                     tauyz     =                 rg_hexa_gll_rhovs2(m,l,k,iel)*(duydz+duzdy)
         !
         !---------->compute viscoelastic stress
   
                 if(LG_VISCO) then
   
                     do imem_var = 1,IG_NRELAX
   
                        tmpx1 = rg_mem_var_exp(imem_var)
   
                        tmpx2 =     rg_hexa_gll_rhovp2(m,l,k,iel)*rg_hexa_gll_wkqp(imem_var,m,l,k,iel)
                        tmpx3 = 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*rg_hexa_gll_wkqs(imem_var,m,l,k,iel)
   
                        tmpx4 = (duxdx+duydy+duzdz)*(tmpx2 - tmpx3)
   
         !
         !------------->anelastic stress at step n+1/2 following s. ma and p. liu (2006) using epsnp1 and unrelaxed material modulus
                        tauxx_n12 = tmpx1*rg_hexa_gll_ksixx(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*duxdx + tmpx4)
                        tauyy_n12 = tmpx1*rg_hexa_gll_ksiyy(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*duydy + tmpx4)
                        tauzz_n12 = tmpx1*rg_hexa_gll_ksizz(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*duzdz + tmpx4)
                        tauxy_n12 = tmpx1*rg_hexa_gll_ksixy(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*0.5*(duxdy+duydx))
                        tauxz_n12 = tmpx1*rg_hexa_gll_ksixz(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*0.5*(duxdz+duzdx))
                        tauyz_n12 = tmpx1*rg_hexa_gll_ksiyz(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*0.5*(duydz+duzdy))
         !        
         !------------->compute final stress at step n+1 according to day and minster (1984) + ma and liu (2006) : tauxx - SUM anelastic stress at step n+1
                        tauxx = tauxx - 0.5*(tauxx_n12 + rg_hexa_gll_ksixx(imem_var,m,l,k,iel))
                        tauyy = tauyy - 0.5*(tauyy_n12 + rg_hexa_gll_ksiyy(imem_var,m,l,k,iel))
                        tauzz = tauzz - 0.5*(tauzz_n12 + rg_hexa_gll_ksizz(imem_var,m,l,k,iel))
                        tauxy = tauxy - 0.5*(tauxy_n12 + rg_hexa_gll_ksixy(imem_var,m,l,k,iel))
                        tauxz = tauxz - 0.5*(tauxz_n12 + rg_hexa_gll_ksixz(imem_var,m,l,k,iel))
                        tauyz = tauyz - 0.5*(tauyz_n12 + rg_hexa_gll_ksiyz(imem_var,m,l,k,iel))
         !        
         !------------->update memory for stress (step n+1/2)
                        rg_hexa_gll_ksixx(imem_var,m,l,k,iel) = tauxx_n12
                        rg_hexa_gll_ksiyy(imem_var,m,l,k,iel) = tauyy_n12
                        rg_hexa_gll_ksizz(imem_var,m,l,k,iel) = tauzz_n12
                        rg_hexa_gll_ksixy(imem_var,m,l,k,iel) = tauxy_n12
                        rg_hexa_gll_ksixz(imem_var,m,l,k,iel) = tauxz_n12
                        rg_hexa_gll_ksiyz(imem_var,m,l,k,iel) = tauyz_n12
   
                     enddo
   
                 endif
         !      
         !---------->store members of integration of the gll node klm
                     intpx1(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxx*dxidx+tauxy*dxidy+tauxz*dxidz)
                     intpx2(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxx*detdx+tauxy*detdy+tauxz*detdz)
                     intpx3(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxx*dzedx+tauxy*dzedy+tauxz*dzedz)
   
                     intpy1(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxy*dxidx+tauyy*dxidy+tauyz*dxidz)
                     intpy2(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxy*detdx+tauyy*detdy+tauyz*detdz)
                     intpy3(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxy*dzedx+tauyy*dzedy+tauyz*dzedz)
   
                     intpz1(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxz*dxidx+tauyz*dxidy+tauzz*dxidz)
                     intpz2(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxz*detdx+tauyz*detdy+tauzz*detdz)
                     intpz3(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxz*dzedx+tauyz*dzedy+tauzz*dzedz)
                  enddo !xi
               enddo    !eta
            enddo       !zeta
   
         !
         !->finish integration for hexa (internal forces at step n+1)
            do k = 1,IG_NGLL
               do l = 1,IG_NGLL
                  do m = 1,IG_NGLL
   
                     tmpx1 = intpx1(1,l,k)*rg_gll_lagrange_deriv(m,1)*rg_gll_weight(1) &
                           + intpx1(2,l,k)*rg_gll_lagrange_deriv(m,2)*rg_gll_weight(2) &
                           + intpx1(3,l,k)*rg_gll_lagrange_deriv(m,3)*rg_gll_weight(3) &
                           + intpx1(4,l,k)*rg_gll_lagrange_deriv(m,4)*rg_gll_weight(4) &
                           + intpx1(5,l,k)*rg_gll_lagrange_deriv(m,5)*rg_gll_weight(5) &
                           + intpx1(6,l,k)*rg_gll_lagrange_deriv(m,6)*rg_gll_weight(6)
   
                     tmpy1 = intpy1(1,l,k)*rg_gll_lagrange_deriv(m,1)*rg_gll_weight(1) &
                           + intpy1(2,l,k)*rg_gll_lagrange_deriv(m,2)*rg_gll_weight(2) &
                           + intpy1(3,l,k)*rg_gll_lagrange_deriv(m,3)*rg_gll_weight(3) &
                           + intpy1(4,l,k)*rg_gll_lagrange_deriv(m,4)*rg_gll_weight(4) &
                           + intpy1(5,l,k)*rg_gll_lagrange_deriv(m,5)*rg_gll_weight(5) &
                           + intpy1(6,l,k)*rg_gll_lagrange_deriv(m,6)*rg_gll_weight(6) 
   
                     tmpz1 = intpz1(1,l,k)*rg_gll_lagrange_deriv(m,1)*rg_gll_weight(1) &
                           + intpz1(2,l,k)*rg_gll_lagrange_deriv(m,2)*rg_gll_weight(2) &
                           + intpz1(3,l,k)*rg_gll_lagrange_deriv(m,3)*rg_gll_weight(3) &
                           + intpz1(4,l,k)*rg_gll_lagrange_deriv(m,4)*rg_gll_weight(4) &
                           + intpz1(5,l,k)*rg_gll_lagrange_deriv(m,5)*rg_gll_weight(5) &
                           + intpz1(6,l,k)*rg_gll_lagrange_deriv(m,6)*rg_gll_weight(6)
   
                     tmpx2 = intpx2(m,1,k)*rg_gll_lagrange_deriv(l,1)*rg_gll_weight(1) &
                           + intpx2(m,2,k)*rg_gll_lagrange_deriv(l,2)*rg_gll_weight(2) &
                           + intpx2(m,3,k)*rg_gll_lagrange_deriv(l,3)*rg_gll_weight(3) &
                           + intpx2(m,4,k)*rg_gll_lagrange_deriv(l,4)*rg_gll_weight(4) &
                           + intpx2(m,5,k)*rg_gll_lagrange_deriv(l,5)*rg_gll_weight(5) &
                           + intpx2(m,6,k)*rg_gll_lagrange_deriv(l,6)*rg_gll_weight(6)
   
                     tmpy2 = intpy2(m,1,k)*rg_gll_lagrange_deriv(l,1)*rg_gll_weight(1) &
                           + intpy2(m,2,k)*rg_gll_lagrange_deriv(l,2)*rg_gll_weight(2) &
                           + intpy2(m,3,k)*rg_gll_lagrange_deriv(l,3)*rg_gll_weight(3) &
                           + intpy2(m,4,k)*rg_gll_lagrange_deriv(l,4)*rg_gll_weight(4) &
                           + intpy2(m,5,k)*rg_gll_lagrange_deriv(l,5)*rg_gll_weight(5) &
                           + intpy2(m,6,k)*rg_gll_lagrange_deriv(l,6)*rg_gll_weight(6)
   
                     tmpz2 = intpz2(m,1,k)*rg_gll_lagrange_deriv(l,1)*rg_gll_weight(1) &
                           + intpz2(m,2,k)*rg_gll_lagrange_deriv(l,2)*rg_gll_weight(2) &
                           + intpz2(m,3,k)*rg_gll_lagrange_deriv(l,3)*rg_gll_weight(3) &
                           + intpz2(m,4,k)*rg_gll_lagrange_deriv(l,4)*rg_gll_weight(4) &
                           + intpz2(m,5,k)*rg_gll_lagrange_deriv(l,5)*rg_gll_weight(5) &
                           + intpz2(m,6,k)*rg_gll_lagrange_deriv(l,6)*rg_gll_weight(6)
   
                     tmpx3 = intpx3(m,l,1)*rg_gll_lagrange_deriv(k,1)*rg_gll_weight(1) &
                           + intpx3(m,l,2)*rg_gll_lagrange_deriv(k,2)*rg_gll_weight(2) &
                           + intpx3(m,l,3)*rg_gll_lagrange_deriv(k,3)*rg_gll_weight(3) &
                           + intpx3(m,l,4)*rg_gll_lagrange_deriv(k,4)*rg_gll_weight(4) &
                           + intpx3(m,l,5)*rg_gll_lagrange_deriv(k,5)*rg_gll_weight(5) &
                           + intpx3(m,l,6)*rg_gll_lagrange_deriv(k,6)*rg_gll_weight(6)
   
                     tmpy3 = intpy3(m,l,1)*rg_gll_lagrange_deriv(k,1)*rg_gll_weight(1) &
                           + intpy3(m,l,2)*rg_gll_lagrange_deriv(k,2)*rg_gll_weight(2) &
                           + intpy3(m,l,3)*rg_gll_lagrange_deriv(k,3)*rg_gll_weight(3) &
                           + intpy3(m,l,4)*rg_gll_lagrange_deriv(k,4)*rg_gll_weight(4) &
                           + intpy3(m,l,5)*rg_gll_lagrange_deriv(k,5)*rg_gll_weight(5) &
                           + intpy3(m,l,6)*rg_gll_lagrange_deriv(k,6)*rg_gll_weight(6)
   
                     tmpz3 = intpz3(m,l,1)*rg_gll_lagrange_deriv(k,1)*rg_gll_weight(1) &
                           + intpz3(m,l,2)*rg_gll_lagrange_deriv(k,2)*rg_gll_weight(2) &
                           + intpz3(m,l,3)*rg_gll_lagrange_deriv(k,3)*rg_gll_weight(3) &
                           + intpz3(m,l,4)*rg_gll_lagrange_deriv(k,4)*rg_gll_weight(4) &
                           + intpz3(m,l,5)*rg_gll_lagrange_deriv(k,5)*rg_gll_weight(5) &
                           + intpz3(m,l,6)*rg_gll_lagrange_deriv(k,6)*rg_gll_weight(6) 
   
                     fac1 = rg_gll_weight(l)*rg_gll_weight(k)
                     fac2 = rg_gll_weight(m)*rg_gll_weight(k)
                     fac3 = rg_gll_weight(m)*rg_gll_weight(l)
   
                     rl_acceleration_gll(1,m,l,k) = rl_acceleration_gll(1,m,l,k) + (fac1*tmpx1 + fac2*tmpx2 + fac3*tmpx3)
                     rl_acceleration_gll(2,m,l,k) = rl_acceleration_gll(2,m,l,k) + (fac1*tmpy1 + fac2*tmpy2 + fac3*tmpy3)
                     rl_acceleration_gll(3,m,l,k) = rl_acceleration_gll(3,m,l,k) + (fac1*tmpz1 + fac2*tmpz2 + fac3*tmpz3)
   
                  enddo
               enddo
            enddo
   
            do k = 1,IG_NGLL        !zeta
               do l = 1,IG_NGLL     !eta
                  do m = 1,IG_NGLL  !xi
   
                     igll                        = ig_hexa_gll_glonum(m,l,k,iel)
   
                     rg_gll_acceleration(1,igll) = rg_gll_acceleration(1,igll) - rl_acceleration_gll(1,m,l,k)
                     rg_gll_acceleration(2,igll) = rg_gll_acceleration(2,igll) - rl_acceleration_gll(2,m,l,k)
                     rg_gll_acceleration(3,igll) = rg_gll_acceleration(3,igll) - rl_acceleration_gll(3,m,l,k)
   
                  enddo
               enddo
            enddo
   
         enddo !loop on hexahedron elements
   
         return
!***********************************************************************************************************************************************************************************
      end subroutine compute_internal_forces_order5
!***********************************************************************************************************************************************************************************
   
!
!
!>@brief
!!This subroutine computes internal forces @f$ \int _{\Omega}  \boldsymbol{\epsilon}(\mathbf{v}) ^{T} \colon \boldsymbol{\tau} \, d\Omega @f$ for spectral-elements of order 6.
!!Stress-strain relationship can be linear elastic (general isotropic fourth-order Hooke's law for continuous media) or viscoelastic (memory variables method).
!>@param elt_start : first hexahedron element of the loop 
!>@param elt_end   : last  hexahedron element of the loop 
!***********************************************************************************************************************************************************************************
      subroutine compute_internal_forces_order6(elt_start,elt_end)
!***********************************************************************************************************************************************************************************

         use mpi

         use mod_global_variables, only :&
                                         IG_NGLL&
                                        ,IG_NDOF&
                                        ,rg_gll_lagrange_deriv&
                                        ,rg_gll_displacement&
                                        ,rg_gll_acceleration&
                                        ,rg_gll_weight&
                                        ,rg_gll_abscissa&
                                        ,ig_nreceiver_hexa&
                                        ,ig_idt&
                                        ,ig_receiver_saving_incr&
                                        ,ig_nhexa&
                                        ,ig_hexa_gll_glonum&
                                        ,rg_hexa_gll_dxidx&
                                        ,rg_hexa_gll_dxidy&
                                        ,rg_hexa_gll_dxidz&
                                        ,rg_hexa_gll_detdx&
                                        ,rg_hexa_gll_detdy&
                                        ,rg_hexa_gll_detdz&
                                        ,rg_hexa_gll_dzedx&
                                        ,rg_hexa_gll_dzedy&
                                        ,rg_hexa_gll_dzedz&
                                        ,rg_hexa_gll_jacobian_det&
                                        ,ig_myrank&
                                        ,ig_nhexa_outer&
                                        ,rg_hexa_gll_rho&
                                        ,rg_hexa_gll_rhovs2&
                                        ,rg_hexa_gll_rhovp2&
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
                                        ,RG_RELAX_COEFF&
                                        ,rg_mem_var_exp&
                                        ,IG_NRELAX&
                                        ,rg_dt&
                                        ,LG_VISCO&
                                        ,ig_nhexa
         
         implicit none
         
         integer, intent(in) :: elt_start
         integer, intent(in) :: elt_end
         
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpx1
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpx2
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpx3
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpy1
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpy2
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpy3
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpz1
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpz2
         real, dimension(IG_NGLL,IG_NGLL,IG_NGLL)         :: intpz3
         real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL) :: rl_displacement_gll
         real, dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL) :: rl_acceleration_gll

         real :: duxdxi
         real :: duxdet
         real :: duxdze
         real :: duydxi
         real :: duydet
         real :: duydze
         real :: duzdxi
         real :: duzdet
         real :: duzdze
         real :: duxdx
         real :: duydy
         real :: duzdz
         real :: duxdy
         real :: duxdz
         real :: duydx
         real :: duydz
         real :: duzdx
         real :: duzdy
         real :: dxidx
         real :: dxidy
         real :: dxidz
         real :: detdx
         real :: detdy
         real :: detdz
         real :: dzedx
         real :: dzedy
         real :: dzedz
         real :: tauxx
         real :: tauyy
         real :: tauzz
         real :: tauxy
         real :: tauxz
         real :: tauyz
         real :: tauxx_n12
         real :: tauyy_n12
         real :: tauzz_n12
         real :: tauxy_n12
         real :: tauxz_n12
         real :: tauyz_n12
         real :: trace_tau
         real :: tmpx1
         real :: tmpx2
         real :: tmpx3
         real :: tmpx4
         real :: tmpy1
         real :: tmpy2
         real :: tmpy3
         real :: tmpz1
         real :: tmpz2
         real :: tmpz3
         real :: fac1
         real :: fac2
         real :: fac3
         
         integer :: iel
         integer :: k
         integer :: l
         integer :: m
         integer :: igll
         integer :: imem_var
         
         do iel = elt_start,elt_end
   
            !
            !------->flush to zero local acceleration
            do k = 1,IG_NGLL
               do l = 1,IG_NGLL
                  do m = 1,IG_NGLL 
   
                     rl_acceleration_gll(1,m,l,k) = 0.0
                     rl_acceleration_gll(2,m,l,k) = 0.0
                     rl_acceleration_gll(3,m,l,k) = 0.0
   
                  enddo
               enddo
            enddo
   
            !
            !------->fill local displacement
            do k = 1,IG_NGLL        !zeta
               do l = 1,IG_NGLL     !eta
                  do m = 1,IG_NGLL  !xi
   
                     igll                         = ig_hexa_gll_glonum(m,l,k,iel)
   
                     rl_displacement_gll(1,m,l,k) = rg_gll_displacement(1,igll)
                     rl_displacement_gll(2,m,l,k) = rg_gll_displacement(2,igll)
                     rl_displacement_gll(3,m,l,k) = rg_gll_displacement(3,igll)
   
                  enddo
               enddo
            enddo
         !
         !
         !******************************************************************************
         !->compute integrale at gll nodes + assemble force in global gll grid for hexa
         !******************************************************************************
            do k = 1,IG_NGLL        !zeta
               do l = 1,IG_NGLL     !eta
                  do m = 1,IG_NGLL  !xi
         !
         !---------->derivative of displacement with respect to local coordinate xi, eta and zeta at the gll node klm
   
   
                     duxdxi = rl_displacement_gll(1,1,l,k)*rg_gll_lagrange_deriv(1,m) &
                            + rl_displacement_gll(1,2,l,k)*rg_gll_lagrange_deriv(2,m) &
                            + rl_displacement_gll(1,3,l,k)*rg_gll_lagrange_deriv(3,m) &
                            + rl_displacement_gll(1,4,l,k)*rg_gll_lagrange_deriv(4,m) &
                            + rl_displacement_gll(1,5,l,k)*rg_gll_lagrange_deriv(5,m) & 
                            + rl_displacement_gll(1,6,l,k)*rg_gll_lagrange_deriv(6,m) &
                            + rl_displacement_gll(1,7,l,k)*rg_gll_lagrange_deriv(7,m) 
   
                     duydxi = rl_displacement_gll(2,1,l,k)*rg_gll_lagrange_deriv(1,m) &
                            + rl_displacement_gll(2,2,l,k)*rg_gll_lagrange_deriv(2,m) &
                            + rl_displacement_gll(2,3,l,k)*rg_gll_lagrange_deriv(3,m) &
                            + rl_displacement_gll(2,4,l,k)*rg_gll_lagrange_deriv(4,m) &
                            + rl_displacement_gll(2,5,l,k)*rg_gll_lagrange_deriv(5,m) &
                            + rl_displacement_gll(2,6,l,k)*rg_gll_lagrange_deriv(6,m) &
                            + rl_displacement_gll(2,7,l,k)*rg_gll_lagrange_deriv(7,m) 
   
                     duzdxi = rl_displacement_gll(3,1,l,k)*rg_gll_lagrange_deriv(1,m) &
                            + rl_displacement_gll(3,2,l,k)*rg_gll_lagrange_deriv(2,m) &
                            + rl_displacement_gll(3,3,l,k)*rg_gll_lagrange_deriv(3,m) &
                            + rl_displacement_gll(3,4,l,k)*rg_gll_lagrange_deriv(4,m) &
                            + rl_displacement_gll(3,5,l,k)*rg_gll_lagrange_deriv(5,m) &
                            + rl_displacement_gll(3,6,l,k)*rg_gll_lagrange_deriv(6,m) &
                            + rl_displacement_gll(3,7,l,k)*rg_gll_lagrange_deriv(7,m)
   
                     duxdet = rl_displacement_gll(1,m,1,k)*rg_gll_lagrange_deriv(1,l) &
                            + rl_displacement_gll(1,m,2,k)*rg_gll_lagrange_deriv(2,l) &
                            + rl_displacement_gll(1,m,3,k)*rg_gll_lagrange_deriv(3,l) &
                            + rl_displacement_gll(1,m,4,k)*rg_gll_lagrange_deriv(4,l) &
                            + rl_displacement_gll(1,m,5,k)*rg_gll_lagrange_deriv(5,l) &
                            + rl_displacement_gll(1,m,6,k)*rg_gll_lagrange_deriv(6,l) &
                            + rl_displacement_gll(1,m,7,k)*rg_gll_lagrange_deriv(7,l)
   
                     duydet = rl_displacement_gll(2,m,1,k)*rg_gll_lagrange_deriv(1,l) &
                            + rl_displacement_gll(2,m,2,k)*rg_gll_lagrange_deriv(2,l) &
                            + rl_displacement_gll(2,m,3,k)*rg_gll_lagrange_deriv(3,l) &
                            + rl_displacement_gll(2,m,4,k)*rg_gll_lagrange_deriv(4,l) &
                            + rl_displacement_gll(2,m,5,k)*rg_gll_lagrange_deriv(5,l) &
                            + rl_displacement_gll(2,m,6,k)*rg_gll_lagrange_deriv(6,l) &
                            + rl_displacement_gll(2,m,7,k)*rg_gll_lagrange_deriv(7,l)
   
                     duzdet = rl_displacement_gll(3,m,1,k)*rg_gll_lagrange_deriv(1,l) &
                            + rl_displacement_gll(3,m,2,k)*rg_gll_lagrange_deriv(2,l) &
                            + rl_displacement_gll(3,m,3,k)*rg_gll_lagrange_deriv(3,l) &
                            + rl_displacement_gll(3,m,4,k)*rg_gll_lagrange_deriv(4,l) &
                            + rl_displacement_gll(3,m,5,k)*rg_gll_lagrange_deriv(5,l) &
                            + rl_displacement_gll(3,m,6,k)*rg_gll_lagrange_deriv(6,l) &
                            + rl_displacement_gll(3,m,7,k)*rg_gll_lagrange_deriv(7,l)
   
                     duxdze = rl_displacement_gll(1,m,l,1)*rg_gll_lagrange_deriv(1,k) &
                            + rl_displacement_gll(1,m,l,2)*rg_gll_lagrange_deriv(2,k) &
                            + rl_displacement_gll(1,m,l,3)*rg_gll_lagrange_deriv(3,k) &
                            + rl_displacement_gll(1,m,l,4)*rg_gll_lagrange_deriv(4,k) &
                            + rl_displacement_gll(1,m,l,5)*rg_gll_lagrange_deriv(5,k) &
                            + rl_displacement_gll(1,m,l,6)*rg_gll_lagrange_deriv(6,k) &
                            + rl_displacement_gll(1,m,l,7)*rg_gll_lagrange_deriv(7,k) 
   
                     duydze = rl_displacement_gll(2,m,l,1)*rg_gll_lagrange_deriv(1,k) &
                            + rl_displacement_gll(2,m,l,2)*rg_gll_lagrange_deriv(2,k) &
                            + rl_displacement_gll(2,m,l,3)*rg_gll_lagrange_deriv(3,k) &
                            + rl_displacement_gll(2,m,l,4)*rg_gll_lagrange_deriv(4,k) &
                            + rl_displacement_gll(2,m,l,5)*rg_gll_lagrange_deriv(5,k) &
                            + rl_displacement_gll(2,m,l,6)*rg_gll_lagrange_deriv(6,k) &
                            + rl_displacement_gll(2,m,l,7)*rg_gll_lagrange_deriv(7,k) 
   
                     duzdze = rl_displacement_gll(3,m,l,1)*rg_gll_lagrange_deriv(1,k) &
                            + rl_displacement_gll(3,m,l,2)*rg_gll_lagrange_deriv(2,k) &
                            + rl_displacement_gll(3,m,l,3)*rg_gll_lagrange_deriv(3,k) &
                            + rl_displacement_gll(3,m,l,4)*rg_gll_lagrange_deriv(4,k) &
                            + rl_displacement_gll(3,m,l,5)*rg_gll_lagrange_deriv(5,k) &
                            + rl_displacement_gll(3,m,l,6)*rg_gll_lagrange_deriv(6,k) &
                            + rl_displacement_gll(3,m,l,7)*rg_gll_lagrange_deriv(7,k)
   
         !      
         !---------->derivative of displacement at step n+1 with respect to global coordinate x, y and z at the gll node klm
                     dxidx = rg_hexa_gll_dxidx  (m,l,k,iel)
                     dxidy = rg_hexa_gll_dxidy  (m,l,k,iel)
                     dxidz = rg_hexa_gll_dxidz  (m,l,k,iel)
                     detdx = rg_hexa_gll_detdx (m,l,k,iel)
                     detdy = rg_hexa_gll_detdy (m,l,k,iel)
                     detdz = rg_hexa_gll_detdz (m,l,k,iel)
                     dzedx = rg_hexa_gll_dzedx(m,l,k,iel)
                     dzedy = rg_hexa_gll_dzedy(m,l,k,iel)
                     dzedz = rg_hexa_gll_dzedz(m,l,k,iel)
   
                     duxdx = duxdxi*dxidx + duxdet*detdx + duxdze*dzedx
                     duxdy = duxdxi*dxidy + duxdet*detdy + duxdze*dzedy
                     duxdz = duxdxi*dxidz + duxdet*detdz + duxdze*dzedz
                     duydx = duydxi*dxidx + duydet*detdx + duydze*dzedx
                     duydy = duydxi*dxidy + duydet*detdy + duydze*dzedy
                     duydz = duydxi*dxidz + duydet*detdz + duydze*dzedz
                     duzdx = duzdxi*dxidx + duzdet*detdx + duzdze*dzedx
                     duzdy = duzdxi*dxidy + duzdet*detdy + duzdze*dzedy
                     duzdz = duzdxi*dxidz + duzdet*detdz + duzdze*dzedz
         !
         !---------->compute elastic stress (elastic simulation) or unrelaxed elastic stress (viscoelastic simulation)
                     trace_tau = (rg_hexa_gll_rhovp2(m,l,k,iel) - 2.0*rg_hexa_gll_rhovs2(m,l,k,iel))*(duxdx+duydy+duzdz)
                     tauxx     = trace_tau + 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*duxdx
                     tauyy     = trace_tau + 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*duydy
                     tauzz     = trace_tau + 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*duzdz
                     tauxy     =                 rg_hexa_gll_rhovs2(m,l,k,iel)*(duxdy+duydx)
                     tauxz     =                 rg_hexa_gll_rhovs2(m,l,k,iel)*(duxdz+duzdx)
                     tauyz     =                 rg_hexa_gll_rhovs2(m,l,k,iel)*(duydz+duzdy)
         !
         !---------->compute viscoelastic stress
   
                 if(LG_VISCO) then
   
                     do imem_var = 1,IG_NRELAX
   
                        tmpx1 = rg_mem_var_exp(imem_var)
   
                        tmpx2 =     rg_hexa_gll_rhovp2(m,l,k,iel)*rg_hexa_gll_wkqp(imem_var,m,l,k,iel)
                        tmpx3 = 2.0*rg_hexa_gll_rhovs2(m,l,k,iel)*rg_hexa_gll_wkqs(imem_var,m,l,k,iel)
   
                        tmpx4 = (duxdx+duydy+duzdz)*(tmpx2 - tmpx3)
   
         !
         !------------->anelastic stress at step n+1/2 following s. ma and p. liu (2006) using epsnp1 and unrelaxed material modulus
                        tauxx_n12 = tmpx1*rg_hexa_gll_ksixx(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*duxdx + tmpx4)
                        tauyy_n12 = tmpx1*rg_hexa_gll_ksiyy(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*duydy + tmpx4)
                        tauzz_n12 = tmpx1*rg_hexa_gll_ksizz(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*duzdz + tmpx4)
                        tauxy_n12 = tmpx1*rg_hexa_gll_ksixy(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*0.5*(duxdy+duydx))
                        tauxz_n12 = tmpx1*rg_hexa_gll_ksixz(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*0.5*(duxdz+duzdx))
                        tauyz_n12 = tmpx1*rg_hexa_gll_ksiyz(imem_var,m,l,k,iel) + (1.0 - tmpx1) * (tmpx3*0.5*(duydz+duzdy))
         !        
         !------------->compute final stress at step n+1 according to day and minster (1984) + ma and liu (2006) : tauxx - SUM anelastic stress at step n+1
                        tauxx = tauxx - 0.5*(tauxx_n12 + rg_hexa_gll_ksixx(imem_var,m,l,k,iel))
                        tauyy = tauyy - 0.5*(tauyy_n12 + rg_hexa_gll_ksiyy(imem_var,m,l,k,iel))
                        tauzz = tauzz - 0.5*(tauzz_n12 + rg_hexa_gll_ksizz(imem_var,m,l,k,iel))
                        tauxy = tauxy - 0.5*(tauxy_n12 + rg_hexa_gll_ksixy(imem_var,m,l,k,iel))
                        tauxz = tauxz - 0.5*(tauxz_n12 + rg_hexa_gll_ksixz(imem_var,m,l,k,iel))
                        tauyz = tauyz - 0.5*(tauyz_n12 + rg_hexa_gll_ksiyz(imem_var,m,l,k,iel))
         !        
         !------------->update memory for stress (step n+1/2)
                        rg_hexa_gll_ksixx(imem_var,m,l,k,iel) = tauxx_n12
                        rg_hexa_gll_ksiyy(imem_var,m,l,k,iel) = tauyy_n12
                        rg_hexa_gll_ksizz(imem_var,m,l,k,iel) = tauzz_n12
                        rg_hexa_gll_ksixy(imem_var,m,l,k,iel) = tauxy_n12
                        rg_hexa_gll_ksixz(imem_var,m,l,k,iel) = tauxz_n12
                        rg_hexa_gll_ksiyz(imem_var,m,l,k,iel) = tauyz_n12
   
                     enddo
   
                 endif
         !      
         !---------->store members of integration of the gll node klm
                     intpx1(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxx*dxidx+tauxy*dxidy+tauxz*dxidz)
                     intpx2(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxx*detdx+tauxy*detdy+tauxz*detdz)
                     intpx3(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxx*dzedx+tauxy*dzedy+tauxz*dzedz)
   
                     intpy1(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxy*dxidx+tauyy*dxidy+tauyz*dxidz)
                     intpy2(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxy*detdx+tauyy*detdy+tauyz*detdz)
                     intpy3(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxy*dzedx+tauyy*dzedy+tauyz*dzedz)
   
                     intpz1(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxz*dxidx+tauyz*dxidy+tauzz*dxidz)
                     intpz2(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxz*detdx+tauyz*detdy+tauzz*detdz)
                     intpz3(m,l,k) = rg_hexa_gll_jacobian_det(m,l,k,iel)*(tauxz*dzedx+tauyz*dzedy+tauzz*dzedz)
                  enddo !xi
               enddo    !eta
            enddo       !zeta
   
         !
         !->finish integration for hexa (internal forces at step n+1)
            do k = 1,IG_NGLL
               do l = 1,IG_NGLL
                  do m = 1,IG_NGLL
   
                     tmpx1 = intpx1(1,l,k)*rg_gll_lagrange_deriv(m,1)*rg_gll_weight(1) &
                           + intpx1(2,l,k)*rg_gll_lagrange_deriv(m,2)*rg_gll_weight(2) &
                           + intpx1(3,l,k)*rg_gll_lagrange_deriv(m,3)*rg_gll_weight(3) &
                           + intpx1(4,l,k)*rg_gll_lagrange_deriv(m,4)*rg_gll_weight(4) &
                           + intpx1(5,l,k)*rg_gll_lagrange_deriv(m,5)*rg_gll_weight(5) &
                           + intpx1(6,l,k)*rg_gll_lagrange_deriv(m,6)*rg_gll_weight(6) &
                           + intpx1(7,l,k)*rg_gll_lagrange_deriv(m,7)*rg_gll_weight(7)
   
                     tmpy1 = intpy1(1,l,k)*rg_gll_lagrange_deriv(m,1)*rg_gll_weight(1) &
                           + intpy1(2,l,k)*rg_gll_lagrange_deriv(m,2)*rg_gll_weight(2) &
                           + intpy1(3,l,k)*rg_gll_lagrange_deriv(m,3)*rg_gll_weight(3) &
                           + intpy1(4,l,k)*rg_gll_lagrange_deriv(m,4)*rg_gll_weight(4) &
                           + intpy1(5,l,k)*rg_gll_lagrange_deriv(m,5)*rg_gll_weight(5) &
                           + intpy1(6,l,k)*rg_gll_lagrange_deriv(m,6)*rg_gll_weight(6) &
                           + intpy1(7,l,k)*rg_gll_lagrange_deriv(m,7)*rg_gll_weight(7)
   
                     tmpz1 = intpz1(1,l,k)*rg_gll_lagrange_deriv(m,1)*rg_gll_weight(1) &
                           + intpz1(2,l,k)*rg_gll_lagrange_deriv(m,2)*rg_gll_weight(2) &
                           + intpz1(3,l,k)*rg_gll_lagrange_deriv(m,3)*rg_gll_weight(3) &
                           + intpz1(4,l,k)*rg_gll_lagrange_deriv(m,4)*rg_gll_weight(4) &
                           + intpz1(5,l,k)*rg_gll_lagrange_deriv(m,5)*rg_gll_weight(5) &
                           + intpz1(6,l,k)*rg_gll_lagrange_deriv(m,6)*rg_gll_weight(6) &
                           + intpz1(7,l,k)*rg_gll_lagrange_deriv(m,7)*rg_gll_weight(7)
   
                     tmpx2 = intpx2(m,1,k)*rg_gll_lagrange_deriv(l,1)*rg_gll_weight(1) &
                           + intpx2(m,2,k)*rg_gll_lagrange_deriv(l,2)*rg_gll_weight(2) &
                           + intpx2(m,3,k)*rg_gll_lagrange_deriv(l,3)*rg_gll_weight(3) &
                           + intpx2(m,4,k)*rg_gll_lagrange_deriv(l,4)*rg_gll_weight(4) &
                           + intpx2(m,5,k)*rg_gll_lagrange_deriv(l,5)*rg_gll_weight(5) &
                           + intpx2(m,6,k)*rg_gll_lagrange_deriv(l,6)*rg_gll_weight(6) &
                           + intpx2(m,7,k)*rg_gll_lagrange_deriv(l,7)*rg_gll_weight(7)
   
                     tmpy2 = intpy2(m,1,k)*rg_gll_lagrange_deriv(l,1)*rg_gll_weight(1) &
                           + intpy2(m,2,k)*rg_gll_lagrange_deriv(l,2)*rg_gll_weight(2) &
                           + intpy2(m,3,k)*rg_gll_lagrange_deriv(l,3)*rg_gll_weight(3) &
                           + intpy2(m,4,k)*rg_gll_lagrange_deriv(l,4)*rg_gll_weight(4) &
                           + intpy2(m,5,k)*rg_gll_lagrange_deriv(l,5)*rg_gll_weight(5) &
                           + intpy2(m,6,k)*rg_gll_lagrange_deriv(l,6)*rg_gll_weight(6) &
                           + intpy2(m,7,k)*rg_gll_lagrange_deriv(l,7)*rg_gll_weight(7)
   
                     tmpz2 = intpz2(m,1,k)*rg_gll_lagrange_deriv(l,1)*rg_gll_weight(1) &
                           + intpz2(m,2,k)*rg_gll_lagrange_deriv(l,2)*rg_gll_weight(2) &
                           + intpz2(m,3,k)*rg_gll_lagrange_deriv(l,3)*rg_gll_weight(3) &
                           + intpz2(m,4,k)*rg_gll_lagrange_deriv(l,4)*rg_gll_weight(4) &
                           + intpz2(m,5,k)*rg_gll_lagrange_deriv(l,5)*rg_gll_weight(5) &
                           + intpz2(m,6,k)*rg_gll_lagrange_deriv(l,6)*rg_gll_weight(6) &
                           + intpz2(m,7,k)*rg_gll_lagrange_deriv(l,7)*rg_gll_weight(7)
   
                     tmpx3 = intpx3(m,l,1)*rg_gll_lagrange_deriv(k,1)*rg_gll_weight(1) &
                           + intpx3(m,l,2)*rg_gll_lagrange_deriv(k,2)*rg_gll_weight(2) &
                           + intpx3(m,l,3)*rg_gll_lagrange_deriv(k,3)*rg_gll_weight(3) &
                           + intpx3(m,l,4)*rg_gll_lagrange_deriv(k,4)*rg_gll_weight(4) &
                           + intpx3(m,l,5)*rg_gll_lagrange_deriv(k,5)*rg_gll_weight(5) &
                           + intpx3(m,l,6)*rg_gll_lagrange_deriv(k,6)*rg_gll_weight(6) &
                           + intpx3(m,l,7)*rg_gll_lagrange_deriv(k,7)*rg_gll_weight(7)
   
                     tmpy3 = intpy3(m,l,1)*rg_gll_lagrange_deriv(k,1)*rg_gll_weight(1) &
                           + intpy3(m,l,2)*rg_gll_lagrange_deriv(k,2)*rg_gll_weight(2) &
                           + intpy3(m,l,3)*rg_gll_lagrange_deriv(k,3)*rg_gll_weight(3) &
                           + intpy3(m,l,4)*rg_gll_lagrange_deriv(k,4)*rg_gll_weight(4) &
                           + intpy3(m,l,5)*rg_gll_lagrange_deriv(k,5)*rg_gll_weight(5) &
                           + intpy3(m,l,6)*rg_gll_lagrange_deriv(k,6)*rg_gll_weight(6) &
                           + intpy3(m,l,7)*rg_gll_lagrange_deriv(k,7)*rg_gll_weight(7)
   
                     tmpz3 = intpz3(m,l,1)*rg_gll_lagrange_deriv(k,1)*rg_gll_weight(1) &
                           + intpz3(m,l,2)*rg_gll_lagrange_deriv(k,2)*rg_gll_weight(2) &
                           + intpz3(m,l,3)*rg_gll_lagrange_deriv(k,3)*rg_gll_weight(3) &
                           + intpz3(m,l,4)*rg_gll_lagrange_deriv(k,4)*rg_gll_weight(4) &
                           + intpz3(m,l,5)*rg_gll_lagrange_deriv(k,5)*rg_gll_weight(5) &
                           + intpz3(m,l,6)*rg_gll_lagrange_deriv(k,6)*rg_gll_weight(6) &
                           + intpz3(m,l,7)*rg_gll_lagrange_deriv(k,7)*rg_gll_weight(7) 
   
                     fac1 = rg_gll_weight(l)*rg_gll_weight(k)
                     fac2 = rg_gll_weight(m)*rg_gll_weight(k)
                     fac3 = rg_gll_weight(m)*rg_gll_weight(l)
   
                     rl_acceleration_gll(1,m,l,k) = rl_acceleration_gll(1,m,l,k) + (fac1*tmpx1 + fac2*tmpx2 + fac3*tmpx3)
                     rl_acceleration_gll(2,m,l,k) = rl_acceleration_gll(2,m,l,k) + (fac1*tmpy1 + fac2*tmpy2 + fac3*tmpy3)
                     rl_acceleration_gll(3,m,l,k) = rl_acceleration_gll(3,m,l,k) + (fac1*tmpz1 + fac2*tmpz2 + fac3*tmpz3)
   
                  enddo
               enddo
            enddo
   
            do k = 1,IG_NGLL        !zeta
               do l = 1,IG_NGLL     !eta
                  do m = 1,IG_NGLL  !xi
   
                     igll                        = ig_hexa_gll_glonum(m,l,k,iel)
   
                     rg_gll_acceleration(1,igll) = rg_gll_acceleration(1,igll) - rl_acceleration_gll(1,m,l,k)
                     rg_gll_acceleration(2,igll) = rg_gll_acceleration(2,igll) - rl_acceleration_gll(2,m,l,k)
                     rg_gll_acceleration(3,igll) = rg_gll_acceleration(3,igll) - rl_acceleration_gll(3,m,l,k)
   
                  enddo
               enddo
            enddo
   
         enddo !loop on hexahedron elements
   
         return
!***********************************************************************************************************************************************************************************
      end subroutine compute_internal_forces_order6
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine computes absorption forces @f$ \int _{\Gamma}  \mathbf{v} ^{T} \cdot \mathbf{T} \, d\Gamma @f$ for any spectral-elements order.
!!A so-called 'P1' explicit paraxial formulation is used to approximate the traction.
!***********************************************************************************************************************************************************************************
   subroutine compute_absorption_forces()
!***********************************************************************************************************************************************************************************

      use mpi

      use mod_global_variables, only : &
                                       rg_gll_velocity&
                                      ,rg_gll_acceleration&
                                      ,rg_gll_acctmp&
                                      ,rg_dt&
                                      ,RG_NEWMARK_GAMMA&
                                      ,IG_NGLL&
                                      ,rg_gll_weight&
                                      ,ig_nquad_parax&
                                      ,rg_quadp_gll_jaco_det&
                                      ,rg_quadp_gll_normal&
                                      ,ig_quadp_gll_glonum&
                                      ,rg_quadp_gll_rhovp&
                                      ,rg_quadp_gll_rhovs

      implicit none

      real          :: vx
      real          :: vy
      real          :: vz
      real          :: jaco
      real          :: nx
      real          :: ny
      real          :: nz
      real          :: tx
      real          :: ty
      real          :: tz
      real          :: vn

      integer       :: iquad
      integer       :: k
      integer       :: l
      integer       :: igll

      do iquad = 1,ig_nquad_parax

         do k = 1,IG_NGLL
            do l = 1,IG_NGLL
    
               jaco   =  rg_quadp_gll_jaco_det(l,k,iquad)
               nx     =  rg_quadp_gll_normal(1,l,k,iquad)
               ny     =  rg_quadp_gll_normal(2,l,k,iquad)
               nz     =  rg_quadp_gll_normal(3,l,k,iquad)
   
               igll   = ig_quadp_gll_glonum(l,k,iquad)
      
               vx     = rg_gll_velocity(1,igll) + rg_dt*(RG_NEWMARK_GAMMA*rg_gll_acctmp(1,igll)) 
               vy     = rg_gll_velocity(2,igll) + rg_dt*(RG_NEWMARK_GAMMA*rg_gll_acctmp(2,igll)) 
               vz     = rg_gll_velocity(3,igll) + rg_dt*(RG_NEWMARK_GAMMA*rg_gll_acctmp(3,igll)) 
   
               vn     = vx*nx+vy*ny+vz*nz
    
               tx     =  rg_quadp_gll_rhovp(l,k,iquad)*vn*nx + rg_quadp_gll_rhovs(l,k,iquad)*(vx-vn*nx)
               ty     =  rg_quadp_gll_rhovp(l,k,iquad)*vn*ny + rg_quadp_gll_rhovs(l,k,iquad)*(vy-vn*ny)
               tz     =  rg_quadp_gll_rhovp(l,k,iquad)*vn*nz + rg_quadp_gll_rhovs(l,k,iquad)*(vz-vn*nz)
       
               rg_gll_acceleration(1,igll) = rg_gll_acceleration(1,igll) - jaco*rg_gll_weight(l)*rg_gll_weight(k)*tx
               rg_gll_acceleration(2,igll) = rg_gll_acceleration(2,igll) - jaco*rg_gll_weight(l)*rg_gll_weight(k)*ty
               rg_gll_acceleration(3,igll) = rg_gll_acceleration(3,igll) - jaco*rg_gll_weight(l)*rg_gll_weight(k)*tz
   
            enddo
         enddo

      enddo
   
      return
!***********************************************************************************************************************************************************************************
      end subroutine compute_absorption_forces
!***********************************************************************************************************************************************************************************

!>@brief
!!This subroutine sets external forces @f$ F^{ext} @f$ of the system @f$ M\ddot{U} + C\dot{U} + KU = F^{ext} @f$ for double couple and single force point sources.
!>@return : filled external force vector : see global variable mod_global_variables::rg_gll_acceleration
subroutine compute_external_force()

   use mpi

   use mod_global_variables, only :& 
                                   IG_NGLL&
                                  ,rg_simu_current_time&
                                  ,rg_gll_acceleration&
                                  ,ig_ndcsource&
                                  ,ig_nsfsource&
                                  ,tg_dcsource&
                                  ,tg_sfsource&
                                  ,rg_dcsource_gll_force&
                                  ,ig_hexa_gll_glonum&
                                  ,ig_idt&
                                  ,rg_dt&
                                  ,rg_dcsource_user_func&
                                  ,rg_sfsource_user_func
   
   use mod_source_function
   
   implicit none
   
   real       :: s
   real       :: fac
   real       :: val
   real, save :: val_dc
   real, save :: val_pf 
   
   integer    :: iso
   integer    :: ipf
   integer    :: k
   integer    :: l
   integer    :: m
   integer    :: igll
   integer    :: idir


   !
   !
   !*****************************************************************************************************
   !double couple point source located at any location
   !*****************************************************************************************************
   do iso = 1,ig_ndcsource

      if (tg_dcsource(iso)%icur == 3) then

         val = gabor(rg_simu_current_time,tg_sfsource(ipf)%shift_time,tg_sfsource(ipf)%rise_time,1.0)
   
      elseif (tg_dcsource(iso)%icur == 4) then

         val = expcos(rg_simu_current_time,tg_dcsource(iso)%shift_time,tg_dcsource(iso)%rise_time)

      elseif (tg_dcsource(iso)%icur == 5) then

         if (ig_idt == 1) val_dc = 0.0
         call ispli3(rg_simu_current_time,tg_dcsource(iso)%rise_time,s)
         val_dc = val_dc + s*rg_dt
         val    = val_dc

      elseif (tg_dcsource(iso)%icur == 6) then

         val = ricker(rg_simu_current_time,tg_dcsource(iso)%shift_time,tg_dcsource(iso)%rise_time,1.0)

      elseif (tg_dcsource(iso)%icur == 7) then

         val = spiexp(rg_simu_current_time,tg_dcsource(iso)%rise_time,1.0)

      elseif (tg_dcsource(iso)%icur == 8) then

         val = fctanh(rg_simu_current_time,tg_dcsource(iso)%shift_time,tg_dcsource(iso)%rise_time,1.0)

      elseif (tg_dcsource(iso)%icur == 9) then

         val = fctanh_dt(rg_simu_current_time,tg_dcsource(iso)%shift_time,tg_dcsource(iso)%rise_time,1.0)

      elseif (tg_dcsource(iso)%icur == 10) then

         val = rg_dcsource_user_func(ig_idt)

      endif
   
      do k = 1,IG_NGLL
         do l = 1,IG_NGLL
            do m = 1,IG_NGLL
               igll                        = ig_hexa_gll_glonum(m,l,k,tg_dcsource(iso)%iel)
               rg_gll_acceleration(1,igll) = rg_gll_acceleration(1,igll) + rg_dcsource_gll_force(1,m,l,k,iso)*val
               rg_gll_acceleration(2,igll) = rg_gll_acceleration(2,igll) + rg_dcsource_gll_force(2,m,l,k,iso)*val
               rg_gll_acceleration(3,igll) = rg_gll_acceleration(3,igll) + rg_dcsource_gll_force(3,m,l,k,iso)*val
            enddo
         enddo
      enddo
   
   enddo


   !
   !
   !*****************************************************************************************************
   !single force point source located at gll nodes
   !*****************************************************************************************************   
   do ipf = 1,ig_nsfsource
   
      if (tg_sfsource(ipf)%icur == 3) then

         val = gabor(rg_simu_current_time,tg_sfsource(ipf)%shift_time,tg_sfsource(ipf)%rise_time,1.0)

      elseif (tg_sfsource(ipf)%icur == 4) then

         val = expcos(rg_simu_current_time,tg_sfsource(ipf)%shift_time,tg_sfsource(ipf)%rise_time)

      elseif (tg_sfsource(ipf)%icur == 5) then

         if (ig_idt == 1) val_pf = 0.0
         call ispli3(rg_simu_current_time,tg_sfsource(ipf)%rise_time,s)
         val_pf = val_pf + s*rg_dt
         val    = val_pf

      elseif (tg_sfsource(ipf)%icur == 6) then

         val = ricker(rg_simu_current_time,tg_sfsource(ipf)%shift_time,tg_sfsource(ipf)%rise_time,1.0)

      elseif (tg_sfsource(ipf)%icur == 7) then

         val = spiexp(rg_simu_current_time,tg_sfsource(ipf)%rise_time,1.0)

      elseif (tg_sfsource(ipf)%icur == 8) then

         val = fctanh(rg_simu_current_time,tg_sfsource(ipf)%shift_time,tg_sfsource(ipf)%rise_time,1.0)

      elseif (tg_sfsource(ipf)%icur == 9) then

         val = fctanh_dt(rg_simu_current_time,tg_sfsource(ipf)%shift_time,tg_sfsource(ipf)%rise_time,1.0)

      elseif (tg_sfsource(ipf)%icur == 10) then

         val = rg_sfsource_user_func(ig_idt)

      endif
   
      igll                           = tg_sfsource(ipf)%iequ
      idir                           = tg_sfsource(ipf)%idir
      fac                            = tg_sfsource(ipf)%fac 
      rg_gll_acceleration(idir,igll) = rg_gll_acceleration(idir,igll) + val*fac
   
   enddo
   
   return

end subroutine compute_external_force

end module mod_solver
