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
!!This file contains a module to compute jacobian matrix.

!>@brief
!!This module contains subroutines to compute jacobian matrix. 
module mod_jacobian
   
   use mpi
   
   implicit none

   private

   public  :: compute_hexa_jacobian
   public  :: compute_quad_jacobian

   contains
   
!
!
!>@brief
!!This subroutine computes jacobian matrix and its determinant at location @f$\xi,\eta,\zeta@f$ in hexahedron element ihexa
!>@param ihexa   : hexahedron element number in cpu myrank
!>@param xisol   : @f$\xi  @f$ local coordinate where to compute jacobian matrix and its determinant
!>@param etsol   : @f$\eta @f$ local coordinate where to compute jacobian matrix and its determinant
!>@param zesol   : @f$\zeta@f$ local coordinate where to compute jacobian matrix and its determinant
!>@param dxidx   : @f$\frac{\partial \xi  }{\partial x}@f$
!>@param dxidy   : @f$\frac{\partial \xi  }{\partial y}@f$
!>@param dxidz   : @f$\frac{\partial \xi  }{\partial z}@f$
!>@param detdx   : @f$\frac{\partial \eta }{\partial x}@f$
!>@param detdy   : @f$\frac{\partial \eta }{\partial y}@f$
!>@param detdz   : @f$\frac{\partial \eta }{\partial z}@f$
!>@param dzedx   : @f$\frac{\partial \zeta}{\partial x}@f$
!>@param dzedy   : @f$\frac{\partial \zeta}{\partial y}@f$
!>@param dzedz   : @f$\frac{\partial \zeta}{\partial z}@f$
!>@param det     : determinant of jacobian matrix
!***********************************************************************************************************************************************************************************
   subroutine compute_hexa_jacobian(ihexa,xisol,etsol,zesol,dxidx,dxidy,dxidz,detdx,detdy,detdz,dzedx,dzedy,dzedz,det)
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only : ig_hexa_gnode_xiloc&
                                      ,ig_hexa_gnode_etloc&
                                      ,ig_hexa_gnode_zeloc&
                                      ,ig_line_nnode&
                                      ,rg_gnode_x&
                                      ,rg_gnode_y&
                                      ,rg_gnode_z&
                                      ,ig_hexa_nnode&
                                      ,ig_hexa_gnode_glonum&
                                      ,ig_myrank

      use mod_lagrange        , only :&
                                      lagrap_geom&
                                     ,lagrad_geom
      
      implicit none
      
      real   , intent( in) :: xisol
      real   , intent( in) :: etsol
      real   , intent( in) :: zesol
      integer, intent( in) :: ihexa
      real   , intent(out) :: dxidx
      real   , intent(out) :: dxidy
      real   , intent(out) :: dxidz
      real   , intent(out) :: detdx
      real   , intent(out) :: detdy
      real   , intent(out) :: detdz
      real   , intent(out) :: dzedx
      real   , intent(out) :: dzedy
      real   , intent(out) :: dzedz
      real   , intent(out) :: det
 
      real   , parameter   :: EPS = 1.0e-8
      real                 :: inv_det
      real                 :: xxi
      real                 :: xet
      real                 :: xze
      real                 :: yxi
      real                 :: yet
      real                 :: yze
      real                 :: zxi
      real                 :: zet
      real                 :: zze
      real                 :: lag_deriv_xi
      real                 :: lag_deriv_et
      real                 :: lag_deriv_ze
      real                 :: lag_xi
      real                 :: lag_et
      real                 :: lag_ze
      real                 :: xgnode
      real                 :: ygnode
      real                 :: zgnode
      integer              :: m
      integer              :: ios

!
!
!---->initialize partial derivatives
      xxi = 0.0
      xet = 0.0
      xze = 0.0
      yxi = 0.0
      yet = 0.0
      yze = 0.0
      zxi = 0.0
      zet = 0.0
      zze = 0.0
   
!
!
!---->compute partial derivatives
      do m = 1,ig_hexa_nnode
   
         lag_deriv_xi = lagrad_geom(ig_hexa_gnode_xiloc(m),xisol,ig_line_nnode)
         lag_deriv_et = lagrad_geom(ig_hexa_gnode_etloc(m),etsol,ig_line_nnode)
         lag_deriv_ze = lagrad_geom(ig_hexa_gnode_zeloc(m),zesol,ig_line_nnode)
         
         lag_xi       = lagrap_geom(ig_hexa_gnode_xiloc(m),xisol,ig_line_nnode)
         lag_et       = lagrap_geom(ig_hexa_gnode_etloc(m),etsol,ig_line_nnode)
         lag_ze       = lagrap_geom(ig_hexa_gnode_zeloc(m),zesol,ig_line_nnode)
         
         xgnode       = rg_gnode_x(ig_hexa_gnode_glonum(m,ihexa))
         ygnode       = rg_gnode_y(ig_hexa_gnode_glonum(m,ihexa))
         zgnode       = rg_gnode_z(ig_hexa_gnode_glonum(m,ihexa))
   
         xxi          = xxi + lag_deriv_xi*lag_et*lag_ze*xgnode
         yxi          = yxi + lag_deriv_xi*lag_et*lag_ze*ygnode
         zxi          = zxi + lag_deriv_xi*lag_et*lag_ze*zgnode
                      
         xet          = xet + lag_xi*lag_deriv_et*lag_ze*xgnode
         yet          = yet + lag_xi*lag_deriv_et*lag_ze*ygnode
         zet          = zet + lag_xi*lag_deriv_et*lag_ze*zgnode
                      
         xze          = xze + lag_xi*lag_et*lag_deriv_ze*xgnode
         yze          = yze + lag_xi*lag_et*lag_deriv_ze*ygnode
         zze          = zze + lag_xi*lag_et*lag_deriv_ze*zgnode
   
      enddo
!
!
!---->compute determinant and its inverse
      det     = xxi*yet*zze - xxi*yze*zet + xet*yze*zxi - xet*yxi*zze + xze*yxi*zet - xze*yet*zxi
      inv_det = 1.0/det
   
      if (det.lt.EPS) then
         write(*,*) "***error in compute_hexa_jacobian: det null or negative, det = ",det
         write(*,*) "***element ",ihexa, "cpu ",ig_myrank
         call mpi_abort(mpi_comm_world,100,ios)
         stop
      endif

!
!
!---->compute partial derivatives
      dxidx = (yet*zze-yze*zet)*inv_det
      dxidy = (xze*zet-xet*zze)*inv_det
      dxidz = (xet*yze-xze*yet)*inv_det
      detdx = (yze*zxi-yxi*zze)*inv_det
      detdy = (xxi*zze-xze*zxi)*inv_det
      detdz = (xze*yxi-xxi*yze)*inv_det
      dzedx = (yxi*zet-yet*zxi)*inv_det
      dzedy = (xet*zxi-xxi*zet)*inv_det
      dzedz = (xxi*yet-xet*yxi)*inv_det

      return
!***********************************************************************************************************************************************************************************
   end subroutine compute_hexa_jacobian
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine computes jacobian matrix and its determinant at location @f$\xi,\eta@f$ in quadrangle element iquad
!>@param iquad   : quadrangle element number in cpu myrank
!>@param xisol   : @f$\xi  @f$ local coordinate where to compute jacobian matrix and its determinant
!>@param etsol   : @f$\eta @f$ local coordinate where to compute jacobian matrix and its determinant
!>@param dxidx   : @f$\frac{\partial \xi  }{\partial x}@f$
!>@param dxidy   : @f$\frac{\partial \xi  }{\partial y}@f$
!>@param detdx   : @f$\frac{\partial \eta }{\partial x}@f$
!>@param detdy   : @f$\frac{\partial \eta }{\partial y}@f$
!>@param det     : determinant of jacobian matrix
!***********************************************************************************************************************************************************************************
   subroutine compute_quad_jacobian(iquad,xisol,etsol,dxidx,dxidy,detdx,detdy,det)
!***********************************************************************************************************************************************************************************
   
      use mpi
   
      use mod_global_variables, only : ig_quad_gnode_xiloc&
                                      ,ig_quad_gnode_etloc&
                                      ,ig_line_nnode&
                                      ,rg_gnode_x&
                                      ,rg_gnode_y&
                                      ,rg_gnode_z&
                                      ,ig_quad_nnode&
                                      ,ig_quadf_gnode_glonum&
                                      ,ig_myrank
      
      use mod_lagrange        , only :&
                                      lagrap_geom&
                                     ,lagrad_geom
      implicit none
      
      real   , intent(out) :: dxidx
      real   , intent(out) :: dxidy
      real   , intent(out) :: detdx
      real   , intent(out) :: detdy
      real   , intent( in) :: xisol
      real   , intent( in) :: etsol
      integer, intent( in) :: iquad
 
      real , parameter     :: EPS = 1.0e-8
      real                 :: xxi
      real                 :: xet
      real                 :: yxi
      real                 :: yet
      real                 :: det
      real                 :: inv_det
      real                 :: lag_deriv_xi
      real                 :: lag_deriv_et
      real                 :: lag_xi
      real                 :: lag_et
      real                 :: xgnode
      real                 :: ygnode
      integer              :: m
      integer              :: ios

!
!
!---->initialize partial derivatives
      xxi = 0.0
      xet = 0.0
      yxi = 0.0
      yet = 0.0
   
!
!
!---->compute partial derivatives
      do m = 1,ig_quad_nnode
   
         lag_deriv_xi = lagrad_geom(ig_quad_gnode_xiloc(m),xisol,ig_line_nnode)
         lag_deriv_et = lagrad_geom(ig_quad_gnode_etloc(m),etsol,ig_line_nnode)
         
         lag_xi       = lagrap_geom(ig_quad_gnode_xiloc(m),xisol,ig_line_nnode)
         lag_et       = lagrap_geom(ig_quad_gnode_etloc(m),etsol,ig_line_nnode)
         
         xgnode       = rg_gnode_x(ig_quadf_gnode_glonum(m,iquad))
         ygnode       = rg_gnode_y(ig_quadf_gnode_glonum(m,iquad))
   
         xxi          = xxi + lag_deriv_xi*lag_et*xgnode
         yxi          = yxi + lag_deriv_xi*lag_et*ygnode
                      
         xet          = xet + lag_xi*lag_deriv_et*xgnode
         yet          = yet + lag_xi*lag_deriv_et*ygnode
                       
      enddo

!
!
!---->compute determinant and its inverse
      det     = xxi*yet - xet*yxi
      inv_det = 1.0/det
   
      if (det.lt.EPS) then
         write(*,*) "***error in compute_quad_jacobian: det null or negative, det = ",det
         write(*,*) "***element ",iquad, "cpu ",ig_myrank
         call mpi_abort(mpi_comm_world,100,ios)
         stop
      endif

!
!
!---->compute partial derivatives
      dxidx = (+yet)*inv_det
      dxidy = (-xet)*inv_det
      detdx = (-yxi)*inv_det
      detdy = (+xxi)*inv_det
   
      return
!***********************************************************************************************************************************************************************************
   end subroutine compute_quad_jacobian
!***********************************************************************************************************************************************************************************

end module mod_jacobian
