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
!!This file contains a module to compute global @f$x,y,z@f$-coordinates of a given point in the physical domain from its local @f$\xi,\eta,\zeta@f$-coordinates

!>@brief
!!This module contains subroutines to compute global @f$x,y,z@f$-coordinates of a given point in the physical domain from its local @f$\xi,\eta,\zeta@f$-coordinates
module mod_coordinate
   
   use mpi
   
   implicit none

   private

   public  :: compute_hexa_point_coord
   public  :: compute_quad_point_coord
   public  :: compute_quad_point_coord_z

   contains

!
!
!>@brief subroutine to compute @f$x,y,z@f$-coordinates of a point knowing to which hexahedron element it belongs and its local coordinates @f$\xi,\eta,\zeta@f$.
!>@param ihexa: hexahedron element number in cpu myrank
!>@param xi   : @f$\xi  @f$ local coordinate where to compute global coordinates
!>@param eta  : @f$\eta @f$ local coordinate where to compute global coordinates
!>@param zeta : @f$\zeta@f$ local coordinate where to compute global coordinates
!>@param x    : @f$x@f$-coordinate
!>@param y    : @f$y@f$-coordinate
!>@param z    : @f$z@f$-coordinate
!***********************************************************************************************************************************************************************************
   subroutine compute_hexa_point_coord(ihexa,xi,eta,zeta,x,y,z)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only :&
                                      ig_hexa_nnode&
                                     ,ig_hexa_gnode_xiloc&
                                     ,ig_hexa_gnode_etloc&
                                     ,ig_hexa_gnode_etloc&
                                     ,ig_hexa_gnode_zeloc&
                                     ,ig_line_nnode&
                                     ,rg_gnode_x&
                                     ,rg_gnode_y&
                                     ,rg_gnode_z&
                                     ,ig_hexa_gnode_glonum

      use mod_lagrange       , only : lagrap_geom

      implicit none

      integer, intent( in) :: ihexa
      real   , intent( in) :: xi
      real   , intent( in) :: eta
      real   , intent( in) :: zeta
      real   , intent(out) :: x
      real   , intent(out) :: y
      real   , intent(out) :: z

      real                 :: lag_xi
      real                 :: lag_et
      real                 :: lag_ze
      real                 :: xgnode
      real                 :: ygnode
      real                 :: zgnode
 
      integer              :: m

      x = 0.0
      y = 0.0
      z = 0.0

      do m = 1,ig_hexa_nnode

         lag_xi = lagrap_geom(ig_hexa_gnode_xiloc(m),xi  ,ig_line_nnode)
         lag_et = lagrap_geom(ig_hexa_gnode_etloc(m),eta ,ig_line_nnode)
         lag_ze = lagrap_geom(ig_hexa_gnode_zeloc(m),zeta,ig_line_nnode)

         xgnode = rg_gnode_x(ig_hexa_gnode_glonum(m,ihexa))
         ygnode = rg_gnode_y(ig_hexa_gnode_glonum(m,ihexa))
         zgnode = rg_gnode_z(ig_hexa_gnode_glonum(m,ihexa))
                
         x      = x + lag_xi*lag_et*lag_ze*xgnode
         y      = y + lag_xi*lag_et*lag_ze*ygnode
         z      = z + lag_xi*lag_et*lag_ze*zgnode

      enddo

      return

!***********************************************************************************************************************************************************************************
   end subroutine compute_hexa_point_coord
!***********************************************************************************************************************************************************************************

!
!
!>@brief subroutine to compute @f$x,y,z@f$-coordinates of a point knowing to which quadrangle element it belongs and its local coordinates @f$\xi,\eta@f$.
!>@param iquad: quadrangle element number in cpu myrank
!>@param xi   : @f$\xi  @f$ local coordinate where to compute global coordinates
!>@param eta  : @f$\eta @f$ local coordinate where to compute global coordinates
!>@param x    : @f$x@f$-coordinate
!>@param y    : @f$y@f$-coordinate
!>@param z    : @f$z@f$-coordinate
!***********************************************************************************************************************************************************************************
   subroutine compute_quad_point_coord(iquad,xi,eta,x,y,z)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only :&
                                      ig_quad_nnode&
                                     ,ig_quad_gnode_xiloc&
                                     ,ig_quad_gnode_etloc&
                                     ,ig_line_nnode&
                                     ,rg_gnode_x&
                                     ,rg_gnode_y&
                                     ,rg_gnode_z&
                                     ,ig_quadf_gnode_glonum

      use mod_lagrange       , only : lagrap_geom

      implicit none

      integer, intent( in) :: iquad
      real   , intent( in) :: xi
      real   , intent( in) :: eta
      real   , intent(out) :: x
      real   , intent(out) :: y
      real   , intent(out) :: z

      real                 :: lag_xi
      real                 :: lag_et
      real                 :: xgnode
      real                 :: ygnode
      real                 :: zgnode
 
      integer              :: m

      x = 0.0
      y = 0.0
      z = 0.0

      do m = 1,ig_quad_nnode

         lag_xi = lagrap_geom(ig_quad_gnode_xiloc(m),xi ,ig_line_nnode)
         lag_et = lagrap_geom(ig_quad_gnode_etloc(m),eta,ig_line_nnode)

         xgnode = rg_gnode_x(ig_quadf_gnode_glonum(m,iquad))
         ygnode = rg_gnode_y(ig_quadf_gnode_glonum(m,iquad))
         zgnode = rg_gnode_z(ig_quadf_gnode_glonum(m,iquad))
                
         x      = x + lag_xi*lag_et*xgnode
         y      = y + lag_xi*lag_et*ygnode
         z      = z + lag_xi*lag_et*zgnode

      enddo

      return

!***********************************************************************************************************************************************************************************
   end subroutine compute_quad_point_coord
!***********************************************************************************************************************************************************************************

!
!
!>@brief function to compute @f$z@f$-coordinates of a point knowing to which quadrangle element it belongs and its local coordinates @f$\xi,\eta@f$.
!>@param iquad: quadrangle element number in cpu myrank
!>@param xi   : @f$\xi @f$ local coordinate where to compute global coordinates
!>@param eta  : @f$\eta@f$ local coordinate where to compute global coordinates
!***********************************************************************************************************************************************************************************
   function compute_quad_point_coord_z(iquad,xi,eta) result(z)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only :&
                                      ig_quad_nnode&
                                     ,rg_gnode_z&
                                     ,ig_quadf_gnode_glonum&
                                     ,ig_line_nnode&
                                     ,ig_quad_gnode_xiloc&
                                     ,ig_quad_gnode_etloc
      
      use mod_lagrange       , only : lagrap_geom

      implicit none
      
      integer, intent(in) :: iquad
      real   , intent(in) :: xi
      real   , intent(in) :: eta
      
      real                :: lag_xi
      real                :: lag_et
      real                :: zgnode
      real                :: z
      integer             :: m

      z = 0.0 

      do m = 1,ig_quad_nnode
   
        lag_xi = lagrap_geom(ig_quad_gnode_xiloc(m),xi ,ig_line_nnode)
        lag_et = lagrap_geom(ig_quad_gnode_etloc(m),eta,ig_line_nnode)
   
        zgnode = rg_gnode_z(ig_quadf_gnode_glonum(m,iquad))
        
        z      = z + lag_xi*lag_et*zgnode
   
      enddo

      return
!***********************************************************************************************************************************************************************************
   end function compute_quad_point_coord_z
!***********************************************************************************************************************************************************************************

end module mod_coordinate
