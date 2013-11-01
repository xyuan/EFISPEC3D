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
!!This file contains a module to compute Lagrange polynomials and its derivatives.

!>@brief
!!This module contains functions to compute Lagrange polynomials and its derivatives; and to interpolate @f$x,y,z@f$-variables at @f$\xi,\eta,\zeta@f$-coordinates. 
module mod_lagrange
   
   use mpi
   
   implicit none

   private

   public :: lagrap
   public :: lagrap_geom
   public :: lagrad
   public :: lagrad_geom
   public :: hexa_lagrange_interp
   public :: quad_lagrange_interp

   contains
   
!  
!
!>@brief function to compute value of order \link mod_global_variables::ig_lagrange_order @f$n@f$ \endlink Lagrange polynomial of the GLL node i at abscissa @f$x@f$: @f$ h_i(x) = \displaystyle\prod_{l=1,\, l\neq i}^{N_{GLL}} \frac{x-x_l}{x_i-x_l} @f$
!>@param x : abscissa where Lagrange polynomial is computed
!>@param i : local position (i.e. node) in the reference domain [-1:1] where lagrange polynomial = 1 (0 at others nodes)
!>@param n : number of GLL nodes in the reference domain [-1:1] (see \link mod_global_variables::ig_ngll \endlink)
!***********************************************************************************************************************************************************************************
   real function lagrap(i,x,n)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only :&
                                      rg_gll_abscissa&
                                     ,rg_gll_abscissa_dist

      implicit none

      real   , intent(in) :: x 
      integer, intent(in) :: i 
      integer, intent(in) :: n 

      integer             :: j

      lagrap = 1.0
      do j = 1,n
         if (j.ne.i) lagrap = lagrap*( (x-rg_gll_abscissa(j))*rg_gll_abscissa_dist(j,i) )
      enddo

!***********************************************************************************************************************************************************************************
   end function lagrap
!***********************************************************************************************************************************************************************************

!
!
!>@brief function to compute value of order @f$n_g@f$ Lagrange polynomial of geometric node i at abscissa @f$x@f$: @f$ h_i(x) = \displaystyle\prod_{l=1,\, l\neq i}^{n_{g}} \frac{x-x_l}{x_i-x_l} @f$
!>@param x : abscissa where Lagrange polynomial is computed
!>@param i : local position (i.e. node) in the reference domain [-1:1] where lagrange polynomial = 1 (0 at others nodes)
!>@param n : number of geometric nodes in the reference domain [-1:1]
!***********************************************************************************************************************************************************************************
   real function lagrap_geom(i,x,n)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only :&
                                      rg_gnode_abscissa_dist&
                                     ,rg_gnode_abscissa

      implicit none

      real   , intent(in) :: x 
      integer, intent(in) :: i 
      integer, intent(in) :: n 

      integer             :: j

      lagrap_geom = 1.0
      do j = 1,n
         if (j.ne.i) lagrap_geom = lagrap_geom*( (x-rg_gnode_abscissa(j))*rg_gnode_abscissa_dist(j,i) )
      enddo

!***********************************************************************************************************************************************************************************
   end function lagrap_geom
!***********************************************************************************************************************************************************************************

!
!
!>@brief function to compute value of the derivative of order \link mod_global_variables::ig_lagrange_order @f$n@f$ \endlink Lagrange polynomial of GLL node i at abscissa @f$x@f$: 
!!@f$ h'_{i}(x) = \displaystyle\sum^{N_{GLL}}_{j \ne i} \frac{1}{x_i - x_j} \displaystyle\prod^{N_{GLL}}_{k \ne i, k \ne j} \frac{x-x_k}{x_i - x_k} @f$
!>@param x : abscissa where the derivative of Lagrange polynomial is computed
!>@param i : local position (i.e. node) in the reference domain [-1:1] where lagrange polynomial = 1 (0 at others nodes). Note: lagrange polynomial = 1 at i location but not its derivative.
!>@param n : number of GLL nodes in the reference domain [-1:1]
!***********************************************************************************************************************************************************************************
   real function lagrad(i,x,n)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : rg_gll_abscissa&
                                      ,rg_gll_abscissa_dist

      implicit none

      real   , intent(in) :: x   
      integer, intent(in) :: i   
      integer, intent(in) :: n   

      real                :: ltmp
      real                :: lprod
      integer             :: j
      integer             :: k

      lagrad = 0.0
      lprod  = 0.0
      do j = 1,n
         if (j.ne.i) then
            lprod = 1.0
            do k = 1,n
               if (k.ne.i) then
                  if(k.eq.j) then
                     ltmp = rg_gll_abscissa_dist(j,i)
                  else
                     ltmp = (x-rg_gll_abscissa(k))*rg_gll_abscissa_dist(k,i)
                  endif
                  lprod = lprod*ltmp
               endif
            enddo
         endif
         lagrad = lagrad + lprod
         lprod  = 0.0
      enddo
!
!***********************************************************************************************************************************************************************************
   end function lagrad
!***********************************************************************************************************************************************************************************

!
!
!>@brief function to compute value of the derivative of order @f$n_g@f$ Lagrange polynomial of geometric node i at abscissa @f$x@f$: 
!!@f$ h'_{i}(x) = \displaystyle\sum^{n_{g}}_{j \ne i} \frac{1}{x_i - x_j} \displaystyle\prod^{n_{g}}_{k \ne i, k \ne j} \frac{x-x_k}{x_i - x_k} @f$
!>@param x : abscissa where the derivative of Lagrange polynomial is computed
!>@param i : local position (i.e. node) in the reference domain [-1:1] where lagrange polynomial = 1 (0 at others nodes). Note: lagrange polynomial = 1 at i location but not its derivative.
!>@param n : number of geometric nodes in the reference domain [-1:1]
!***********************************************************************************************************************************************************************************
   real function lagrad_geom(i,x,n)
!***********************************************************************************************************************************************************************************
     
      use mod_global_variables, only :&
                                      rg_gnode_abscissa&
                                     ,rg_gnode_abscissa_dist

      implicit none

      real   , intent(in) :: x   
      integer, intent(in) :: i   
      integer, intent(in) :: n   

      real                :: ltmp
      real                :: lprod
      integer             :: j
      integer             :: k

      lagrad_geom = 0.0
      lprod  = 0.0
      do j = 1,n
         if (j.ne.i) then
            lprod = 1.0
            do k = 1,n
               if (k.ne.i) then
                  if(k.eq.j) then
                     ltmp = rg_gnode_abscissa_dist(j,i)
                  else
                     ltmp =  (x-rg_gnode_abscissa(k))*rg_gnode_abscissa_dist(k,i)
                  endif
                  lprod = lprod*ltmp
               endif
            enddo
         endif
         lagrad_geom = lagrad_geom + lprod
         lprod  = 0.0
      enddo
!
!***********************************************************************************************************************************************************************************
   end function lagrad_geom
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine interpolates GLL nodes @f$x,y,z@f$-values of a hexahedron element using Lagrange polynomials.
!>@param gll_values : GLL nodes @f$x,y,z@f$-values of a hexahedron element
!>@param lag        : pre-computed values of Lagrange polynomials at location @f$\xi,\eta,\zeta@f$ in a hexahedron element
!>@param x          : interpolated @f$x@f$-values at location @f$\xi,\eta,\zeta@f$
!>@param y          : interpolated @f$y@f$-values at location @f$\xi,\eta,\zeta@f$
!>@param z          : interpolated @f$z@f$-values at location @f$\xi,\eta,\zeta@f$
!***********************************************************************************************************************************************************************************
   subroutine hexa_lagrange_interp(gll_values,lag,x,y,z)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only :&
                                      IG_NGLL&
                                     ,IG_NDOF

      implicit none

      real, intent( in), dimension(IG_NDOF,IG_NGLL,IG_NGLL,IG_NGLL) :: gll_values
      real, intent( in), dimension(        IG_NGLL,IG_NGLL,IG_NGLL) :: lag
      real, intent(out)                                             :: x
      real, intent(out)                                             :: y
      real, intent(out)                                             :: z
                                                                       
      integer                                                       :: k
      integer                                                       :: l
      integer                                                       :: m

 
      x = 0.0
      y = 0.0
      z = 0.0

      do k = 1,IG_NGLL
         do l = 1,IG_NGLL
            do m = 1,IG_NGLL

               x = x + gll_values(1,m,l,k)*lag(m,l,k)
               y = y + gll_values(2,m,l,k)*lag(m,l,k)
               z = z + gll_values(3,m,l,k)*lag(m,l,k)
 
            enddo
         enddo
      enddo

      return
!***********************************************************************************************************************************************************************************
   end subroutine hexa_lagrange_interp
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!This subroutine interpolates GLL nodes @f$x,y,z@f$-values of a quadrangle element using Lagrange polynomials.
!>@param gll_values : GLL nodes @f$x,y,z@f$-values of a quadrangle element
!>@param lag        : pre-computed values of Lagrange polynomials at a location @f$\xi,\eta@f$ in a quadrangle element
!>@param x          : @f$x@f$-value interpolation
!>@param y          : @f$y@f$-value interpolation
!>@param z          : @f$z@f$-value interpolation
!***********************************************************************************************************************************************************************************
   subroutine quad_lagrange_interp(gll_values,lag,x,y,z)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only :&
                                      IG_NGLL&
                                     ,IG_NDOF

      implicit none

      real, intent( in), dimension(IG_NDOF,IG_NGLL,IG_NGLL) :: gll_values
      real, intent( in), dimension(        IG_NGLL,IG_NGLL) :: lag
      real, intent(out)                                     :: x
      real, intent(out)                                     :: y
      real, intent(out)                                     :: z
                                                               
      integer                                               :: k
      integer                                               :: l

 
      x   = 0.0
      y   = 0.0
      z   = 0.0

      do k = 1,IG_NGLL
         do l = 1,IG_NGLL

            x   = x + gll_values(1,l,k)*lag(l,k)
            y   = y + gll_values(2,l,k)*lag(l,k)
            z   = z + gll_values(3,l,k)*lag(l,k)

         enddo
      enddo

      return
!***********************************************************************************************************************************************************************************
   end subroutine quad_lagrange_interp
!***********************************************************************************************************************************************************************************

end module mod_lagrange
