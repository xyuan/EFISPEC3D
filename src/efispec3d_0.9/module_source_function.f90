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
!!This file contains a module to compute source functions' time history.

!>@brief
!!This module contains functions and subroutines to compute tsource functions's time history. 
module mod_source_function
   
   use mpi
   
   implicit none

   private

   public  :: gabor
   public  :: expcos
   public  :: ispli3
   public  :: ricker
   public  :: spiexp
   public  :: fctanh
   public  :: fctanh_dt
   public  :: interp_linear
   private :: rvec_ascends_strictly
   private :: rvec_bracket

   contains



!
!
!>@brief function @f$N^{\circ}3@f$ to compute real part of Gabor wavelet : @f$ f(t) = a \times \exp \left[ \frac{t-ts}{2\sigma} \right]^{2} \times \cos \left[ \frac{2 \pi (t-ts)}{\lambda}  \right] @f$
!>@param t  : current time of computation @f$ t_{simu} @f$
!>@param ts : time shift of the function
!>@param l  : parameter @f$\lambda@f$. By default: @f$\sigma = \lambda/20@f$.
!>@param a  : amplitude
!***********************************************************************************************************************************************************************************
   real function gabor(t,ts,l,a)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : RG_PI

      implicit none

      real, intent(in) :: t
      real, intent(in) :: ts
      real, intent(in) :: l
      real, intent(in) :: a

      real             :: s

      s = l/20.0
      
      gabor = a*exp(-((t-ts)**2)/(2.0*s**2))* ( cos(2.0*RG_PI*(t-ts)/l) )

!***********************************************************************************************************************************************************************************
   end function gabor
!***********************************************************************************************************************************************************************************

!
!
!>@brief function @f$N^{\circ}4@f$ to compute euroseistest project source function for case can1 : @f$ f(t) = \exp \left[ - \frac{\omega(t-ts)}{\gamma} \right]^{2} \times \cos[\omega(t-ts)+\theta]  @f$
!>@param t  : current time of computation @f$ t_{simu} @f$
!>@param ts : time shift of the function
!>@param fp : pseudo-frequency of the function
!***********************************************************************************************************************************************************************************
   real function expcos(t,ts,fp)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : RG_PI

      implicit none

      real, intent(in)    :: t 
      real, intent(in)    :: ts
      real, intent(in)    :: fp

      real                :: ome
      real                :: gam
      real                :: the

      ome    = 2.0*RG_PI*fp
      gam    = 2.0
      the    = RG_PI/2.0
      expcos = exp(-(ome*(t-ts)/gam)**2)*cos(ome*(t-ts)+the)

!***********************************************************************************************************************************************************************************
   end function expcos
!***********************************************************************************************************************************************************************************

!
!
!>@brief function @f$N^{\circ}5@f$ to compute order 3 spline : see <a href="http://en.wikipedia.org/wiki/Spline_%28mathematics%29" target="_blank">Wikipedia</a>
!>@param s  : spline value
!>@param t  : current time of computation @f$ t_{simu} @f$
!>@param du : source duration
!***********************************************************************************************************************************************************************************
   subroutine ispli3(t,du,s)
!***********************************************************************************************************************************************************************************
      implicit none

      real, intent(inout) :: s  
      real, intent(in)    :: t  
      real, intent(in)    :: du 

      real                :: dtsp
      real                :: abst
      real                :: t0

      t0   = du/2.0
      dtsp = du/4.0

      abst = abs(t-t0)
      if (abst >= 2.0*dtsp) then
         s = 0.0
      elseif ( (abst < (2.0*dtsp)) .and. (abst >= dtsp) ) then
         s = ((2.0*dtsp-abst)**3)/(6.0*dtsp**3)
      else
         s = (3.0*abst**3 - 6.0*dtsp*abst**2 + 4.0*dtsp**3)/(6.0*dtsp**3)
      endif
      s = s/dtsp

      return
!***********************************************************************************************************************************************************************************
   end subroutine ispli3
!***********************************************************************************************************************************************************************************
   
!
!
!>@brief function @f$N^{\circ}6@f$ to compute order 2 Ricker wavelet : @f$ f(t) = 2 a \left[ \pi^2 \frac{(t-ts)^2}{tp^2} -0.5 \right] \times \exp  \left[ -\pi^2 \frac{(t-ts)^2}{tp^2} \right]  @f$
!>@param t  : current time of computation @f$ t_{simu} @f$
!>@param ts : time shift of the function
!>@param tp : pseudo-period of the function
!>@param a  : amplitude factor
!***********************************************************************************************************************************************************************************
   real function ricker(t,ts,tp,a)
!***********************************************************************************************************************************************************************************

      use mod_global_variables, only : RG_PI

      implicit none

      real, intent(in) :: t  
      real, intent(in) :: ts 
      real, intent(in) :: tp 
      real, intent(in) :: a 

      ricker = 2.0*a*((RG_PI**2*(t-ts)**2)/(tp**2)-0.5)*exp((-RG_PI**2*(t-ts)**2)/(tp**2)) ! max peak at -1.00
!
!***********************************************************************************************************************************************************************************
   end function ricker
!***********************************************************************************************************************************************************************************

!
!
!>@brief function @f$N^{\circ}7@f$ to compute spice project exponential source function : @f$ f(t) = a \left[1-\left(1+\frac{t}{tp}\right) \right] \times \exp \left[ -\frac{t}{tp} \right] @f$
!>@param t  : current time of computation @f$ t_{simu} @f$
!>@param tp : pseudo-period of the function
!>@param a  : amplitude factor
!***********************************************************************************************************************************************************************************
   real function spiexp(t,tp,a)
!***********************************************************************************************************************************************************************************

      implicit none

      real, intent(in) :: t  
      real, intent(in) :: tp 
      real, intent(in) :: a 

      spiexp = a*(1.0-(1.0+t/tp)*exp(-t/tp))

!***********************************************************************************************************************************************************************************
   end function spiexp
!***********************************************************************************************************************************************************************************

!
!
!>@brief function @f$N^{\circ}8@f$ to compute hyperbolic tangent function (also called 'Bouchon pulse') : @f$ f(t) = 0.5 a  \left[ 1 + \tanh \left( \frac{t-ts}{2/5 tp} \right)  \right]  @f$
!>@param t  : current time of computation @f$ t_{simu} @f$
!>@param ts : time shift of the function
!>@param tp : pseudo-period of the function
!>@param a  : amplitude factor
!***********************************************************************************************************************************************************************************
   real function fctanh(t,ts,tp,a)
!***********************************************************************************************************************************************************************************

      implicit none

      real, intent(in) :: t  
      real, intent(in) :: ts 
      real, intent(in) :: tp 
      real, intent(in) :: a 

      real             :: tau
      real             :: tt0   

      tau    = 2.0*tp
!     tt0    = 2.0*tau
      tt0    = ts
      fctanh = a*0.5*( 1.0 + tanh((t-tt0)/(tau/5.0)) )

!***********************************************************************************************************************************************************************************
   end function fctanh
!***********************************************************************************************************************************************************************************

!
!
!>@brief function @f$N^{\circ}9@f$ to compute first derivative of hyperbolic tangent function : @f$ f(t) = 2 a \times f_c \left[ 1 - \tanh^{2} \left( 4 f_c (t-ts)  \right)   \right]  @f$ with @f$f_c = \frac{1}{8/5 tp}@f$
!>@param t  : current time of computation @f$ t_{simu} @f$
!>@param ts : time shift of the function
!>@param tp : pseudo-period of the function
!>@param a  : amplitude factor
!***********************************************************************************************************************************************************************************
   real function fctanh_dt(t,ts,tp,a)
!***********************************************************************************************************************************************************************************

      implicit none

      real, intent(in) :: t 
      real, intent(in) :: tp
      real, intent(in) :: ts
      real, intent(in) :: a

      real             :: tau
      real             :: fc
      real             :: tt0   

      tau       = 2.0*tp
      fc        = 1.0/(4.0*(tau/5.0))
!     tt0       = 2.0*tau
      tt0       = ts
      fctanh_dt = +a*2.0*fc*(1.0-(tanh(4.0*fc*(t-tt0)))**2)

!***********************************************************************************************************************************************************************************
   end function fctanh_dt
!***********************************************************************************************************************************************************************************

!
!
!
!>@brief INTERP_LINEAR applies piecewise linear interpolation to data.
!!
!!From a space of DIM_NUM dimensions, we are given a sequence of
!!DATA_NUM points, which are presumed to be successive samples
!!from a curve of points P.
!!
!!We are also given a parameterization of this data, that is,
!!an associated sequence of DATA_NUM values of a variable T.
!!The values of T are assumed to be strictly increasing.
!!
!!Thus, we have a sequence of values P(T), where T is a scalar,
!!and each value of P is of dimension DIM_NUM.
!!
!!We are then given INTERP_NUM values of T, for which values P
!!are to be produced, by linear interpolation of the data we are given.
!!
!!Note that the user may request extrapolation.  This occurs whenever
!!a T_INTERP value is less than the minimum T_DATA or greater than the
!!maximum T_DATA.  In that case, linear extrapolation is used.
!! 
!!Licensing: This code is distributed under the GNU LGPL license.
!!
!! @version 03 December 2007
!!
!! @author John Burkardt
!!
!>@param dim_num                    : the spatial dimension
!>@param data_num                   : the number of data points.
!>@param t_data(data_num)           : values of the independent variable at the sample points. Values of T_DATA must be strictly increasing.
!>@param p_data(dim_num,data_num)   : values of the dependent variables at the sample points.
!>@param interp_num                 : number of points at which interpolation is to be done.
!>@param t_interp(interp_num)       : values of the independent variable at the interpolation points.
!>@param p_interp(dim_num,data_num) : interpolated values of the dependent variables at the interpolation points.
!***********************************************************************************************************************************************************************************
   subroutine interp_linear( dim_num, data_num, t_data, p_data, interp_num, t_interp, p_interp )
!***********************************************************************************************************************************************************************************

   use mod_global_variables, only : error_stop

   implicit none
 
   integer, intent( in) :: data_num
   integer, intent( in) :: dim_num
   integer, intent( in) :: interp_num
   real   , intent( in) :: t_data(data_num)
   real   , intent( in) :: t_interp(interp_num)
   real   , intent( in) :: p_data(dim_num,data_num)
   real   , intent(out) :: p_interp(dim_num,interp_num)
 
   real                 :: t
   integer              :: interp
   integer              :: left
   integer              :: right

   character(len=255)   :: info
 
   if ( .not. rvec_ascends_strictly ( data_num, t_data ) ) then
     write (info,'(a)') 'error in subroutine interp_linear: Independent variable array T_DATA is not strictly increasing.'
     call error_stop(info)
   endif
 
   do interp = 1, interp_num
 
      t = t_interp(interp)
!!
!!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!!  nearest to, TVAL.
!!
      call rvec_bracket ( data_num, t_data, t, left, right )
 
      p_interp(1:dim_num,interp) = &
       ( ( t_data(right) - t                ) * p_data(1:dim_num,left)   &
       + (                 t - t_data(left) ) * p_data(1:dim_num,right) ) &
       / ( t_data(right)     - t_data(left) )
 
   enddo

   return
!***********************************************************************************************************************************************************************************
   end subroutine interp_linear
!***********************************************************************************************************************************************************************************

!
!
!
!>@brief This function determines if an R8VEC is strictly ascending.
!
!>@license GNU LGPL
!
!>@author John Burkardt
!
!>@param n : size of the array
!>@param x : array to be examined
!>@return TRUE if the entries of x strictly ascend.
!***********************************************************************************************************************************************************************************
   logical function rvec_ascends_strictly(n,x)
!***********************************************************************************************************************************************************************************

   implicit none
 
   integer, intent(in)               :: n
   real   , intent(in), dimension(n) :: x

   integer             :: i
 
   do i = 1, n - 1
      if ( x(i+1) <= x(i) ) then
         rvec_ascends_strictly = .false.
         return
      endif
   enddo
 
   rvec_ascends_strictly = .true.
 
   return
!***********************************************************************************************************************************************************************************
   end function rvec_ascends_strictly 
!***********************************************************************************************************************************************************************************

!
!
!>@brief This subroutine searches a sorted R8VEC for successive brackets of a value.
!
!>@license GNU LGPL
!
!>@author John Burkardt
!
!>@param  n     : size of input array
!>@param  x     : an array sorted into ascending order
!>@param  xval  : a value to be bracketed
!>@param  left  : x(left) <= xval <= x(right). xval < x(1), when left = 1, right = 2. x(n) < xval, when left = n-1, right = n
!>@param  right : x(left) <= xval <= x(right). xval < x(1), when left = 1, right = 2. x(n) < xval, when left = n-1, right = n
!***********************************************************************************************************************************************************************************
   subroutine rvec_bracket(n,x,xval,left,right)
!***********************************************************************************************************************************************************************************

      implicit none
      
      integer, intent( in)               :: n 
      real   , intent( in), dimension(n) :: x
      real   , intent( in)               :: xval
      integer, intent(out)               :: left
      integer, intent(out)               :: right
      
      integer :: i

      do i = 2,n-1

         if ( xval < x(i) ) then
            left  = i - 1
            right = i
            return
         endif

      enddo
      
      left  = n - 1
      right = n
      
      return
!***********************************************************************************************************************************************************************************
   end subroutine rvec_bracket
!***********************************************************************************************************************************************************************************

end module mod_source_function
