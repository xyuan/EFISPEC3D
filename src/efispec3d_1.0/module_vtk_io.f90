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
!!This file written by S. Zaghi (<a href="https://github.com/szaghi/Lib_VTK_IO" target="_blank">https://github.com/szaghi/Lib_VTK_IO</a>) 
!!contains a module to write ascii or binary VTK files.
!>@author S. Zaghi (see original sources at <a href="https://github.com/szaghi/Lib_VTK_IO" target="_blank">https://github.com/szaghi/Lib_VTK_IO</a>). 
!!Original sources have been modified under GNU GPL V3 license terms to suit EFISPEC3D purpose.

!>@brief
!!This module programmed by S. Zaghi (<a href="https://github.com/szaghi/Lib_VTK_IO" target="_blank">https://github.com/szaghi/Lib_VTK_IO</a>) 
!!contains functions to write ascii or binary VTK files.
module mod_vtk_io

   use mod_global_variables, only : get_newunit

   implicit none

   private

   !functions for VTK XML
   public :: VTK_INI_XML
   public :: VTK_GEO_XML
   public :: VTK_CON_XML
   public :: VTK_DAT_XML
   public :: VTK_VAR_XML
   public :: VTK_END_XML

   !portable kind-precision
   public :: R4P, FR4P
   public :: R_P, FR_P
   public :: I4P, FI4P
   public :: I1P, FI1P
   public :: I_P, FI_P

   !----------------------------------------------------------------------------------------------------------------------------------
   
!>@brief
!!overloading of VTK_GEO_XML
   interface VTK_GEO_XML
     module procedure VTK_GEO_XML_STRG_R4, & ! real(R4P) StructuredGrid
                      VTK_GEO_XML_RECT_R4, & ! real(R4P) RectilinearGrid
                      VTK_GEO_XML_UNST_R4, & ! real(R4P) UnstructuredGrid
                      VTK_GEO_XML_CLOSEP     ! closing tag "Piece" function
   endinterface

!>@brief
!!overloading of VTK_VAR_XML
   interface VTK_VAR_XML
     module procedure VTK_VAR_XML_SCAL_R4, & ! real   (R4P) scalar
                      VTK_VAR_XML_SCAL_I4, & ! integer(I4P) scalar
                      VTK_VAR_XML_VECT_R4, & ! real   (R4P) vectorial
                      VTK_VAR_XML_VECT_I4    ! integer(I4P) vectorial
   endinterface

  
!******************************** 
!Real precision definitions
!******************************** 

!>6 digits, range @f$[\pm 10^{-37},\pm 10^{+37}-1]@f$
   integer      , parameter :: R4P = selected_real_kind(6,37)

!>default real precision
   integer      , parameter :: R_P = R4P

   
!******************************** 
!Integer precision definitions
!******************************** 

!>range @f$[-2^{31},+2^{31}-1]@f$
   integer      , parameter :: I4P = selected_int_kind(9)

!>range @f$[-2^{7},+2^{7}-1]@f$
   integer      , parameter :: I1P = selected_int_kind(2)

!>default integer precision
   integer      , parameter :: I_P = I4P

   
!******************************** 
!Real output formats
!******************************** 

!>R4P   output format
   character( 9), parameter :: FR4P  = '(E14.6E2)' 

!>R\_P  output format
   character(10), parameter :: FR_P  = '(E23.15E3)'
   

!******************************** 
!Integer output formats
!******************************** 

!>I4P  output format
   character(5) , parameter :: FI4P  = '(I12)'

!>I1P  output format
   character(4) , parameter :: FI1P  = '(I5)' 

!>I\_P output format
   character(5) , parameter :: FI_P  = '(I12)'
  

!******************************** 
!private variables
!******************************** 

!>max number of characters os static string
   integer(I4P) , parameter :: maxlen = 500         

!>end-character for binary-record finalize
   character(1) , parameter :: end_rec = char(10)    

!>ascii-output-format parameter identifier
   integer(I4P) , parameter :: f_out_ascii = 0           

!>binary-output-format parameter identifier
   integer(I4P) , parameter :: f_out_binary = 1           

!>current output-format (initialized to ascii format)
   integer(I4P) :: f_out = f_out_ascii 

!>internal logical unit
   integer(I4P) :: Unit_VTK                   

!>internal logical unit for raw binary XML append file
   integer(I4P) :: Unit_VTK_Append            

!>number of byte to be written/read
   integer(I4P) :: N_Byte                     

!>offset pointer
   integer(I4P) :: ioffset                    

!>indent pointer
   integer(I4P) :: indent                     

!>prototype of I4P integer
   integer(I4P) :: Tipo_I4                    

!>prototype of I1P integer
   integer(I1P) :: Tipo_I1                    

!>prototype of R4P integer
   real(R4P) :: Tipo_R4                    

!>mesh topology
   character(len=maxlen) :: topology                   

!----------------------------------------------------------------------------------------------------------------------------------

   contains

!
!
!>@brief
!!The VTK\_INI\_XML function is used for initializing file. This function must be the first to be called.
!>@param output_format : output format: ASCII or BINARY
!>@param filename      : file name 
!>@param mesh_topology : VTK mesh topology (e.g., RectilinearGrid, StructuredGrid, etc.)
!>@param nx1           : initial node of x axis
!>@param nx2           : final   node of x axis
!>@param ny1           : initial node of y axis
!>@param ny2           : final   node of y axis
!>@param nz1           : initial node of z axis
!>@param nz2           : final   node of z axis
!***********************************************************************************************************************************************************************************
   function VTK_INI_XML(output_format,filename,mesh_topology,nx1,nx2,ny1,ny2,nz1,nz2) result(E_IO)
!***********************************************************************************************************************************************************************************

      implicit none
    
      !--------------------------------------------------------------------------------------------------------------------------------
      character(*), intent(IN)           :: output_format ! output format: ASCII or BINARY
      character(*), intent(IN)           :: filename      ! file name
      character(*), intent(IN)           :: mesh_topology ! mesh topology
      integer(I4P), intent(IN), optional :: nx1,nx2       ! initial and final nodes of x axis
      integer(I4P), intent(IN), optional :: ny1,ny2       ! initial and final nodes of y axis
      integer(I4P), intent(IN), optional :: nz1,nz2       ! initial and final nodes of z axis
    
      integer(I4P)                       :: E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      character(len=maxlen)              :: s_buffer      ! buffer string
      !--------------------------------------------------------------------------------------------------------------------------------
    
      topology = trim(mesh_topology)
    
      Unit_VTK = get_newunit()
    
      select case(trim(output_format))
    
         case('ASCII')
           f_out = f_out_ascii
           open(unit=Unit_VTK,file=trim(filename),form='FORMATTED',access='SEQUENTIAL',action='WRITE',iostat=E_IO)
           ! writing header of file
           write(unit=Unit_VTK,fmt='(A)'                ,iostat=E_IO)'<?xml version="1.0"?>'
           write(unit=Unit_VTK,fmt='(A)'                ,iostat=E_IO)'<VTKFile type="'//trim(topology)//'" version="0.1" byte_order="BigEndian">'
           indent = 2
           select case(trim(topology))
           case('RectilinearGrid','StructuredGrid')
             write(unit=Unit_VTK,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<'//trim(topology)//' WholeExtent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
           case('UnstructuredGrid')
             write(unit=Unit_VTK,fmt='(A)'              ,iostat=E_IO)repeat(' ',indent)//'<'//trim(topology)//'>'
           endselect
           indent = indent + 2
         
         case('BINARY')
           f_out = f_out_binary
           open(unit=Unit_VTK,file=trim(filename),form='UNFORMATTED',access='SEQUENTIAL',action='WRITE',convert='BIG_ENDIAN',buffered='YES',recordtype='STREAM',iostat=E_IO)
           ! writing header of file
           write(unit=Unit_VTK,                     iostat=E_IO)'<?xml version="1.0"?>'//end_rec
           write(unit=Unit_VTK,                     iostat=E_IO)'<VTKFile type="'//trim(topology)//'" version="0.1" byte_order="BigEndian">'//end_rec
           indent = 2
           select case(trim(topology))
           case('RectilinearGrid','StructuredGrid')
             write(s_buffer,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<'//trim(topology)//' WholeExtent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
           case('UnstructuredGrid')
             write(s_buffer,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<'//trim(topology)//'>'
           endselect
           write(unit=Unit_VTK,                     iostat=E_IO)trim(s_buffer)//end_rec
           indent = indent + 2
    
           !opening the SCRATCH file used for appending raw binary data
           Unit_VTK_Append=get_newunit()
           open(unit=Unit_VTK_Append,form='UNFORMATTED',access='SEQUENTIAL',action='WRITE',convert='BIG_ENDIAN',recordtype='STREAM',buffered='YES',status='SCRATCH',iostat=E_IO)
           ioffset = 0 ! initializing offset puntator
    
      endselect
    
      return

!***********************************************************************************************************************************************************************************
   endfunction VTK_INI_XML
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!The VTK\_GEO\_XML\_STRG\_R4 function is used for saving mesh with a VTK StructuredGrid R4P topology.
!>@param nx1           : initial node of x axis
!>@param nx2           : final   node of x axis
!>@param ny1           : initial node of y axis
!>@param ny2           : final   node of y axis
!>@param nz1           : initial node of z axis
!>@param nz2           : final   node of z axis
!>@param NN            : total number of nodes
!>@param X             : @f$x@f$-coordinates
!>@param Y             : @f$y@f$-coordinates
!>@param Z             : @f$z@f$-coordinates
!***********************************************************************************************************************************************************************************
   function VTK_GEO_XML_STRG_R4(nx1,nx2,ny1,ny2,nz1,nz2,NN,X,Y,Z) result(E_IO)
!***********************************************************************************************************************************************************************************

      implicit none
    
      !--------------------------------------------------------------------------------------------------------------------------------
      integer(I4P), intent(IN) :: nx1,nx2  ! initial and final nodes of x axis
      integer(I4P), intent(IN) :: ny1,ny2  ! initial and final nodes of y axis
      integer(I4P), intent(IN) :: nz1,nz2  ! initial and final nodes of z axis
      integer(I4P), intent(IN) :: NN       ! number of all nodes
      real   (R4P), intent(IN) :: X(1:NN)  ! x coordinates
      real   (R4P), intent(IN) :: Y(1:NN)  ! y coordinates
      real   (R4P), intent(IN) :: Z(1:NN)  ! z coordinates
    
      integer(I4P)             :: E_IO     ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      integer(I4P)             :: n1       ! counter
    
      character(len=maxlen)    :: s_buffer ! buffer string
      !--------------------------------------------------------------------------------------------------------------------------------

    
      select case(f_out)
    
         case(f_out_ascii)
           write(unit=Unit_VTK,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
           indent = indent + 2
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<Points>'
           indent = indent + 2
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" NumberOfComponents="3" Name="Point" format="ascii">'
           write(unit=Unit_VTK,fmt='(3'//FR4P//')',    iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</DataArray>'
           indent = indent - 2
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</Points>'
    
         case(f_out_binary)
           write(s_buffer,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
           indent = indent + 2
           write(unit=Unit_VTK,                   iostat=E_IO)trim(s_buffer)//end_rec
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'<Points>'//end_rec
           indent = indent + 2
           write(s_buffer,fmt='(I8)',             iostat=E_IO)ioffset
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" NumberOfComponents="3" Name="Point" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = 3*NN*sizeof(Tipo_R4)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,            iostat=E_IO)N_Byte,'R4',3*NN
           write(unit=Unit_VTK_Append,            iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
           indent = indent - 2
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</Points>'//end_rec
    
      endselect
    
      return
!***********************************************************************************************************************************************************************************
   endfunction VTK_GEO_XML_STRG_R4
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!The VTK\_GEO\_XML\_RECT\_R4 function is used for saving mesh with a VTK RectilinearGrid R4P topology.
!>@param nx1           : initial node of x axis
!>@param nx2           : final   node of x axis
!>@param ny1           : initial node of y axis
!>@param ny2           : final   node of y axis
!>@param nz1           : initial node of z axis
!>@param nz2           : final   node of z axis
!>@param X             : @f$x@f$-coordinates
!>@param Y             : @f$y@f$-coordinates
!>@param Z             : @f$z@f$-coordinates
!***********************************************************************************************************************************************************************************
   function VTK_GEO_XML_RECT_R4(nx1,nx2,ny1,ny2,nz1,nz2,X,Y,Z) result(E_IO)
!***********************************************************************************************************************************************************************************

      implicit none
    
      !--------------------------------------------------------------------------------------------------------------------------------
      integer(I4P), intent(IN) :: nx1,nx2    ! initial and final nodes of x axis
      integer(I4P), intent(IN) :: ny1,ny2    ! initial and final nodes of y axis
      integer(I4P), intent(IN) :: nz1,nz2    ! initial and final nodes of z axis
      real   (R4P), intent(IN) :: X(nx1:nx2) ! x coordinates
      real   (R4P), intent(IN) :: Y(ny1:ny2) ! y coordinates
      real   (R4P), intent(IN) :: Z(nz1:nz2) ! z coordinates
    
      integer(I4P)             :: E_IO       ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      integer(I4P)             :: n1         ! counter
    
      character(len=maxlen)    :: s_buffer   ! buffer string
      !--------------------------------------------------------------------------------------------------------------------------------
    
      select case(f_out)
    
         case(f_out_ascii)
           write(unit=Unit_VTK,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
           indent = indent + 2
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<Coordinates>'
           indent = indent + 2
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="X" format="ascii">'
           write(unit=Unit_VTK,fmt=FR4P,               iostat=E_IO)(X(n1),n1=nx1,nx2)
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</DataArray>'
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="Y" format="ascii">'
           write(unit=Unit_VTK,fmt=FR4P,               iostat=E_IO)(Y(n1),n1=ny1,ny2)
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</DataArray>'
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="Z" format="ascii">'
           write(unit=Unit_VTK,fmt=FR4P,               iostat=E_IO)(Z(n1),n1=nz1,nz2)
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</DataArray>'
           indent = indent - 2
           write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</Coordinates>'
    
         case(f_out_binary)
           write(s_buffer,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
           indent = indent + 2
           write(unit=Unit_VTK,                   iostat=E_IO)trim(s_buffer)//end_rec
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'<Coordinates>'//end_rec
           indent = indent + 2
           write(s_buffer,fmt='(I8)',             iostat=E_IO)ioffset
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="X" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = (nx2-nx1+1)*sizeof(Tipo_R4)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,            iostat=E_IO)N_Byte,'R4',nx2-nx1+1
           write(unit=Unit_VTK_Append,            iostat=E_IO)(X(n1),n1=nx1,nx2)
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
           write(s_buffer,fmt='(I8)',             iostat=E_IO)ioffset
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="Y" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = (ny2-ny1+1)*sizeof(Tipo_R4)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,            iostat=E_IO)N_Byte,'R4',ny2-ny1+1
           write(unit=Unit_VTK_Append,            iostat=E_IO)(Y(n1),n1=ny1,ny2)
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
           write(s_buffer,fmt='(I8)',             iostat=E_IO)ioffset
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="Z" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = (nz2-nz1+1)*sizeof(Tipo_R4)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,            iostat=E_IO)N_Byte,'R4',nz2-nz1+1
           write(unit=Unit_VTK_Append,            iostat=E_IO)(Z(n1),n1=nz1,nz2)
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
           indent = indent - 2
           write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</Coordinates>'//end_rec
    
      endselect
    
      return
!***********************************************************************************************************************************************************************************
   endfunction VTK_GEO_XML_RECT_R4
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!The VTK\_GEO\_XML\_UNST\_R4 function is used for saving mesh with a VTK UnstructuredGrid R4P topology.
!>@param NN            : total number of nodes
!>@param NC            : total number of cells
!>@param X             : @f$x@f$-coordinates
!>@param Y             : @f$y@f$-coordinates
!>@param Z             : @f$z@f$-coordinates
!***********************************************************************************************************************************************************************************
   function VTK_GEO_XML_UNST_R4(NN,NC,X,Y,Z) result(E_IO)
!***********************************************************************************************************************************************************************************

      implicit none
    
      !--------------------------------------------------------------------------------------------------------------------------------
      integer(I4P), intent(IN) :: NN       ! number of nodes
      integer(I4P), intent(IN) :: NC       ! number of cells
      real   (R4P), intent(IN) :: X(1:NN)  ! x coordinates
      real   (R4P), intent(IN) :: Y(1:NN)  ! y coordinates
      real   (R4P), intent(IN) :: Z(1:NN)  ! z coordinates
    
      integer(I4P)             :: E_IO     ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      integer(I4P)             :: n1       ! counter
    
      character(len=maxlen)    :: s_buffer ! buffer string
      !--------------------------------------------------------------------------------------------------------------------------------
    
      select case(f_out)
    
         case(f_out_ascii)
           write(unit=Unit_VTK,fmt='(A,'//FI4P//',A,'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece NumberOfPoints="',NN,'" NumberOfCells="',NC,'">'
           indent = indent + 2
           write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//'<Points>'
           indent = indent + 2
           write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" NumberOfComponents="3" Name="Point" format="ascii">'
           write(unit=Unit_VTK,fmt='(3'//FR4P//')',                iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
           write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//'</DataArray>'
           indent = indent - 2
           write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//'</Points>'
    
         case(f_out_binary)
           write(s_buffer,fmt='(A,'//FI4P//',A,'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece NumberOfPoints="',NN,'" NumberOfCells="',NC,'">'
           indent = indent + 2
           write(unit=Unit_VTK,                               iostat=E_IO)trim(s_buffer)//end_rec
           write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//'<Points>'//end_rec
           indent = indent + 2
           write(s_buffer,fmt='(I8)',                         iostat=E_IO)ioffset
           write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" NumberOfComponents="3" Name="Point" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = 3*NN*sizeof(Tipo_R4)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,                        iostat=E_IO)N_Byte,'R4',3*NN
           write(unit=Unit_VTK_Append,                        iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
           write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
           indent = indent - 2
           write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//'</Points>'//end_rec
    
      endselect
    
      return
!***********************************************************************************************************************************************************************************
   endfunction VTK_GEO_XML_UNST_R4
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!The VTK\_GEO\_XML\_CLOSE function is used for closing mesh block data.
!***********************************************************************************************************************************************************************************
   function VTK_GEO_XML_CLOSEP() result(E_IO)
!***********************************************************************************************************************************************************************************

      implicit none
    
      !--------------------------------------------------------------------------------------------------------------------------------
      integer(I4P) :: E_IO ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      !--------------------------------------------------------------------------------------------------------------------------------
    
      indent = indent - 2
    
      select case(f_out)
    
         case(f_out_ascii)
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</Piece>'
    
         case(f_out_binary)
           write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'</Piece>'//end_rec
    
      endselect
    
      return
!***********************************************************************************************************************************************************************************
   endfunction VTK_GEO_XML_CLOSEP
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!The VTK\_CON\_XML function must be called when unstructured grid topology is used. It saves the connectivity of the unstructured mesh.
!>@param NC        : total number of cells 
!>@param connect   : mesh connectivity array
!>@param offset    : cell offset
!>@param cell_type : VTK cell type
!***********************************************************************************************************************************************************************************
   function VTK_CON_XML(NC,connect,offset,cell_type) result(E_IO)
!***********************************************************************************************************************************************************************************

      implicit none
    
      !--------------------------------------------------------------------------------------------------------------------------------
      integer(I4P), intent(IN) :: NC              ! number of cells
      integer(I4P), intent(IN) :: connect(:)      ! mesh connectivity
      integer(I4P), intent(IN) :: offset(1:NC)    ! cell offset
      integer(I1P), intent(IN) :: cell_type(1:NC) ! VTK cell type
    
      integer(I4P)             :: E_IO            ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      integer(I4P)             :: n1              ! counter
    
      character(len=maxlen)    :: s_buffer        ! buffer string
      !--------------------------------------------------------------------------------------------------------------------------------

      select case(f_out)
    
         case(f_out_ascii)
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<Cells>'
           indent = indent + 2
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="connectivity" format="ascii">'
           write(unit=Unit_VTK,fmt=FI4P, iostat=E_IO)(connect(n1),n1=1,size(connect))
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="offsets" format="ascii">'
           write(unit=Unit_VTK,fmt=FI4P, iostat=E_IO)(offset(n1),n1=1,NC)
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int8" Name="types" format="ascii">'
           write(unit=Unit_VTK,fmt=FI1P, iostat=E_IO)(cell_type(n1),n1=1,NC)
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
           indent = indent - 2
           write(unit=Unit_VTK,fmt='(A)', iostat=E_IO)repeat(' ',indent)//'</Cells>'
    
         case(f_out_binary)
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<Cells>'//end_rec
           indent = indent + 2
           write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="connectivity" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = size(connect)*sizeof(Tipo_I4)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I4',size(connect)
           write(unit=Unit_VTK_Append,iostat=E_IO)(connect(n1),n1=1,size(connect))
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
           write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="offsets" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = NC*sizeof(Tipo_I4)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I4',NC
           write(unit=Unit_VTK_Append,iostat=E_IO)(offset(n1),n1=1,NC)
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
           write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int8" Name="types" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = NC*sizeof(Tipo_I1)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I1',NC
           write(unit=Unit_VTK_Append,iostat=E_IO)(cell_type(n1),n1=1,NC)
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
           indent = indent - 2
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</Cells>'//end_rec
    
      endselect
    
      return
!***********************************************************************************************************************************************************************************
   endfunction VTK_CON_XML
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!The VTK\_DAT\_XML function opens or closes CellData/PointData structures.
!>@param var_location     : location of saving variables: CELL for cell-centered, NODE for node-centered
!>@param var_block_action : variables block action: OPEN or CLOSE block
!***********************************************************************************************************************************************************************************
   function VTK_DAT_XML(var_location,var_block_action) result(E_IO)
!***********************************************************************************************************************************************************************************

      implicit none
    
      !--------------------------------------------------------------------------------------------------------------------------------
      character(*), intent(IN) :: var_location     ! location of saving variables: CELL for cell-centered, NODE for node-centered
      character(*), intent(IN) :: var_block_action ! variables block action: OPEN or CLOSE block
      integer(I4P)             :: E_IO             ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      !--------------------------------------------------------------------------------------------------------------------------------

      select case(f_out)
    
         case(f_out_ascii)
           select case(trim(var_location))
              case('CELL')
                select case(trim(var_block_action))
                   case('OPEN')
                     write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<CellData>'
                     indent = indent + 2
                   case('CLOSE')
                     indent = indent - 2
                     write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</CellData>'
                endselect
              case('NODE')
                select case(trim(var_block_action))
                   case('OPEN')
                     write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<PointData>'
                     indent = indent + 2
                   case('CLOSE')
                     indent = indent - 2
                     write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</PointData>'
                endselect
           endselect
    
         case(f_out_binary)
           select case(trim(var_location))
              case('CELL')
                select case(trim(var_block_action))
                   case('OPEN')
                     write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'<CellData>'//end_rec
                     indent = indent + 2
                   case('CLOSE')
                     indent = indent - 2
                     write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'</CellData>'//end_rec
                endselect
              case('NODE')
                select case(trim(var_block_action))
                   case('OPEN')
                     write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'<PointData>'//end_rec
                     indent = indent + 2
                   case('CLOSE')
                     indent = indent - 2
                     write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'</PointData>'//end_rec
                endselect
           endselect
    
      endselect
    
      return
!***********************************************************************************************************************************************************************************
  endfunction VTK_DAT_XML
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!The VTK\_VAR\_XML\_SCAL\_R4 function saves R4P scalar variables.
!>@param NC_NN   : number of cells or nodes
!>@param varname : variable name
!>@param var     : values to be saved
!***********************************************************************************************************************************************************************************
   function VTK_VAR_XML_SCAL_R4(NC_NN,varname,var) result(E_IO)
!***********************************************************************************************************************************************************************************
    
      implicit none
    
      !--------------------------------------------------------------------------------------------------------------------------------
      integer(I4P), intent(IN) :: NC_NN        ! number of cells or nodes
      character(*), intent(IN) :: varname      ! variable name
      real   (R4P), intent(IN) :: var(1:NC_NN) ! variable to be saved

      integer(I4P)             :: E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      integer(I4P)             :: n1           ! counter
    
      character(len=maxlen)    :: s_buffer     ! buffer string
      !--------------------------------------------------------------------------------------------------------------------------------
    
      select case(f_out)
    
         case(f_out_ascii)
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="'//trim(varname)//'" NumberOfComponents="1" format="ascii">'
           write(unit=Unit_VTK,fmt=FR4P, iostat=E_IO)var
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    
         case(f_out_binary)
           write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="'//trim(varname)//'" NumberOfComponents="1" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = NC_NN*sizeof(Tipo_R4)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'R4',NC_NN
           write(unit=Unit_VTK_Append,iostat=E_IO)(var(n1),n1=1,NC_NN)
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    
      endselect
    
      return
!***********************************************************************************************************************************************************************************
   endfunction VTK_VAR_XML_SCAL_R4
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!The VTK\_VAR\_XML\_SCAL\_I4 function saves I4P scalar variables.
!>@param NC_NN   : number of cells or nodes
!>@param varname : variable name
!>@param var     : values to be saved
!***********************************************************************************************************************************************************************************
   function VTK_VAR_XML_SCAL_I4(NC_NN,varname,var) result(E_IO)
!***********************************************************************************************************************************************************************************

      implicit none
    
      !--------------------------------------------------------------------------------------------------------------------------------
      integer(I4P), intent(IN) :: NC_NN        ! number of cells or nodes
      character(*), intent(IN) :: varname      ! variable name
      integer(I4P), intent(IN) :: var(1:NC_NN) ! variable to be saved
    
      integer(I4P)             :: E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      integer(I4P)             :: n1           ! counter
    
      character(len=maxlen)    :: s_buffer     ! buffer string
      !--------------------------------------------------------------------------------------------------------------------------------
    
      select case(f_out)
    
         case(f_out_ascii)
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="'//trim(varname)//'" NumberOfComponents="1" format="ascii">'
           write(unit=Unit_VTK,fmt=FI4P, iostat=E_IO)var
           write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    
         case(f_out_binary)
           write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="'//trim(varname)//'" NumberOfComponents="1" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = NC_NN*sizeof(Tipo_I4)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I4',NC_NN
           write(unit=Unit_VTK_Append,iostat=E_IO)(var(n1),n1=1,NC_NN)
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    
      endselect
    
      return
!***********************************************************************************************************************************************************************************
  endfunction VTK_VAR_XML_SCAL_I4
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!The VTK\_VAR\_XML\_VECT\_R4 function saves R4P vectorial variables.
!>@param NC_NN   : number of cells or nodes
!>@param varname : variable name
!>@param varX    : @f$x@f$-component values to be saved
!>@param varY    : @f$y@f$-component values to be saved
!>@param varZ    : @f$z@f$-component values to be saved
!***********************************************************************************************************************************************************************************
   function VTK_VAR_XML_VECT_R4(NC_NN,varname,varX,varY,varZ) result(E_IO)
!***********************************************************************************************************************************************************************************

      implicit none

      !--------------------------------------------------------------------------------------------------------------------------------
      integer(I4P), intent(IN) :: NC_NN         ! number of cells or nodes
      character(*), intent(IN) :: varname       ! variable name
      real   (R4P), intent(IN) :: varX(1:NC_NN) ! x component
      real   (R4P), intent(IN) :: varY(1:NC_NN) ! y component
      real   (R4P), intent(IN) :: varZ(1:NC_NN) ! z component

      integer(I4P)             :: E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      integer(I4P)             :: n1            ! counter

      character(len=maxlen)    :: s_buffer      ! buffer string
      !--------------------------------------------------------------------------------------------------------------------------------

      select case(f_out)

         case(f_out_ascii)
           write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="'//trim(varname)//'" NumberOfComponents="3" format="ascii">'
           write(unit=Unit_VTK,fmt='(3'//FR4P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
           write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'</DataArray>'

         case(f_out_binary)
           write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="'//trim(varname)//'" NumberOfComponents="3" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = 3*NC_NN*sizeof(Tipo_R4)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'R4',3*NC_NN
           write(unit=Unit_VTK_Append,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec

      endselect

      return
!***********************************************************************************************************************************************************************************
  endfunction VTK_VAR_XML_VECT_R4
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!The VTK\_VAR\_XML\_VECT\_I4 function saves I4P vectorial variables.
!>@param NC_NN   : number of cells or nodes
!>@param varname : variable name
!>@param varX    : @f$x@f$-component values to be saved
!>@param varY    : @f$y@f$-component values to be saved
!>@param varZ    : @f$z@f$-component values to be saved
!***********************************************************************************************************************************************************************************
   function VTK_VAR_XML_VECT_I4(NC_NN,varname,varX,varY,varZ) result(E_IO)
!***********************************************************************************************************************************************************************************

      implicit none
    
      !--------------------------------------------------------------------------------------------------------------------------------
      integer(I4P), intent(IN) :: NC_NN         ! number of cells or nodes
      character(*), intent(IN) :: varname       ! variable name
      integer(I4P), intent(IN) :: varX(1:NC_NN) ! x component
      integer(I4P), intent(IN) :: varY(1:NC_NN) ! y component
      integer(I4P), intent(IN) :: varZ(1:NC_NN) ! z component
    
      integer(I4P)             :: E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      integer(I4P)             :: n1            ! counter
    
      character(len=maxlen)    :: s_buffer      ! buffer string
      !--------------------------------------------------------------------------------------------------------------------------------
    
      select case(f_out)
    
         case(f_out_ascii)
           write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="'//trim(varname)//'" NumberOfComponents="3" format="ascii">'
           write(unit=Unit_VTK,fmt='(3'//FI4P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
           write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    
         case(f_out_binary)
           write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="'//trim(varname)//'" NumberOfComponents="3" format="appended" offset="',trim(s_buffer),'">'//end_rec
           N_Byte  = 3*NC_NN*sizeof(Tipo_I4)
           ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
           write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I4',3*NC_NN
           write(unit=Unit_VTK_Append,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
           write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    
      endselect
    
      return
!***********************************************************************************************************************************************************************************
   endfunction VTK_VAR_XML_VECT_I4
!***********************************************************************************************************************************************************************************

!
!
!>@brief
!!The VTK\_END\_XML function finalizes opened files.
!***********************************************************************************************************************************************************************************
   function VTK_END_XML() result(E_IO)
!***********************************************************************************************************************************************************************************

      implicit none
    
      !--------------------------------------------------------------------------------------------------------------------------------
      integer(I4P)              :: E_IO      ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
      character(2)              :: var_type  ! var\_type = R4,I4
      real   (R4P), allocatable :: v_R4(:)   ! R4 vector for IO in AppendData
      integer(I4P), allocatable :: v_I4(:)   ! I4 vector for IO in AppendData
      integer(I1P), allocatable :: v_I1(:)   ! I4 vector for IO in AppendData
      integer(I4P)              :: N_v       ! vector dimension
      integer(I4P)              :: n1        ! counter
      !--------------------------------------------------------------------------------------------------------------------------------
    
      select case(f_out)
    
      case(f_out_ascii)
    
        indent = indent - 2
        write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</'//trim(topology)//'>'
        write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'</VTKFile>'
    
      case(f_out_binary)
    
        indent = indent - 2
        write(unit  =Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</'//trim(topology)//'>'//end_rec
        write(unit  =Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<AppendedData encoding="raw">'//end_rec
        write(unit  =Unit_VTK,       iostat=E_IO)'_'
        endfile(unit=Unit_VTK_Append,iostat=E_IO)
        rewind(unit =Unit_VTK_Append,iostat=E_IO)
    
        do
    
          read(unit=Unit_VTK_Append,iostat=E_IO,end=100)N_Byte,var_type,N_v
    
          select case(var_type)
    
             case('R4')
               allocate(v_R4(1:N_v))
               read (unit=Unit_VTK_Append,iostat=E_IO)(v_R4(n1),n1=1,N_v)
               write(unit=Unit_VTK,       iostat=E_IO)N_Byte,(v_R4(n1),n1=1,N_v)
               deallocate(v_R4)
    
             case('I4')
               allocate(v_I4(1:N_v))
               read (unit=Unit_VTK_Append,iostat=E_IO)(v_I4(n1),n1=1,N_v)
               write(unit=Unit_VTK,       iostat=E_IO)N_Byte,(v_I4(n1),n1=1,N_v)
               deallocate(v_I4)
    
             case('I1')
               allocate(v_I1(1:N_v))
               read (unit=Unit_VTK_Append,iostat=E_IO)(v_I1(n1),n1=1,N_v)
               write(unit=Unit_VTK,       iostat=E_IO)N_Byte,(v_I1(n1),n1=1,N_v)
               deallocate(v_I1)
    
          endselect
    
        enddo
    
        100 continue
        write(unit=Unit_VTK,iostat=E_IO)end_rec
        write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'</AppendedData>'//end_rec
        write(unit=Unit_VTK,iostat=E_IO)'</VTKFile>'//end_rec

        !->closing AppendData file
        close(unit=Unit_VTK_Append,iostat=E_IO)

      endselect
    
      close(unit=Unit_VTK,iostat=E_IO)
    
      return
!***********************************************************************************************************************************************************************************
   endfunction VTK_END_XML
!***********************************************************************************************************************************************************************************

end module mod_vtk_io
