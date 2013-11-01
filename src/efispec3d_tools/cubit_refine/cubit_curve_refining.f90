!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                  !
!Author : Florent DE MARTIN, July 2013                                             !
!                                                                                  !
!Contact: f.demartin at brgm.fr                                                    !
!                                                                                  !
!This program is free software: you can redistribute it and/or modify it under     !
!the terms of the GNU General Public License as published by the Free Software     !
!Foundation, either version 3 of the License, or (at your option) any later        !
!version.                                                                          !
!                                                                                  !
!This program is distributed in the hope that it will be useful, but WITHOUT ANY   !
!WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   !
!PARTICULAR PURPOSE. See the GNU General Public License for more details.          !
!                                                                                  !
!You should have received a copy of the GNU General Public License along with      !
!this program. If not, see http://www.gnu.org/licenses/.                           !
!                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cubit_curve_refining
   
   implicit none

   integer, parameter   :: NBLOCKMAX         = 50
   integer, parameter   :: N_GEOM_NODE_HEXA  = 8
   integer, parameter   :: N_GEOM_NODE_LOC   = 2
   integer, parameter   :: N_HEXA_INSIDE_MAX = 1000000
   integer, parameter   :: N_HEXA_MAX        = 1000000
   
   real, allocatable    :: x_geom_node(:) !>X coordinates of geometrical nodes
   real, allocatable    :: y_geom_node(:) !>Y coordinates of geometrical nodes
   real, allocatable    :: z_geom_node(:) !>Z coordinates of geometrical nodes

   real                 :: x_node_hexa(N_GEOM_NODE_HEXA)
   real                 :: y_node_hexa(N_GEOM_NODE_HEXA)
   real                 :: z_node_hexa(N_GEOM_NODE_HEXA)

   real                 :: xtmp
   real                 :: ytmp
   real                 :: ztmp

   real                 :: zmin

   real                 :: x_centroid
   real                 :: y_centroid
   real                 :: z_centroid

   real, allocatable    :: xcurve(:),ycurve(:)


   integer              :: get_local_node_from_global
   integer, allocatable :: geom_node_local_to_global(:)
   integer, allocatable :: hexa_local_to_global(:)

   integer              :: xi_loc_gnode_hexa(N_GEOM_NODE_HEXA)
   integer              :: et_loc_gnode_hexa(N_GEOM_NODE_HEXA)
   integer              :: ze_loc_gnode_hexa(N_GEOM_NODE_HEXA)

   integer              :: hexa_inside(N_HEXA_INSIDE_MAX)

   real                 :: co_loc_gnode_hexa(N_GEOM_NODE_LOC)
   real                 :: abscissa(N_GEOM_NODE_LOC,N_GEOM_NODE_LOC)

   integer              :: ndata_curve
   integer              :: n_hexa_inside
   integer              :: number_of_nodes
   integer              :: nblock
   integer              :: iblock
   integer              :: len_code
   integer              :: inode
   integer              :: ihexa
   integer              :: number_of_hexa
   integer              :: local_node

   integer              :: is_in
   integer              :: is_out

   integer              :: number_of_hexa_block(NBLOCKMAX)
   integer              :: material_of_element(NBLOCKMAX)
   
   integer              :: node_hexa(N_GEOM_NODE_HEXA)
   
   character(len= 90)   :: prefix
   character(len= 90)   :: file_inp
   character(len= 90)   :: file_cur
   character(len= 90)   :: czmin
   character(len= 90)   :: myformat
   character(len= 20)   :: code1
   character(len= 20)   :: code2
   character(len= 20)   :: code3
   character(len= 20)   :: code4
   character(len= 15)   :: cn
   character(len=  1)   :: type_of_element(NBLOCKMAX)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   interface

      subroutine set_element_convention(xiloc,etloc,zeloc,coordloc,absci)
         implicit none
         integer, intent(out) :: xiloc(:)
         integer, intent(out) :: etloc(:)
         integer, intent(out) :: zeloc(:)
         real   , intent(out) :: coordloc(:)
         real   , intent(out) :: absci(:,:)
      end subroutine set_element_convention
 
      subroutine read_curve(f,n,x,y)
         implicit none
         character(len= 90), intent( in) :: f
         integer           , intent(out) :: n
         real, allocatable , intent(out) :: x(:)
         real, allocatable , intent(out) :: y(:)
      end subroutine read_curve

      subroutine get_hexa_centroid(xiloc,etloc,zeloc,absci,coloc,xn,yn,zn,x,y,z)
         implicit none
         integer, intent( in) :: xiloc(:)
         integer, intent( in) :: etloc(:)
         integer, intent( in) :: zeloc(:)
         real   , intent( in) :: absci(:,:)
         real   , intent( in) :: coloc(:)
         real   , intent( in) :: xn(:)
         real   , intent( in) :: yn(:)
         real   , intent( in) :: zn(:)
         real   , intent(out) :: x
         real   , intent(out) :: y
         real   , intent(out) :: z
      end subroutine get_hexa_centroid

   end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   call getarg(1,file_inp)
   call getarg(2,file_cur)
   call getarg(3,czmin)
   read(czmin,'(F10.0)') zmin

   prefix = file_inp(1:len_trim(file_inp)-4)

   write(*,*) "inp   file = ",trim(file_inp)
   write(*,*) "curve file = ",trim(file_cur)
   write(*,*) "zmin       = ",zmin

!
!->initialize variable 
   call set_element_convention(xi_loc_gnode_hexa,et_loc_gnode_hexa,ze_loc_gnode_hexa,co_loc_gnode_hexa,abscissa)
   
   number_of_nodes    = 0
   number_of_hexa_block(:) = 0
   nblock             = 0
   number_of_hexa     = 0

!
!->read file that contains close curve
   call read_curve(file_cur,ndata_curve,xcurve,ycurve)
   
!
!->first pass on inp file to count number of nodes and number of hexa
   open(unit=1,file=trim(adjustl(file_inp)),status='old')
   do while (.true.)
   
      read(unit=1,fmt=*) code1
   
      select case (trim(adjustl(code1)))
   
         case("*NODE")
            write(*,'(A)') "Counting number of nodes..."
            do while(.true.)
               read(unit=1,fmt=*) code2
               if (trim(adjustl(code2(1:1))) /= "*") then 
                  number_of_nodes = number_of_nodes + 1
               else
                  exit
               endif
            enddo
            close(2)
            write(*,'(A,I15)') "Number of nodes = ",number_of_nodes
            write(*,*) " "
            backspace(1)
   
         case("*ELEMENT")

            backspace(1)
            read(unit=1,fmt=*) code2,code3,code4

            nblock = nblock + 1

            write(*,'(A,I6)') "Counting number of elements in nblock ",nblock

            len_code                = len_trim(code4)
            type_of_element(nblock) = code4(len_code-2:len_code-2)

            if (type_of_element(nblock) == "l") read(code4(len_code-1:len_code),'(I2)') material_of_element(nblock)
            write(*,'(2A)')   "Name of the nblock ",code4(len_code-2:len_code)

            if (type_of_element(nblock) == "l") write(*,'(A,I2)') "Material number of the nblock ",material_of_element(nblock)

            do while(.true.)

               read(unit=1,fmt=*) code2
               if (trim(adjustl(code2(1:1))) /= "*") then 
                  number_of_hexa_block(nblock) = number_of_hexa_block(nblock) + 1
               else
                  exit
               endif

            enddo
   
            write(*,'(A,I15)') "Number of elements = ",number_of_hexa_block(nblock)

            if (type_of_element(nblock) == "l") then 
               number_of_hexa = number_of_hexa + number_of_hexa_block(nblock)
            endif

            write(*,*) " "
            backspace(1)
   
         case("*SOLID")
            write(*,*) "Searching hexa inside closed curve..."
            exit
   
         case default
   
      end select 
   
   enddo
   rewind(1)
   
   
!
!->Arrays allocation
   allocate(x_geom_node(number_of_nodes))
   allocate(y_geom_node(number_of_nodes))
   allocate(z_geom_node(number_of_nodes))
   allocate(geom_node_local_to_global(number_of_nodes))

!
!->second pass on inp file to fill x,y,z coordinates and hexa connectivity
!->find hexa inside closed curve

   iblock        = 0
   n_hexa_inside = 0

   do while (.true.)
   
      read(unit=1,fmt=*) code1
   
      select case (trim(adjustl(code1)))
   
         case("*NODE")

            do inode = 1,number_of_nodes

               read (unit=1,fmt=*) geom_node_local_to_global(inode),xtmp,ytmp,ztmp

               x_geom_node(inode) = xtmp
               y_geom_node(inode) = ytmp
               z_geom_node(inode) = ztmp

            enddo
   
         case("*ELEMENT")
   
            iblock = iblock + 1
   
            if (type_of_element(iblock) == "l") then

               allocate(hexa_local_to_global(number_of_hexa_block(iblock)))
   
               do ihexa = 1,number_of_hexa_block(iblock)

   
                  read (unit=1,fmt=*) hexa_local_to_global(ihexa),(node_hexa(inode),inode=1,N_GEOM_NODE_HEXA)

                  do inode = 1,N_GEOM_NODE_HEXA

                     
                     local_node = get_local_node_from_global(node_hexa(inode),geom_node_local_to_global,number_of_nodes)

                     x_node_hexa(inode) = x_geom_node(local_node)
                     y_node_hexa(inode) = y_geom_node(local_node)
                     z_node_hexa(inode) = z_geom_node(local_node)

                  enddo

                  call get_hexa_centroid(                  &
                                         xi_loc_gnode_hexa &
                                        ,et_loc_gnode_hexa &
                                        ,ze_loc_gnode_hexa &
                                        ,abscissa          &
                                        ,co_loc_gnode_hexa &
                                        ,x_node_hexa       &
                                        ,y_node_hexa       &
                                        ,z_node_hexa       &
                                        ,x_centroid        &
                                        ,y_centroid        &
                                        ,z_centroid)

                  if (z_centroid >= zmin) then
                 
                     call in_poly_curve(x_centroid,y_centroid,xcurve,ycurve,ndata_curve,is_in,is_out)
                 
                     if(is_in >= 0) then
                 
                        n_hexa_inside = n_hexa_inside + 1
                 
                        if (n_hexa_inside > N_HEXA_INSIDE_MAX) then
                           stop "Code stopped: Too much hexa inside curve"
                        else
                           hexa_inside(n_hexa_inside) = hexa_local_to_global(ihexa)
                        endif
                 
                     endif
                 
                  endif

               enddo

               deallocate(hexa_local_to_global)

            endif
   
         case("*SOLID")
            write(*,*) "Done"
            exit
   
         case default
   
      end select

   enddo
   close(1)

   write(*,*) "Number of hexa found inside curve = ",n_hexa_inside

!
!->write hexa inside closed curve for cubit
   write(cn,'(i15)') n_hexa_inside
   myformat = "(a,"//trim(adjustl(cn))//"(i0,1X))"

   open (unit=2,file="curve_refining.jou")
   write(unit=2,fmt =trim(adjustl(myformat))) "refine hex ",(hexa_inside(ihexa),ihexa=1,n_hexa_inside)
   close(2)

   stop

end program cubit_curve_refining

integer function get_local_node_from_global(target_num,a,n)

   implicit none

   integer,               intent(in) :: n
   integer, dimension(n), intent(in) :: a
   integer,               intent(in) :: target_num

   integer                           :: i

   do i = 1,n

      if (a(i) == target_num) then
         get_local_node_from_global = i
         exit
      endif

   enddo

end function get_local_node_from_global
