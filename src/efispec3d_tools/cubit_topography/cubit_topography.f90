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

program cubit_topography
   
   implicit none

   real, allocatable    :: x_geom_node(:) !>X coordinates of geometrical nodes
   real, allocatable    :: y_geom_node(:) !>Y coordinates of geometrical nodes
   real, allocatable    :: z_geom_node(:) !>Z coordinates of geometrical nodes

   real                 :: x_node
   real                 :: y_node
   real                 :: z_node

   real                 :: xmin
   real                 :: xmax
   real                 :: ymin
   real                 :: ymax
   real                 :: zmin
   real                 :: zmax

   real                 :: x0_topo
   real                 :: y0_topo
   real                 :: dx_topo
   real                 :: dy_topo
   real                 :: z_max_apply_topo
   real                 :: z_min_apply_topo

   real, allocatable    :: x_topo(:)
   real, allocatable    :: y_topo(:)
   real, allocatable    :: z_topo(:,:)

   integer              :: nx_topo
   integer              :: ny_topo
   integer              :: number_of_nodes
   integer              :: inode
   integer              :: ipos_node
   integer              :: ios
   integer              :: narg

   character(len= 90)   :: prefix
   character(len= 90)   :: line
   character(len= 90)   :: file_inp
   character(len= 90)   :: file_topo
   character(len= 90)   :: prog_name
   character(len= 20)   :: code1
   character(len= 20)   :: code2
   character(len= 20)   :: czmax_hexa
   character(len= 20)   :: czmin_hexa
   character(len=  1)   :: ctmp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   interface

      subroutine read_topography(fname,x,y,z,dx,dy,x0,y0,nx,ny)

         implicit none

         character(len=90), intent( in) :: fname
         real, allocatable, intent(out) :: x(:)
         real, allocatable, intent(out) :: y(:)
         real, allocatable, intent(out) :: z(:,:)
         real             , intent(out) :: dx
         real             , intent(out) :: dy
         real             , intent(out) :: x0
         real             , intent(out) :: y0
         integer          , intent(out) :: nx
         integer          , intent(out) :: ny

      end subroutine read_topography


      subroutine compute_z(x_node,y_node,z_node,x_topo,y_topo,z_topo,nx_topo,ny_topo,z_max,z_min)
      
         implicit none
      
         real             , intent(in   ) :: x_node
         real             , intent(in   ) :: y_node
         real             , intent(inout) :: z_node
         integer          , intent(in   ) :: nx_topo
         integer          , intent(in   ) :: ny_topo
         real             , intent(in   ) :: x_topo(nx_topo)
         real             , intent(in   ) :: y_topo(ny_topo)
         real             , intent(in   ) :: z_topo(nx_topo,ny_topo)
         real             , intent(in   ) :: z_min  
         real             , intent(in   ) :: z_max  

      end subroutine compute_z
 
   end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   narg = command_argument_count()
   if (narg /= 4) then
      call get_command_argument(0,prog_name)
      write(*,*) "Usage : "//trim(adjustl(prog_name))//" file.inp file.dem zmax zmin"
      stop
   endif
 
   call get_command_argument(1,file_inp)
   call get_command_argument(2,file_topo)
   call get_command_argument(3,czmax_hexa)
   call get_command_argument(4,czmin_hexa)


   read(czmax_hexa,*) z_max_apply_topo
   read(czmin_hexa,*) z_min_apply_topo

   write(*,*) ""
   write(*,*) "Topography will be applyied on nodes between"
   write(*,*) "z_max = ",z_max_apply_topo
   write(*,*) "z_min = ",z_min_apply_topo
   write(*,*) ""

   prefix = file_inp(1:len_trim(file_inp)-4)

!
!->read file that contains topography
   call read_topography(file_topo,x_topo,y_topo,z_topo,dx_topo,dy_topo,x0_topo,y0_topo,nx_topo,ny_topo)
   
!
!->first pass on inp file to count number of nodes
   number_of_nodes = 0
   open(unit=1,file=trim(adjustl(file_inp)),status='old')
   do while (.true.)
   
      read(unit=1,fmt=*) code1
   
      select case (trim(adjustl(code1)))
   
         case("*NODE")
            write(*,*) "Counting number of nodes..."
            do while(.true.)
               read(unit=1,fmt=*) code2
               if (trim(adjustl(code2(1:1))) /= "*") then 
                  number_of_nodes = number_of_nodes + 1
               else
                  exit
               endif
            enddo
            close(2)
            write(*,*) "Number of nodes = ",number_of_nodes
            write(*,*) " "
            backspace(1)
   
         case("*ELEMENT")
            write(*,*) "Reading geometric node coordinates..."
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

!
!->second pass in inp file to read geometric node coordinates
   do while (.true.)
   
      read(unit=1,fmt=*) code1
   
      select case (trim(adjustl(code1)))
   
         case("*NODE")

            do inode = 1,number_of_nodes

               read (unit=1,fmt=*) ipos_node,x_node,y_node,z_node

               x_geom_node(ipos_node) = x_node
               y_geom_node(ipos_node) = y_node
               z_geom_node(ipos_node) = z_node

            enddo
   
         case("*ELEMENT")
            write(*,*) "Applying topography..."
            exit
   
         case default
   
      endselect

   enddo
   rewind(1)

!
!->compute new z that include topo
   do inode = 1,number_of_nodes

      x_node = x_geom_node(inode)
      y_node = y_geom_node(inode)
      z_node = z_geom_node(inode)

      if ( (z_node <= z_max_apply_topo*1.0001) .and. (z_node >= z_min_apply_topo*1.0001) ) then
         call compute_z(x_node,y_node,z_node,x_topo,y_topo,z_topo,nx_topo,ny_topo,z_max_apply_topo,z_min_apply_topo)
      endif

      z_geom_node(inode) = z_node


   enddo
!
!->find min/max
   xmin      = minval(x_geom_node(:))
   xmax      = maxval(x_geom_node(:))
   ymin      = minval(y_geom_node(:))
   ymax      = maxval(y_geom_node(:))
   zmin      = minval(z_geom_node(:))
   zmax      = maxval(z_geom_node(:))

   write(*,*) "xmin = ",xmin
   write(*,*) "xmax = ",xmax
   write(*,*) "ymin = ",ymin
   write(*,*) "ymax = ",ymax
   write(*,*) "zmin = ",zmin
   write(*,*) "zmax = ",zmax

!
!->write inp file with topography
   open(unit=2,file=trim(adjustl(prefix))//"_tmp.inp")
   do

      read(unit=1,fmt='(a)',iostat=ios) line
      if (ios /= 0) exit

      select case (trim(adjustl(line(1:5))))

         case ("*NODE")
            
            write(unit=2,fmt='(a)') trim(line)
            do inode = 1,number_of_nodes
               !
               !write new inp file with topography
               write(unit=2,fmt='(I8,",",1X,ES15.6,",",1X,ES15.6,",",1X,ES15.6)') inode,x_geom_node(inode),y_geom_node(inode),z_geom_node(inode)
               !
               !read unsed node in old inp file
               read (unit=1,fmt='(a)') ctmp
            enddo

         case default
            write(unit=2,fmt='(a)') trim(line)

      endselect
   enddo
   close(2)
   close(1)
!
!->write journal file to import inp file
   open(unit=4,file=trim(adjustl(prefix))//"_topo.jou")
   write(unit=4,fmt='(3a)') "import abaqus '",trim(adjustl(prefix))//"_tmp.inp","'"

   write(unit=4,fmt='( a)') "delete nodeset all"
   write(unit=4,fmt='( a)') "delete group all"
   write(unit=4,fmt='( a)') "delete block all"

   write(unit=4,fmt='( a)') "compress ids all"

   write(unit=4,fmt='( a)') "block 1 hex all"
   write(unit=4,fmt='( a)') "block 1 name 'l01'"
   write(unit=4,fmt='(5(a,E15.7))') "block 2 face with x_coord >=",xmax," or x_coord <=",xmin," or y_coord >=",ymax," or y_coord <=",ymin," or z_coord <=",zmin
   write(unit=4,fmt='( a)') "block 2 name 'prx'"
   write(unit=4,fmt='( a)') "skin hex all make block 3"
   write(unit=4,fmt='( a)') "block 3 name 'fsu'"

   write(unit=4,fmt='(3a)') "export abaqus '",trim(adjustl(prefix))//"_topo.inp","' block all dimension 3 overwrite"
   close(4)

   stop

end program cubit_topography
