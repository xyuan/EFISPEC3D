subroutine set_element_convention(xiloc,etloc,zeloc,coordloc,absci)

   implicit none
      
   integer, intent(out) :: xiloc(:)
   integer, intent(out) :: etloc(:)
   integer, intent(out) :: zeloc(:)
   real   , intent(out) :: coordloc(:)
   real   , intent(out) :: absci(:,:)

   integer              :: i
   integer              :: j

!
!->local coordinate of geometric nodes
   coordloc(1) = +1.0
   coordloc(2) = -1.0
!
!->local position of geometric nodes (see cubit convention)
   xiloc(1) = 1
   etloc(1) = 1
   zeloc(1) = 1

   xiloc(2) = 1
   etloc(2) = 2
   zeloc(2) = 1

   xiloc(3) = 2
   etloc(3) = 2
   zeloc(3) = 1

   xiloc(4) = 2
   etloc(4) = 1
   zeloc(4) = 1

   xiloc(5) = 1
   etloc(5) = 1
   zeloc(5) = 2

   xiloc(6) = 1
   etloc(6) = 2
   zeloc(6) = 2

   xiloc(7) = 2
   etloc(7) = 2
   zeloc(7) = 2

   xiloc(8) = 2
   etloc(8) = 1
   zeloc(8) = 2

!
!->compute inverse of distance between geometric nodes

   do i = 1,size(coordloc)
      do j = 1,size(coordloc)

         absci(j,i) = 0.0

         if (i /= j) absci(j,i) = 1.0/(coordloc(i) - coordloc(j))

      enddo
   enddo

   return

end subroutine set_element_convention
