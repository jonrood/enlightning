c ===========================================================================
c Copyright (C) 2017 Jon Rood.
c
c This file is part of Enlightning source code.
c 
c Enlightning source code is free software; you can redistribute it
c and/or modify it under the terms of the GNU General Public License as
c published by the Free Software Foundation; either version 3 of the License,
c or (at your option) any later version.
c 
c Enlightning source code is distributed in the hope that it will be
c useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c 
c You should have received a copy of the GNU General Public License
c along with Enlightning; if not, see <http://www.gnu.org/licenses/>.
c ===========================================================================

define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

c This function calculates the smoothness or gradient
c in order to tag cells that require refinement.

      subroutine tag_cells(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  tags,
     &  w_n, H0,
     &  refine_tag_val,
     &  tolerance,
     &  nequ)
c***********************************************************************
c***********************************************************************     
      implicit none
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      integer refine_tag_val
      integer nequ
      REAL    tolerance(0:nequ-1)
      REAL    tmp1,tmp2,xi,ep
      REAL    rx,ry
c
c variables in 2d cell indexed         
      integer 
     &     tags(CELL2d(ifirst,ilast,0))
      REAL
     &     w_n(CELL2dVECG(ifirst,ilast,gcw),0:nequ-1),
     &     H0(CELL2dVECG(ifirst,ilast,gcw))
c
      integer i,j
c
c***********************************************************************     

      do j=ifirst1,ilast1
         do i=ifirst0,ilast0
            tags(i,j) = 0
          end do
       end do

      xi=1
      ep=0.9*tolerance(0)*(xi**2)/(1-0.9*(tolerance(0)**2))

      do j=ifirst1,ilast1
         do i=ifirst0,ilast0
            tmp1 = H0(i,j)-H0(i-1,j)
            tmp2 = H0(i+1,j)-H0(i,j)
            rx = (2*abs(tmp1*tmp2)+ep)/(tmp1**2+tmp2**2+ep)
            if (rx .lt. 0.9999)
     &         tags(i,j) = refine_tag_val
          end do
       end do

      do j=ifirst1,ilast1
         do i=ifirst0,ilast0
            tmp1 = H0(i,j)-H0(i,j-1)
            tmp2 = H0(i,j+1)-H0(i,j)
            ry = (2*abs(tmp1*tmp2)+ep)/(tmp1**2+tmp2**2+ep)
            if (ry .lt. 0.9999)
     &         tags(i,j) = refine_tag_val
          end do
       end do

      xi=1
      ep=0.9*tolerance(1)*(xi**2)/(1-0.9*(tolerance(1)**2))

      do j=ifirst1,ilast1
         do i=ifirst0,ilast0
            tmp1 = w_n(i,j,0)-w_n(i-1,j,0)
            tmp2 = w_n(i+1,j,0)-w_n(i,j,0)
            rx = (2*abs(tmp1*tmp2)+ep)/(tmp1**2+tmp2**2+ep)
            if (rx .lt. 0.9999)
     &         tags(i,j) = refine_tag_val
          end do
       end do

      do j=ifirst1,ilast1
         do i=ifirst0,ilast0
            tmp1 = w_n(i,j,0)-w_n(i,j-1,0)
            tmp2 = w_n(i,j+1,0)-w_n(i,j,0)
            ry = (2*abs(tmp1*tmp2)+ep)/(tmp1**2+tmp2**2+ep)
            if (ry .lt. 0.9999)
     &         tags(i,j) = refine_tag_val
          end do
       end do

      return
      end

