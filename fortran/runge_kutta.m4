c ===========================================================================
c Copyright (C) 2012 Jon Rood.
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
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

c This function performs the runge kutta.

      subroutine runge_kutta(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  dt, alpha_1, alpha_2, beta,
     &  w_n,
     &  ww,
     &  KK,
     &  nequ)
c***********************************************************************
      implicit none
cinclude(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      integer nequ

      REAL    dt, alpha_1, alpha_2, beta
c
c variables in 2d cell indexed
      REAL
     &     w_n(CELL2dVECG(ifirst,ilast,gcw),0:nequ-1),
     &     ww(CELL2d(ifirst,ilast,0),0:nequ-1),
     &     KK(CELL2d(ifirst,ilast,0),0:nequ-1)
c
c***********************************************************************
c***********************************************************************
c
      integer i,j,k

      do k=0,nequ-1
         do j=ifirst1,ilast1
            do i=ifirst0,ilast0
               w_n(i,j,k) = 
     &            alpha_1*ww(i,j,k) + 
     &            alpha_2*w_n(i,j,k) +
     &            beta*dt*KK(i,j,k)
            end do
         end do
      end do

      return
      end
