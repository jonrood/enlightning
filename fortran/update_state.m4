c ===========================================================================
c Copyright (C) 2015 Jon Rood.
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

c This function updates the pressure and temperature state variables.

      subroutine update_state(
     &  dxy, xylo,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  w_n,
     &  pres,
     &  temp,
     &  rho0,
     &  c,
     &  mygamma,
     &  c_p,
     &  p0,
     &  T0,
     &  c_lr,
     &  c_v,
     &  R_tilde,
     &  M,
     &  R,
     &  nequ)
c***********************************************************************
      implicit none
cinclude(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      REAL dxy(0:1), xylo(0:1)
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      REAL rho0
      REAL c
      REAL mygamma
      REAL c_p
      REAL p0
      REAL T0
      REAL c_lr
      REAL c_v
      REAL R_tilde
      REAL M
      REAL R
      integer nequ
      REAL xc(0:1)
c
c variables in 2d cell indexed
      REAL
     &     w_n(CELL2dVECG(ifirst,ilast,gcw),0:nequ-1),
     &     pres(CELL2dVECG(ifirst,ilast,gcw)),
     &     temp(CELL2dVECG(ifirst,ilast,gcw))
c
c***********************************************************************
c***********************************************************************
c
      integer i,j

c      write(*,100) 'c_p=',c_p
c      write(*,100) 'c=',c
c      write(*,100) 'mygamma=',mygamma
c      write(*,100) 'rho0=',rho0
c      write(*,100) 'p0=',p0
c      write(*,100) 'T0=',T0
c  100 format (A,E8.3)

      do j=ifirst1,ilast1
         xc(1) = xylo(1)+dxy(1)*(dble(j-ifirst1)+0.5)
         do i=ifirst0,ilast0
            pres(i,j) = ((c+c_lr*xc(1))**2)*((w_n(i,j,0)-rho0) +
     &         ((mygamma-1)/(2*rho0))*((w_n(i,j,0)-rho0)**2) +
     &         (w_n(i,j,0)/c_p) *
c     &         (rho0/c_p) *
     &         (w_n(i,j,3)/w_n(i,j,0)))+p0
c            pres(i,j) = (w_n(i,j,0)*R_tilde*temp(i,j))/M +
c     &         ((w_n(i,j,0)**2)/(M**2))*(R_tilde*temp(i,j) *
c     &         0.03649-136100)+((w_n(i,j,0)**3)/(M**3))*(R_tilde *
c     &         temp(i,j)*(0.03649**2))

c            temp(i,j) = 
c     &         (temp(i,j)/c_p)*(w_n(i,j,3)/w_n(i,j,0)) +
c     &         (1/(w_n(i,j,0)*c_p))*(pres(i,j)-p0)+T0
             temp(i,j) =
     &          T0*exp((w_n(i,j,3)/w_n(i,j,0) -
     &          R*log(rho0/w_n(i,j,0)))/c_v)
          end do
       end do

      return
      end
