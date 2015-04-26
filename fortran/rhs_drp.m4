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

c This function calculates the RHS using a DRP scheme.

      subroutine rhs(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  xylo,
     &  dxy,
     &  w_n,
     &  H0,
     &  KK,
     &  pres,
     &  temp,
     &  alpha,
     &  kappa,
     &  mu,
     &  muB,
     &  R,
     &  T_star_n,
     &  T_star_o,
     &  tau_n,
     &  tau_o,
     &  tolerance,
     &  level,
     &  c,
     &  c_lr,
     &  nequ)
c***********************************************************************
      implicit none
cinclude(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      REAL alpha
      REAL kappa
      REAL mu
      REAL muB
      REAL R
      REAL T_star_n
      REAL T_star_o
      REAL tau_n
      REAL tau_o
      REAL tolerance(0:2)
      integer level
      REAL c
      REAL c_lr
      integer nequ

      REAL dxy(0:1)
      REAL xylo(0:1)
c
c variables in 2d cell indexed
      REAL
     &     w_n(CELL2dVECG(ifirst,ilast,gcw),0:nequ-1),
     &     H0(CELL2d(ifirst,ilast,0)),
     &     KK(CELL2d(ifirst,ilast,0),0:nequ-1),
     &     pres(CELL2dVECG(ifirst,ilast,gcw)),
     &     temp(CELL2dVECG(ifirst,ilast,gcw))
c
c***********************************************************************
c***********************************************************************     
c
      integer k,q
      integer ix,i,j
      integer ineq,ii,jj
      integer tag
      REAL t_prime_x,t_prime_y,t_dprime_x,t_dprime_y
      REAL u_prime_x,u_prime_y,v_prime_y,v_prime_x
      REAL u_dprime_x,u_dprime_y,v_dprime_xy
      REAL u_dprime_yx,v_dprime_x,v_dprime_y
      REAL phi_xx,phi_xy_yx,phi_yy
      REAL phi_xx_prime_x,phi_xy_prime_y
      REAL phi_yx_prime_x,phi_yy_prime_y
      REAL F(0:7-1,0:nequ-1),G(0:7-1,0:nequ-1),H(0:nequ-1)
      REAL Mn,Mo,c_vn,c_vo,A_n,A_o
      REAL T_n,T_o
      REAL sigma
      REAL a_1,a_2,a_3
      REAL tmp1,tmp2

      a_1=0.79926643
      a_2=-0.18941314
      a_3=0.02651995

      do j=ifirst1,ilast1
         do i=ifirst0,ilast0
            ix = 3
            do q=-3,3
               F(ix+q,0) = w_n(i+q,j,0) *
     &            (w_n(i+q,j,1)/w_n(i+q,j,0))
               G(ix+q,0) = w_n(i,j+q,0) *
     &            (w_n(i,j+q,2)/w_n(i,j+q,0))

               F(ix+q,1) = w_n(i+q,j,0) *
     &            ((w_n(i+q,j,1)/w_n(i+q,j,0))**2) +
     &            pres(i+q,j)
               G(ix+q,1) = w_n(i,j+q,0) *
     &            (w_n(i,j+q,2)/w_n(i,j+q,0)) *
     &            (w_n(i,j+q,1)/w_n(i,j+q,0))

               F(ix+q,2) = w_n(i+q,j,0) *
     &            (w_n(i+q,j,1)/w_n(i+q,j,0)) *
     &            (w_n(i+q,j,2)/w_n(i+q,j,0))
               G(ix+q,2) = w_n(i,j+q,0) *
     &            ((w_n(i,j+q,2)/w_n(i,j+q,0))**2) +
     &            pres(i,j+q)

               F(ix+q,3) = w_n(i+q,j,0) *
     &            (w_n(i+q,j,1)/w_n(i+q,j,0)) *
     &            (w_n(i+q,j,3)/w_n(i+q,j,0))
               G(ix+q,3) = w_n(i,j+q,0) *
     &            (w_n(i,j+q,2)/w_n(i,j+q,0)) *
     &            (w_n(i,j+q,3)/w_n(i,j+q,0))

               F(ix+q,4) = w_n(i+q,j,0) *
     &            (w_n(i+q,j,1)/w_n(i+q,j,0)) *
     &            (w_n(i+q,j,4)/w_n(i+q,j,0))
               G(ix+q,4) = w_n(i,j+q,0) *
     &            (w_n(i,j+q,2)/w_n(i,j+q,0)) *
     &            (w_n(i,j+q,4)/w_n(i,j+q,0))
             end do

            H(0)=H0(i,j)

            u_dprime_x=((w_n(i-1,j,1)/w_n(i-1,j,0)) -
     &         2*(w_n(i,j,1)/w_n(i,j,0)) +
     &         (w_n(i+1,j,1)/w_n(i+1,j,0))) /
     &         (dxy(0)**2)
            u_dprime_y=((w_n(i,j-1,1)/w_n(i,j-1,0)) -
     &         2*(w_n(i,j,1)/w_n(i,j,0)) +
     &         (w_n(i,j+1,1)/w_n(i,j+1,0))) /
     &         (dxy(1)**2)
            v_dprime_xy=((w_n(i+1,j+1,2)/w_n(i+1,j+1,0)) -
     &         (w_n(i+1,j-1,2)/w_n(i+1,j-1,0)) -
     &         (w_n(i-1,j+1,2)/w_n(i-1,j+1,0)) +
     &         (w_n(i-1,j-1,2)/w_n(i-1,j-1,0))) /
     &         (4*dxy(0)*dxy(1))
            phi_xx_prime_x=(4.0/3.0)*u_dprime_x -
     &         (2.0/3.0)*v_dprime_xy
            phi_xy_prime_y=u_dprime_y+v_dprime_xy
            H(1)=muB*(u_dprime_x+v_dprime_xy) +
     &         mu*(phi_xx_prime_x+phi_xy_prime_y)
c            H(1)=0

            v_dprime_y=((w_n(i,j-1,2)/w_n(i,j-1,0)) -
     &         2*(w_n(i,j,2)/w_n(i,j,0)) +
     &         (w_n(i,j+1,2)/w_n(i,j+1,0))) /
     &         (dxy(1)**2)
            v_dprime_x=((w_n(i-1,j,2)/w_n(i-1,j,0)) -
     &         2*(w_n(i,j,2)/w_n(i,j,0)) +
     &         (w_n(i+1,j,2)/w_n(i+1,j,0))) /
     &         (dxy(0)**2)
            u_dprime_yx=((w_n(i+1,j+1,1)/w_n(i+1,j+1,0)) -
     &         (w_n(i+1,j-1,1)/w_n(i+1,j-1,0)) -
     &         (w_n(i-1,j+1,1)/w_n(i-1,j+1,0)) +
     &         (w_n(i-1,j-1,1)/w_n(i-1,j-1,0))) /
     &         (4*dxy(0)*dxy(1))
            phi_yx_prime_x=u_dprime_yx+v_dprime_x
            phi_yy_prime_y=(4.0/3.0)*v_dprime_y-(2.0/3.0)*u_dprime_yx
            H(2)=muB*(v_dprime_y+u_dprime_yx) +
     &         mu*(phi_yx_prime_x+phi_yy_prime_y)
c            H(2)=0

            u_prime_x=((w_n(i+1,j,1)/w_n(i+1,j,0)) -
     &         (w_n(i-1,j,1)/w_n(i-1,j,0)))/(2*dxy(0))
            u_prime_y=((w_n(i,j+1,1)/w_n(i,j+1,0)) -
     &         (w_n(i,j-1,1)/w_n(i,j-1,0)))/(2*dxy(1))
            v_prime_y=((w_n(i,j+1,2)/w_n(i,j+1,0)) -
     &         (w_n(i,j-1,2)/w_n(i,j-1,0)))/(2*dxy(1))
            v_prime_x=((w_n(i+1,j,2)/w_n(i+1,j,0)) -
     &         (w_n(i-1,j,2)/w_n(i-1,j,0)))/(2*dxy(0))
            t_prime_x=(temp(i+1,j)-temp(i-1,j))/(2*dxy(0))
            t_prime_y=(temp(i,j+1)-temp(i,j-1))/(2*dxy(1))
            t_dprime_x=(temp(i-1,j)-2*temp(i,j) +
     &         temp(i+1,j))/(dxy(0)**2)
            t_dprime_y=(temp(i,j-1)-2*temp(i,j) +
     &         temp(i,j+1))/(dxy(1)**2)
            T_n=w_n(i,j,4)/w_n(i,j,0)
            Mn=(temp(i,j)-T_n)/tau_n
            tmp1=T_star_n/(T_n)
            c_vn=0.78*R*(tmp1**2)*exp(-tmp1)
            A_n=((temp(i,j)/T_n)-1)*c_vn
            phi_xx=(4.0/3.0)*u_prime_x-(2.0/3.0)*v_prime_y
            phi_xy_yx=u_prime_y+v_prime_x
            phi_yy=(4.0/3.0)*v_prime_y-(2.0/3.0)*u_prime_x
            sigma=(muB/temp(i,j))*((u_prime_x+v_prime_y)**2) +
     &         (mu/(2*temp(i,j))) *
     &         (phi_xx**2+2*(phi_xy_yx**2)+phi_yy**2) +
     &         (kappa/(temp(i,j)**2)) *
     &         (t_prime_x**2+t_prime_y**2) +
     &         (w_n(i,j,0)/temp(i,j))
     &         *(A_n*Mn)
            H(3)=sigma-
     &         (w_n(i,j,0)/T_n)*c_vn*Mn +
     &         (-((kappa*(t_prime_x**2))/(temp(i,j)**2)) +
     &         ((kappa*t_dprime_x)/temp(i,j)) -
     &         ((kappa*(t_prime_y**2))/(temp(i,j)**2)) +
     &         ((kappa*t_dprime_y)/temp(i,j)))
c            H(3)=0

            H(4)=(w_n(i,j,0)/tau_n)*(temp(i,j)-T_n)

            ix=3
            do k=0,nequ-1
               KK(i,j,k)=
     &            (-a_3*F(ix-3,k)-a_2*F(ix-2,k)-a_1*F(ix-1,k) +
     &            a_1*F(ix+1,k)+a_2*F(ix+2,k)+a_3*F(ix+3,k)) /
     &            (-dxy(0)) +
     &            (-a_3*G(ix-3,k)-a_2*G(ix-2,k)-a_1*G(ix-1,k) +
     &            a_1*G(ix+1,k)+a_2*G(ix+2,k)+a_3*G(ix+3,k)) /
     &            (-dxy(1))+H(k)
            end do
         end do
      end do

      return
      end
