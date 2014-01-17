/*
 ===========================================================================
 Copyright (C) 2012 Jon Rood.
 
 This file is part of Enlightning source code.
 
 Enlightning source code is free software; you can redistribute it
 and/or modify it under the terms of the GNU General Public License as
 published by the Free Software Foundation; either version 3 of the License,
 or (at your option) any later version.
 
 Enlightning source code is distributed in the hope that it will be
 useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Enlightning; if not, see <http://www.gnu.org/licenses/>.
 ===========================================================================
 */

#define SQUARE(x) ((x)*(x))  //faster way to square values that using pow()
#define KSIZE 4              //number of equations in set of equations

float wenoM(float fgm1, float fg, float fgp1);
float wenoP(float fgm1, float fg, float fgp1);

float wenoM(float fgm1, float fg, float fgp1)
{
    float alpha0 = (1.0/3.0)/SQUARE(10e-6+SQUARE(fg-fgm1));
    float alpha1 = (2.0/3.0)/SQUARE(10e-6+SQUARE(fgp1-fg));
    
    float sum = alpha0+alpha1;
    
    return ((alpha0/sum)*(0.5*(-fgm1+3*fg))+(alpha1/sum)*(0.5*(fg+fgp1)));
}

float wenoP(float fgm1, float fg, float fgp1)
{
    float alpha0 = (2.0/3.0)/SQUARE(10e-6+SQUARE(fg-fgm1));
    float alpha1 = (1.0/3.0)/SQUARE(10e-6+SQUARE(fgp1-fg));
    
    float sum = alpha0+alpha1;
    
    return ((alpha0/sum)*(0.5*(fg+fgm1))+(alpha1/sum)*(0.5*(-fgp1+3*fg)));
}

__kernel void rhs(__global float *w_n,
                  __global float *K,
                  __global float *p,
                  __global float *H0,
                  const int Y,
                  const int X,
                  const float dy,
                  const float dx,
                  const float alpha)
{
    int i = get_global_id(0);
    int j = get_global_id(1);
    
    if((i < 2) || (i >= X-2))
        return;
    if(( j<2 ) || (j >= Y-2))
        return;
    
    float F[KSIZE][5], G[KSIZE][5], H[KSIZE];
    int ix = 2;
    float p_prime_x = (p[(j*X)+i+1]-p[(j*X)+i-1])/(2*dx);
    float p_prime_y = (p[((j+1)*X)+i]-p[((j-1)*X)+i])/(2*dy);
    
    H[0] = H0[(j*X)+i];
    H[1] = -p_prime_x;
    H[2] = -p_prime_y;
    H[3] = 0;
    
    for (int q = -2; q <= 2; q++) {
        F[0][ix+q] = w_n[(0*X*Y)+(j*X)+i+q]*(w_n[(1*X*Y)+(j*X)+i+q]
            /w_n[(0*X*Y)+(j*X)+i+q]);
        F[1][ix+q] = w_n[(0*X*Y)+(j*X)+i+q]*SQUARE(w_n[(1*X*Y)+(j*X)+i+q]
            /w_n[(0*X*Y)+(j*X)+i+q]);
        F[2][ix+q] = w_n[(0*X*Y)+(j*X)+i+q]*(w_n[(1*X*Y)+(j*X)+i+q]
            /w_n[(0*X*Y)+(j*X)+i+q])*(w_n[(2*X*Y)+(j*X)+i+q]/w_n[(0*X*Y)+(j*X)+i+q]);
        F[3][ix+q] = w_n[(0*X*Y)+(j*X)+i+q]*(w_n[(1*X*Y)+(j*X)+i+q]
            /w_n[(0*X*Y)+(j*X)+i+q])*(w_n[(3*X*Y)+(j*X)+i+q]/w_n[(0*X*Y)+(j*X)+i+q]);
        
        G[0][ix+q] = w_n[(0*X*Y)+((j+q)*X)+i]*(w_n[(2*X*Y)+((j+q)*X)+i]
            /w_n[(0*X*Y)+((j+q)*X)+i]);
        G[1][ix+q] = w_n[(0*X*Y)+((j+q)*X)+i]*(w_n[(2*X*Y)+((j+q)*X)+i]
            /w_n[(0*X*Y)+((j+q)*X)+i])*(w_n[(1*X*Y)+((j+q)*X)+i]/w_n[(0*X*Y)+((j+q)*X)+i]);
        G[2][ix+q] = w_n[(0*X*Y)+((j+q)*X)+i]*SQUARE(w_n[(2*X*Y)+((j+q)*X)+i]
            /w_n[(0*X*Y)+((j+q)*X)+i]);
        G[3][ix+q] = w_n[(0*X*Y)+((j+q)*X)+i]*(w_n[(2*X*Y)+((j+q)*X)+i]
            /w_n[(0*X*Y)+((j+q)*X)+i])*(w_n[(3*X*Y)+((j+q)*X)+i]/w_n[(0*X*Y)+((j+q)*X)+i]);
    }
    
    for (int k = 0; k < KSIZE; k++) {
        int ii, jj;
        float fgm1, fg, fgp1, vm1f, vm2f, vp1f, vp2f, vm1g, vm2g, vp1g, vp2g;
        
        //x
        ii = i-1;
        ix = 2-1;
        fgm1 = 0.5*(F[k][ix-1]+alpha*w_n[(k*X*Y)+(j*X)+ii-1]);
        fg   = 0.5*(F[k][ix  ]+alpha*w_n[(k*X*Y)+(j*X)+ii  ]);
        fgp1 = 0.5*(F[k][ix+1]+alpha*w_n[(k*X*Y)+(j*X)+ii+1]);
        
        vm1f=wenoM(fgm1, fg, fgp1);
        
        ii++;
        ix++;
        fgm1 = 0.5*(F[k][ix-1]-alpha*w_n[(k*X*Y)+(j*X)+ii-1]);
        fg   = 0.5*(F[k][ix  ]-alpha*w_n[(k*X*Y)+(j*X)+ii  ]);
        fgp1 = 0.5*(F[k][ix+1]-alpha*w_n[(k*X*Y)+(j*X)+ii+1]);
        
        vp1f = wenoP(fgm1, fg, fgp1);
        
        fgm1 = 0.5*(F[k][ix-1]+alpha*w_n[(k*X*Y)+(j*X)+ii-1]);
        fg   = 0.5*(F[k][ix  ]+alpha*w_n[(k*X*Y)+(j*X)+ii  ]);
        fgp1 = 0.5*(F[k][ix+1]+alpha*w_n[(k*X*Y)+(j*X)+ii+1]);
        
        vm2f = wenoM(fgm1, fg, fgp1);
        
        ii++;
        ix++;
        fgm1 = 0.5*(F[k][ix-1]-alpha*w_n[(k*X*Y)+(j*X)+ii-1]);
        fg   = 0.5*(F[k][ix  ]-alpha*w_n[(k*X*Y)+(j*X)+ii  ]);
        fgp1 = 0.5*(F[k][ix+1]-alpha*w_n[(k*X*Y)+(j*X)+ii+1]);
        
        vp2f = wenoP(fgm1, fg, fgp1);
        
        //y
        jj = j-1;
        ix = 2-1;
        fgm1 = 0.5*(G[k][ix-1]+alpha*w_n[(k*X*Y)+((jj-1)*X)+i]);
        fg   = 0.5*(G[k][ix  ]+alpha*w_n[(k*X*Y)+((jj  )*X)+i]);
        fgp1 = 0.5*(G[k][ix+1]+alpha*w_n[(k*X*Y)+((jj+1)*X)+i]);
        
        vm1g = wenoM(fgm1, fg, fgp1);
        
        jj++;
        ix++;
        fgm1 = 0.5*(G[k][ix-1]-alpha*w_n[(k*X*Y)+((jj-1)*X)+i]);
        fg   = 0.5*(G[k][ix  ]-alpha*w_n[(k*X*Y)+((jj  )*X)+i]);
        fgp1 = 0.5*(G[k][ix+1]-alpha*w_n[(k*X*Y)+((jj+1)*X)+i]);
        
        vp1g = wenoP(fgm1, fg, fgp1);
        
        fgm1 = 0.5*(G[k][ix-1]+alpha*w_n[(k*X*Y)+((jj-1)*X)+i]);
        fg   = 0.5*(G[k][ix  ]+alpha*w_n[(k*X*Y)+((jj  )*X)+i]);
        fgp1 = 0.5*(G[k][ix+1]+alpha*w_n[(k*X*Y)+((jj+1)*X)+i]);
        
        vm2g = wenoM(fgm1, fg, fgp1);
        
        jj++;
        ix++;
        fgm1 = 0.5*(G[k][ix-1]-alpha*w_n[(k*X*Y)+((jj-1)*X)+i]);
        fg   = 0.5*(G[k][ix  ]-alpha*w_n[(k*X*Y)+((jj  )*X)+i]);
        fgp1 = 0.5*(G[k][ix+1]-alpha*w_n[(k*X*Y)+((jj+1)*X)+i]);
        
        vp2g = wenoP(fgm1, fg, fgp1);
        
        K[(k*X*Y)+(j*X)+i]=(vm1f+vp1f-vm2f-vp2f)/dx+(vm1g+vp1g-vm2g-vp2g)/dy+H[k];
    }
}
