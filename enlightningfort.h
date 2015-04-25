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

extern "C" {
void SAMRAI_F77_FUNC(rhs, RHS) (const int&,                       // ifirst(0)
                                const int&,                       // ilast(0)
                                const int&,                       // ifirst(1)
                                const int&,                       // ilast(1)
                                const int&,                       // d_nghosts(0)
                                const int&,                       // d_nghosts(1)
                                const double *,                   // xlo
                                const double *,                   // dx
                                const double *,                   // w_n
                                const double *,                   // source_grid
                                double *,                         // rhs
                                const double *,                   // pressure_grid
                                const double *,                   // temperature_grid
                                const double&,                    // weno_alpha
                                const double&,                    // kappa
                                const double&,                    // mu
                                const double&,                    // muB
                                const double&,                    // R
                                const double&,                    // T_star_n
                                const double&,                    // T_star_o
                                const double&,                    // tau_n
                                const double&,                    // tau_o
                                const double *,                   // tag_tolerance
                                const int&,                       // level
                                const double&,                    // c
                                const double&,                    // c_lr
                                const int&);                      // NEQU

void SAMRAI_F77_FUNC(runge_kutta, RUNGE_KUTTA) (const int&,       // ifirst(0)
                                                const int&,       // ilast(0)
                                                const int&,       // ifirst(1)
                                                const int&,       // ilast(1)
                                                const int&,       // d_nghosts(0)
                                                const int&,       // d_nghosts(1)
                                                const double&,    // dt
                                                const double&,    // alpha_1
                                                const double&,    // alpha_2
                                                const double&,    // beta
                                                double *,         // w_n
                                                const double *,   // w
                                                const double *,   // rhs
                                                const int&);      // NEQU

void SAMRAI_F77_FUNC(update_state, UPDATE_STATE) (const double *, // dx
                                                  const double *, // xlo
                                                  const int&,     // ifirst(0)
                                                  const int&,     // ilast(0)
                                                  const int&,     // ifirst(1)
                                                  const int&,     // ilast(1)
                                                  const int&,     // d_nghosts(0)
                                                  const int&,     // d_nghosts(1)
                                                  const double *, // w_n
                                                  double *,       // pressure_grid
                                                  double *,       // temperature_grid
                                                  const double&,  // rho0
                                                  const double&,  // c
                                                  const double&,  // gamma1
                                                  const double&,  // c_p
                                                  const double&,  // p0
                                                  const double&,  // T0
                                                  const double&,  // c_lr
                                                  const double&,  // c_v
                                                  const double&,  // R_tilde
                                                  const double&,  // M
                                                  const double&,  // R
                                                  const int&);    // NEQU

void SAMRAI_F77_FUNC(tag_cells, TAG_CELLS) (const int&,           // ifirst(0)
                                            const int&,           // ilast(0)
                                            const int&,           // ifirst(1)
                                            const int&,           // ilast(1)
                                            const int&,           // d_nghosts(0)
                                            const int&,           // d_nghosts(1)
                                            int *,                // tags
                                            const double *,       // w_n
                                            const double *,       // source_grid
                                            const int&,           // true
                                            const double *,       // tag_tolerance
                                            const int&);          // NEQU
}
