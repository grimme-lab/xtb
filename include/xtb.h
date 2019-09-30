/* extended tight binding program package
 * Copyright (C) 2019  Stefan Grimme (xtb@thch.uni-bonn.de)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#ifdef __cplusplus
namespace xtb {
#else
#include <stdbool.h>
#endif

typedef struct _PEEQ_options {
   int prlevel;
   int parallel;
   double acc;
   double etemp;
   bool grad;
   bool ccm;
} PEEQ_options;

typedef struct _SCC_options {
   int prlevel;
   int parallel;
   double acc;
   double etemp;
   bool grad;
   bool restart;
   int maxiter;
   char solvent[20];
} SCC_options;

#ifdef __cplusplus
extern "C" {
#endif

extern int
GFN0_PBC_calculation(const int* natoms, const int* attyp, const double* charge,
      const double* coord, const double* lattice, const bool* pbc,
      const PEEQ_options* opt, const char* output,
      double* energy, double* grad, double* glat);

extern int
GFN2_calculation(const int* natoms, const int* attyp, const double* charge,
      const double* coord, const SCC_options* opt, const char* output,
      double* energy, double* grad, double* dipole, double* q,
      double* dipm, double* qp, double* wbo);

extern int
GFN2_QMMM_calculation(const int* natoms, const int* attyp, const double* charge,
      const double* coord, const SCC_options* opt, const char* output,
      const int* npc, const double* pc_q, const int* pc_at, const double* pc_gam,
      const double* pc_coord, double* energy, double* grad, double* pc_grad);

extern int
GFN1_calculation(const int* natoms, const int* attyp, const double* charge,
      const double* coord, const SCC_options* opt, const char* output,
      double* energy, double* grad);

extern int
GFN1_QMMM_calculation(const int* natoms, const int* attyp, const double* charge,
      const double* coord, const SCC_options* opt, const char* output,
      const int* npc, const double* pc_q, const int* pc_at, const double* pc_gam,
      const double* pc_coord, double* energy, double* grad, double* pc_grad);

extern int
GFN0_calculation(const int* natoms, const int* attyp, const double* charge,
      const double* coord, const PEEQ_options* opt, const char* output,
      double* energy, double* grad);

#ifdef __cplusplus
}
}
#endif
