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
   char solvent[20];
} PEEQ_options;

typedef struct _SCC_options {
   int prlevel;
   int parallel;
   double acc;
   double etemp;
   bool grad;
   bool restart;
   bool ccm;
   int maxiter;
   char solvent[20];
} SCC_options;

#ifdef __cplusplus
extern "C" {
#endif

extern int
GFN0_PBC_calculation(
      const int* natoms,
      const int* attyp,
      const double* charge,
      const int* uhf,
      const double* coord,
      const double* lattice,
      const bool* pbc,
      const PEEQ_options* opt,
      const char* output,
      double* energy,
      double* grad,
      double* stress,
      double* glat);

extern int
GFN2_calculation(
      const int* natoms,
      const int* attyp,
      const double* charge,
      const int* uhf,
      const double* coord,
      const SCC_options* opt,
      const char* output,
      double* energy,
      double* grad,
      double* dipole,
      double* q,
      double* dipm,
      double* qp,
      double* wbo);

extern int
GFN2_QMMM_calculation(
      const int* natoms,
      const int* attyp,
      const double* charge,
      const int* uhf,
      const double* coord,
      const SCC_options* opt,
      const char* output,
      const int* npc,
      const double* pc_q,
      const int* pc_at,
      const double* pc_gam,
      const double* pc_coord,
      double* energy,
      double* grad,
      double* pc_grad);

extern int
GFN1_calculation(
      const int* natoms,
      const int* attyp,
      const double* charge,
      const int* uhf,
      const double* coord,
      const SCC_options* opt,
      const char* output,
      double* energy,
      double* grad,
      double* dipole,
      double* q,
      double* wbo);

extern int
GFN1_QMMM_calculation(
      const int* natoms,
      const int* attyp,
      const double* charge,
      const int* uhf,
      const double* coord,
      const SCC_options* opt,
      const char* output,
      const int* npc,
      const double* pc_q,
      const int* pc_at,
      const double* pc_gam,
      const double* pc_coord,
      double* energy,
      double* grad,
      double* pc_grad);

extern int
GFN0_calculation(
      const int* natoms,
      const int* attyp,
      const double* charge,
      const int* uhf,
      const double* coord,
      const PEEQ_options* opt,
      const char* output,
      double* energy,
      double* grad);

extern int
GBSA_model_preload(
      const double* epsv,
      const double* smass,
      const double* rhos,
      const double* c1,
      const double* rprobe,
      const double* gshift,
      const double* soset,
      const double* dum,
      const double* gamscale,
      const double* sx,
      const double tmp);

extern int
GBSA_calculation(
      const int* natoms,
      const int* attyp,
      const double* coord,
      const char* solvent,
      const int* reference,
      const double* temperature,
      const int* method,
      const int* grid_size,
      const char* file,
      double* brad,
      double* sasa);

#ifdef __cplusplus
}
}
#endif
