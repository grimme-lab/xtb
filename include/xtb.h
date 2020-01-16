/* extended tight binding program package
 * Copyright (C) 2019-2020  Stefan Grimme (xtb@thch.uni-bonn.de)
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

typedef struct PEEQ_options {
   int prlevel;
   int parallel;
   double acc;
   double etemp;
   bool grad;
   bool ccm;
   char solvent[20];
} PEEQ_options;

typedef struct SCC_options {
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

/// opaque molecular structure type (allocated, manipulated and freed by xtb)
typedef void* xTB_molecule;

/// opaque calculation parameters (allocated and freed by xtb)
typedef void* xTB_parameters; // dummy

/// opaque tight binding basisset (allocated and freed by xtb)
typedef void* xTB_basisset;

/// opaque tight binding wavefunction (allocated and freed by xtb)
typedef void* xTB_wavefunction;

/// opaque external potential data (allocated, manipulated and freed by xtb)
typedef void* xTB_external; // dummy

#ifdef __cplusplus
extern "C" {
#endif

/// Constructor for the molecular structure type, returns a nullptr in case
/// the generation fails. All values have to be provided.
extern xTB_molecule
new_xTB_molecule(
      const int* natoms,
      const int* numbers, // [natoms]
      const double* positions, // [3*natoms3]
      const double* charge,
      const int* uhf,
      const double* lattice, // [3*3]
      const bool* periodic); // [3]

/// Updates the molecular structure type. This routine checks for the allocation
/// state of the opaque pointer before copying positions and lattice.
extern int
update_xTB_molecule(
      xTB_molecule mol,
      const double* positions, // [3*natoms]
      const double* lattice); // [3*3]

/// Deconstructor for the molecular structure type.
extern void
delete_xTB_molecule(xTB_molecule mol);

/// Constructor for the tight binding wavefunction.
extern xTB_wavefunction
new_xTB_wavefunction(const xTB_molecule mol, const xTB_basisset basis);

/// Return dimensions of this tight binding wavefunction.
extern void
dimensions_xTB_wavefunction(
      const xTB_wavefunction wfn,
      int* nat,
      int* nao,
      int* nsh);

/// Query the wavefunction properties, pass NULL/nullptr in case you want
/// a query to be ignored.
/// In case the wavefunction is not allocated nothing is returned.
extern void
query_xTB_wavefunction(
      const xTB_wavefunction wfn,
      double* charges, // [nat]
      double* dipoles, // [3*nat]
      double* quadrupoles, // [6*nat]
      double* bond_orders, // [nat*nat]
      double* hl_gap,
      double* orbital_energies); // [nao]

/// Deconstructor for the tight binding basisset.
extern void
delete_xTB_wavefunction(xTB_wavefunction wfn);

/// Constructor for the tight binding basisset.
extern xTB_basisset
new_xTB_basisset(const xTB_basisset basis);

/// Return dimensions of this tight binding basisset.
/// For not allocated input zero is returned, every argument is optional.
extern void
dimensions_xTB_basisset(
      const xTB_basisset basis,
      int* nat,
      int* nbf,
      int* nao,
      int* nsh);

/// Deconstructor for the tight binding wavefunction.
extern void
delete_xTB_basisset(xTB_basisset basis);

/// Load GFN-xTB parametrisation and create a lock, the file name is optional.
/// The lock is maintained by the shared libary. This call is secured with
/// an omp critical mutex called xtb_load_api.
extern int
load_xTB_parameters(const int* gfn, const char* filename);

/// Check if xTB was correctly initialized.
extern bool
check_xTB_init(void);

/// Check if xTB was locked on the correct parametrisation.
extern bool
check_xTB_lock(const int* gfn);

/// Release the lock on the xTB parametrisation.
extern void
reset_xTB_lock(void);

/// Unsafe xTB calculation, this calculation mode requires to setup at least
/// the parametrisation and the molecular structure data.
///
/// In case you pass null-pointers instead of the basisset and wavefunction,
/// those will be constructed on-the-fly and deleted afterwards.
/// Properties like energy, gradient and stress tensor are returned directly,
/// while other properties can be obtained from querying the wavefunction.
extern int
xTB_calculation(
      const xTB_molecule mol,
      const xTB_parameters param, // not used
      const xTB_basisset basis, // optional
      const xTB_wavefunction wfn, // optional
      const xTB_external pcem, // not used
      const SCC_options* opt,
      const char* output,
      double* energy,
      double* gradient, // [3*nat]
      double* stress); // [3*3]

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
      const double* tmp);

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
