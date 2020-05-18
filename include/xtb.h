/* This file is part of xtb.
 *
 * Copyright (C) 2019-2020  Sebastian Ehlert
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#define XTB_API_ENTRY
#define XTB_API_CALL
#define XTB_API_SUFFIX__VERSION_6_3

/// Define proprocessor to allow to check for specific API features
#define XTB_VERSION_6_3   1

/// Possible print levels for API calls
#define XTB_VERBOSITY_FULL    2
#define XTB_VERBOSITY_MINIMAL 1
#define XTB_VERBOSITY_MUTED   0

#ifdef __cplusplus
extern "C" {
#else
#include <stdbool.h>
#endif

/*
 * Opaque pointers to Fortran objects
**/

/// Calculation environment class
typedef struct _xtb_TEnvironment* xtb_TEnvironment;

/// Molecular structure data class
typedef struct _xtb_TMolecule* xtb_TMolecule;

/// Single point calculator class
typedef struct _xtb_TCalculator* xtb_TCalculator;

/// Single point results class
typedef struct _xtb_TResults* xtb_TResults;

/*
 * Calculation environment
**/

/// Create new xtb calculation environment object
extern XTB_API_ENTRY xtb_TEnvironment XTB_API_CALL
xtb_newEnvironment(void) XTB_API_SUFFIX__VERSION_6_3;

/// Delete a xtb calculation environment object
extern XTB_API_ENTRY void XTB_API_CALL
xtb_delEnvironment(xtb_TEnvironment* /* env */) XTB_API_SUFFIX__VERSION_6_3;

/// Check current status of calculation environment
extern XTB_API_ENTRY int XTB_API_CALL
xtb_checkEnvironment(xtb_TEnvironment /* env */) XTB_API_SUFFIX__VERSION_6_3;

/// Show and empty error stack
extern XTB_API_ENTRY void XTB_API_CALL
xtb_showEnvironment(xtb_TEnvironment /* env */,
                    const char* /* message */) XTB_API_SUFFIX__VERSION_6_3;

/// Bind output from this environment
extern XTB_API_ENTRY void XTB_API_CALL
xtb_setOutput(xtb_TEnvironment /* env */,
              const char* /* filename */) XTB_API_SUFFIX__VERSION_6_3;

/// Release output unit from this environment
extern XTB_API_ENTRY void XTB_API_CALL
xtb_releaseOutput(xtb_TEnvironment /* env */) XTB_API_SUFFIX__VERSION_6_3;

/// Set verbosity of calculation output
extern XTB_API_ENTRY void XTB_API_CALL
xtb_setVerbosity(xtb_TEnvironment /* env */,
                 int /* verbosity */) XTB_API_SUFFIX__VERSION_6_3;

/*
 * Molecular structure data class
**/

/// Create new molecular structure data
extern XTB_API_ENTRY xtb_TMolecule XTB_API_CALL
xtb_newMolecule(xtb_TEnvironment /* env */,
                const int* /* natoms */,
                const int* /* numbers [natoms] */,
                const double* /* positions [natoms][3] */,
                const double* /* charge in e */,
                const int* /* uhf */,
                const double* /* lattice [3][3] */,
                const bool* /* periodic [3] */) XTB_API_SUFFIX__VERSION_6_3;

/// Delete molecular structure data
extern XTB_API_ENTRY void XTB_API_CALL
xtb_delMolecule(xtb_TMolecule* /* mol */) XTB_API_SUFFIX__VERSION_6_3;

/// Update coordinates and lattice parameters
extern XTB_API_ENTRY void XTB_API_CALL
xtb_updateMolecule(xtb_TEnvironment /* env */,
                   xtb_TMolecule /* mol */,
                   const double* /* positions [natoms][3] */,
                   const double* /* lattice [3][3] */) XTB_API_SUFFIX__VERSION_6_3;

/*
 * Singlepoint calculator
**/

/// Create new calculator object
extern XTB_API_ENTRY xtb_TCalculator XTB_API_CALL
xtb_newCalculator(void) XTB_API_SUFFIX__VERSION_6_3;

/// Delete calculator object
extern XTB_API_ENTRY void XTB_API_CALL
xtb_delCalculator(xtb_TCalculator* /* calc */) XTB_API_SUFFIX__VERSION_6_3;

/// Load GFN0-xTB calculator
extern XTB_API_ENTRY void XTB_API_CALL
xtb_loadGFN0xTB(xtb_TEnvironment /* env */,
                xtb_TMolecule /* mol */,
                xtb_TCalculator /* calc */,
                char* /* filename */) XTB_API_SUFFIX__VERSION_6_3;

/// Load GFN1-xTB calculator
extern XTB_API_ENTRY void XTB_API_CALL
xtb_loadGFN1xTB(xtb_TEnvironment /* env */,
                xtb_TMolecule /* mol */,
                xtb_TCalculator /* calc */,
                char* /* filename */) XTB_API_SUFFIX__VERSION_6_3;

/// Load GFN2-xTB calculator
extern XTB_API_ENTRY void XTB_API_CALL
xtb_loadGFN2xTB(xtb_TEnvironment /* env */,
                xtb_TMolecule /* mol */,
                xtb_TCalculator /* calc */,
                char* /* filename */) XTB_API_SUFFIX__VERSION_6_3;

/// Load GFN-FF calculator
extern XTB_API_ENTRY void XTB_API_CALL
xtb_loadGFNFF(xtb_TEnvironment /* env */,
              xtb_TMolecule /* mol */,
              xtb_TCalculator /* calc */,
              char* /* filename */) XTB_API_SUFFIX__VERSION_6_3;

/// Add a solvation model to calculator (requires loaded parametrisation)
extern XTB_API_ENTRY void XTB_API_CALL
xtb_setSolvent(xtb_TEnvironment /* env */,
               xtb_TCalculator /* calc */,
               char* /* solvent */,
               int* /* state */,
               double* /* temp */,
               int* /* grid */) XTB_API_SUFFIX__VERSION_6_3;

/// Unset the solvation model
extern XTB_API_ENTRY void XTB_API_CALL
xtb_releaseSolvent(xtb_TEnvironment /* env */,
                   xtb_TCalculator /* calc */) XTB_API_SUFFIX__VERSION_6_3;

/// Add a external charge potential to calculator (only supported in GFN1/2-xTB)
extern XTB_API_ENTRY void XTB_API_CALL
xtb_setExternalCharges(xtb_TEnvironment /* env */,
                       xtb_TCalculator /* calc */,
                       int* /* n */,
                       int* /* numbers [n] */,
                       double* /* charges [n] */,
                       double* /* positions [n][3] */) XTB_API_SUFFIX__VERSION_6_3;

/// Unset the external charge potential
extern XTB_API_ENTRY void XTB_API_CALL
xtb_releaseExternalCharges(xtb_TEnvironment /* env */,
                           xtb_TCalculator /* calc */) XTB_API_SUFFIX__VERSION_6_3;

/// Perform singlepoint calculation
extern XTB_API_ENTRY void XTB_API_CALL
xtb_singlepoint(xtb_TEnvironment /* env */,
                xtb_TMolecule /* mol */,
                xtb_TCalculator /* calc */,
                xtb_TResults /* res */) XTB_API_SUFFIX__VERSION_6_3;

/*
 * Calculation results
**/

/// Create new singlepoint results object
extern XTB_API_ENTRY xtb_TResults XTB_API_CALL
xtb_newResults(void) XTB_API_SUFFIX__VERSION_6_3;

/// Delete singlepoint results object
extern XTB_API_ENTRY void XTB_API_CALL
xtb_delResults(xtb_TResults* /* res */) XTB_API_SUFFIX__VERSION_6_3;

/// Query singlepoint results object for energy
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getEnergy(xtb_TEnvironment /* env */,
              xtb_TResults /* res */,
              double* /* energy */) XTB_API_SUFFIX__VERSION_6_3;

/// Query singlepoint results object for gradient
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getGradient(xtb_TEnvironment /* env */,
                xtb_TResults /* res */,
                double* /* gradient [natoms][3] */) XTB_API_SUFFIX__VERSION_6_3;

/// Query singlepoint results object for virial
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getVirial(xtb_TEnvironment /* env */,
              xtb_TResults /* res */,
              double* /* virial [3][3] */) XTB_API_SUFFIX__VERSION_6_3;

/// Query singlepoint results object for dipole
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getDipole(xtb_TEnvironment /* env */,
              xtb_TResults /* res */,
              double* /* dipole [3] */) XTB_API_SUFFIX__VERSION_6_3;

/// Query singlepoint results object for partial charges
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getCharges(xtb_TEnvironment /* env */,
               xtb_TResults /* res */,
               double* /* charges [natoms] */) XTB_API_SUFFIX__VERSION_6_3;

/// Query singlepoint results object for bond orders
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getBondOrders(xtb_TEnvironment /* env */,
                  xtb_TResults /* res */,
                  double* /* wbo [natoms][natoms] */) XTB_API_SUFFIX__VERSION_6_3;

#ifdef __cplusplus
}
#endif
