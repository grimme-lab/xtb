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
#define XTB_API_SUFFIX__VERSION_6_3_0

/// Define proprocessor to allow to check for specific API features
#define XTB_API_VERSION 10000
#define XTB_VERSION_6_3_0   1

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
 * Query for semantic API version
**/

/// Returns API version as 10000 * major + 100 * minor + 1 * patch
extern XTB_API_ENTRY int XTB_API_CALL
xtb_getAPIVersion() XTB_API_SUFFIX__VERSION_6_3_0;

/*
 * Calculation environment
**/

/// Create new xtb calculation environment object
extern XTB_API_ENTRY xtb_TEnvironment XTB_API_CALL
xtb_newEnvironment(void) XTB_API_SUFFIX__VERSION_6_3_0;

/// Delete a xtb calculation environment object
extern XTB_API_ENTRY void XTB_API_CALL
xtb_delEnvironment(xtb_TEnvironment* /* env */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Check current status of calculation environment
extern XTB_API_ENTRY int XTB_API_CALL
xtb_checkEnvironment(xtb_TEnvironment /* env */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Show and empty error stack
extern XTB_API_ENTRY void XTB_API_CALL
xtb_showEnvironment(xtb_TEnvironment /* env */,
                    const char* /* message */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Return and empty error stack
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getError(xtb_TEnvironment /* env */,
             char* /* buffer */,
             const int* /* buffersize */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Bind output from this environment
extern XTB_API_ENTRY void XTB_API_CALL
xtb_setOutput(xtb_TEnvironment /* env */,
              const char* /* filename */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Release output unit from this environment
extern XTB_API_ENTRY void XTB_API_CALL
xtb_releaseOutput(xtb_TEnvironment /* env */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Set verbosity of calculation output
extern XTB_API_ENTRY void XTB_API_CALL
xtb_setVerbosity(xtb_TEnvironment /* env */,
                 int /* verbosity */) XTB_API_SUFFIX__VERSION_6_3_0;

/*
 * Molecular structure data class
**/

/// Create new molecular structure data (quantities in Bohr)
extern XTB_API_ENTRY xtb_TMolecule XTB_API_CALL
xtb_newMolecule(xtb_TEnvironment /* env */,
                const int* /* natoms */,
                const int* /* numbers [natoms] */,
                const double* /* positions [natoms][3] */,
                const double* /* charge in e */,
                const int* /* uhf */,
                const double* /* lattice [3][3] */,
                const bool* /* periodic [3] */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Delete molecular structure data
extern XTB_API_ENTRY void XTB_API_CALL
xtb_delMolecule(xtb_TMolecule* /* mol */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Update coordinates and lattice parameters (quantities in Bohr)
extern XTB_API_ENTRY void XTB_API_CALL
xtb_updateMolecule(xtb_TEnvironment /* env */,
                   xtb_TMolecule /* mol */,
                   const double* /* positions [natoms][3] */,
                   const double* /* lattice [3][3] */) XTB_API_SUFFIX__VERSION_6_3_0;

/*
 * Singlepoint calculator
**/

/// Create new calculator object
extern XTB_API_ENTRY xtb_TCalculator XTB_API_CALL
xtb_newCalculator(void) XTB_API_SUFFIX__VERSION_6_3_0;

/// Delete calculator object
extern XTB_API_ENTRY void XTB_API_CALL
xtb_delCalculator(xtb_TCalculator* /* calc */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Load GFN0-xTB calculator
extern XTB_API_ENTRY void XTB_API_CALL
xtb_loadGFN0xTB(xtb_TEnvironment /* env */,
                xtb_TMolecule /* mol */,
                xtb_TCalculator /* calc */,
                char* /* filename */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Load GFN1-xTB calculator
extern XTB_API_ENTRY void XTB_API_CALL
xtb_loadGFN1xTB(xtb_TEnvironment /* env */,
                xtb_TMolecule /* mol */,
                xtb_TCalculator /* calc */,
                char* /* filename */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Load GFN2-xTB calculator
extern XTB_API_ENTRY void XTB_API_CALL
xtb_loadGFN2xTB(xtb_TEnvironment /* env */,
                xtb_TMolecule /* mol */,
                xtb_TCalculator /* calc */,
                char* /* filename */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Load GFN-FF calculator
extern XTB_API_ENTRY void XTB_API_CALL
xtb_loadGFNFF(xtb_TEnvironment /* env */,
              xtb_TMolecule /* mol */,
              xtb_TCalculator /* calc */,
              char* /* filename */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Add a solvation model to calculator (requires loaded parametrisation)
extern XTB_API_ENTRY void XTB_API_CALL
xtb_setSolvent(xtb_TEnvironment /* env */,
               xtb_TCalculator /* calc */,
               char* /* solvent */,
               int* /* state */,
               double* /* temp */,
               int* /* grid */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Unset the solvation model
extern XTB_API_ENTRY void XTB_API_CALL
xtb_releaseSolvent(xtb_TEnvironment /* env */,
                   xtb_TCalculator /* calc */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Add a external charge potential to calculator (only supported in GFN1/2-xTB)
extern XTB_API_ENTRY void XTB_API_CALL
xtb_setExternalCharges(xtb_TEnvironment /* env */,
                       xtb_TCalculator /* calc */,
                       int* /* n */,
                       int* /* numbers [n] */,
                       double* /* charges [n] */,
                       double* /* positions [n][3] */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Unset the external charge potential
extern XTB_API_ENTRY void XTB_API_CALL
xtb_releaseExternalCharges(xtb_TEnvironment /* env */,
                           xtb_TCalculator /* calc */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Set numerical accuracy of calculator in the range of 1000 to 0.0001
extern XTB_API_ENTRY void XTB_API_CALL
xtb_setAccuracy(xtb_TEnvironment /* env */,
                xtb_TCalculator /* calc */,
                double /* accuracy */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Set maximum number of iterations for self-consistent TB calculators
extern XTB_API_ENTRY void XTB_API_CALL
xtb_setMaxIter(xtb_TEnvironment /* env */,
               xtb_TCalculator /* calc */,
               int /* iterations */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Set electronic temperature for level filling in tight binding calculators in K
extern XTB_API_ENTRY void XTB_API_CALL
xtb_setElectronicTemp(xtb_TEnvironment /* env */,
                      xtb_TCalculator /* calc */,
                      double /* temperature */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Perform singlepoint calculation
extern XTB_API_ENTRY void XTB_API_CALL
xtb_singlepoint(xtb_TEnvironment /* env */,
                xtb_TMolecule /* mol */,
                xtb_TCalculator /* calc */,
                xtb_TResults /* res */) XTB_API_SUFFIX__VERSION_6_3_0;

/*
 * Calculation results
**/

/// Create new singlepoint results object
extern XTB_API_ENTRY xtb_TResults XTB_API_CALL
xtb_newResults(void) XTB_API_SUFFIX__VERSION_6_3_0;

/// Delete singlepoint results object
extern XTB_API_ENTRY void XTB_API_CALL
xtb_delResults(xtb_TResults* /* res */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Create copy from a singlepoint results object
extern XTB_API_ENTRY xtb_TResults XTB_API_CALL
xtb_copyResults(xtb_TResults /* res */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Query singlepoint results object for energy in Hartree
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getEnergy(xtb_TEnvironment /* env */,
              xtb_TResults /* res */,
              double* /* energy */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Query singlepoint results object for gradient in Hartree / Bohr
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getGradient(xtb_TEnvironment /* env */,
                xtb_TResults /* res */,
                double* /* gradient [natoms][3] */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Query singlepoint results object for virial in Hartree
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getVirial(xtb_TEnvironment /* env */,
              xtb_TResults /* res */,
              double* /* virial [3][3] */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Query singlepoint results object for dipole in e Bohr
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getDipole(xtb_TEnvironment /* env */,
              xtb_TResults /* res */,
              double* /* dipole [3] */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Query singlepoint results object for partial charges in e
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getCharges(xtb_TEnvironment /* env */,
               xtb_TResults /* res */,
               double* /* charges [natoms] */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Query singlepoint results object for bond orders
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getBondOrders(xtb_TEnvironment /* env */,
                  xtb_TResults /* res */,
                  double* /* wbo [natoms][natoms] */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Query singlepoint results object for the number of basis functions
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getNao(xtb_TEnvironment /* env */,
           xtb_TResults /* res */,
           int* /* nao */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Query singlepoint results object for orbital energies in Hartree [nao]
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getOrbitalEigenvalues(xtb_TEnvironment /* env */,
                          xtb_TResults /* res */,
                          double* /* emo */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Query singlepoint results object for occupation numbers [nao]
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getOrbitalOccupations(xtb_TEnvironment /* env */,
                          xtb_TResults /* res */,
                          double* /* focc */) XTB_API_SUFFIX__VERSION_6_3_0;

/// Query singlepoint results object for orbital coefficients [nao][nao]
extern XTB_API_ENTRY void XTB_API_CALL
xtb_getOrbitalCoefficients(xtb_TEnvironment /* env */,
                           xtb_TResults /* res */,
                           double* /* c */) XTB_API_SUFFIX__VERSION_6_3_0;

#ifdef __cplusplus
}
#endif
