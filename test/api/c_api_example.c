#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "xtb.h"

#define check(x, ...)                                                          \
  _Generic((x), int : check_int, double : check_double)(x, __VA_ARGS__)

static inline bool check_int(int actual, int expected, const char *msg) {
  if (expected == actual) {
    return true;
  }
  fprintf(stderr, "FAIL: %s: expected %d, got %d\n", msg, expected, actual);
  return false;
}

static inline bool check_double(double actual, double expected, double tol,
                                const char *msg) {
  if (fabs(expected - actual) < tol) {
    return true;
  }
  fprintf(stderr, "FAIL: %s: expected %.10f, got %.10f\n", msg, expected, actual);
  return false;
}

int testFirst() {
  // first molecule, various solvent models
  double* q;
  char* buffer;
  double* wbo;
  double* hess;

  int buffersize = 512;
  int tester = 0;
  int const natoms = 7;
  int const natsq = 49;
  int const attyp[7] = {6, 6, 6, 1, 1, 1, 1};
  double const coord[3 * 7] = {
      +0.00000000000000, +0.00000000000000, -1.79755622305860,
      +0.00000000000000, +0.00000000000000, +0.95338756106749,
      +0.00000000000000, +0.00000000000000, +3.22281255790261,
      -0.96412815539807, -1.66991895015711, -2.53624948351102,
      -0.96412815539807, +1.66991895015711, -2.53624948351102,
      +1.92825631079613, +0.00000000000000, -2.53624948351102,
      +0.00000000000000, +0.00000000000000, +5.23010455462158};

  xtb_TEnvironment env = NULL;
  xtb_TMolecule mol = NULL;
  xtb_TCalculator calc = NULL;
  xtb_TResults res = NULL;
  double energy;
  double solv_energy;
  double dipole[3];
  q = (double*) malloc(natoms * sizeof(double));
  wbo = (double*) malloc(natsq * sizeof(double));
  buffer = (char*) malloc(buffersize *sizeof(char));
  hess = (double*) malloc(9*natsq * sizeof(double));
  char solvent[] = "h2o";
  char gbsa[] = "gbsa";
  char alpb[] = "alpb";
  char cosmo[] = "cosmo";
  char cpcmx[] = "cpcmx";

  if (!check(XTB_API_VERSION, xtb_getAPIVersion(), "API version does not match"))
     goto error;

  env = xtb_newEnvironment();
  calc = xtb_newCalculator();
  res = xtb_newResults();
  mol = xtb_newMolecule(env, &natoms, attyp, coord, NULL, NULL, NULL, NULL);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_getEnergy(env, res, &energy);
  if (xtb_checkEnvironment(env)) {
    xtb_getError(env, buffer, &buffersize);
    printf("Error message is:\n%s\n", buffer);
  }

  xtb_setVerbosity(env, XTB_VERBOSITY_FULL);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_loadGFN2xTB(env, mol, calc, NULL);
  xtb_setAccuracy(env, calc, 1.0);
  xtb_setElectronicTemp(env, calc, 300.0);
  xtb_setMaxIter(env, calc, 30);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_singlepoint(env, mol, calc, res);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_getEnergy(env, res, &energy);
  xtb_getCharges(env, res, q);
  xtb_getDipole(env, res, dipole);
  xtb_getBondOrders(env, res, wbo);
  if (xtb_checkEnvironment(env))
    goto error;

  if (!check(energy, -8.3824793849585, 1.0e-9, "Energy does not match"))
    goto error;
  if (!check(q[5], 0.05184019996829, 1.0e-8, "Charge does not match"))
    goto error;
  if (!check(dipole[2], -0.298279305689518, 1.0e-6, "Dipole does not match"))
    goto error;
  if (!check(wbo[9], 2.89823984265213, 1.0e-8, "Bond order does not match"))
    goto error;

  // GBSA
  xtb_setSolvent(env, calc, solvent, NULL, NULL, NULL, gbsa);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_singlepoint(env, mol, calc, res);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_getEnergy(env, res, &energy);
  xtb_getCharges(env, res, q);
  xtb_getDipole(env, res, dipole);
  xtb_getBondOrders(env, res, wbo);
  if (xtb_checkEnvironment(env))
    goto error;

  if (!check(energy, -8.38393864716134, 1.0e-9, "GBSA Energy does not match"))
    goto error;
  if (!check(q[5], 0.06090868805034, 1.0e-8, "GBSA Charge does not match"))
    goto error;
  if (!check(dipole[2], -0.35455233974705, 1.0e-6, "GBSA Dipole does not match"))
    goto error;
  if (!check(wbo[9], +2.89453979224265, 1.0e-8, "GBSA Bond order does not match"))
    goto error;

  // ALPB
  xtb_releaseSolvent(env, calc);
  xtb_setSolvent(env, calc, solvent, NULL, NULL, NULL, alpb);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_singlepoint(env, mol, calc, res);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_getEnergy(env, res, &energy);
  xtb_getCharges(env, res, q);
  xtb_getDipole(env, res, dipole);
  xtb_getBondOrders(env, res, wbo);
  if (xtb_checkEnvironment(env))
    goto error;

  if (!check(energy, -8.384076843892, 1.0e-9, "ALPB Energy does not match"))
    goto error;
  if (!check(q[5], 0.0644849340, 1.0e-8, "ALPB Charge does not match"))
    goto error;
  if (!check(dipole[2], -0.3641866008, 1.0e-6, "ALPB Dipole does not match"))
    goto error;
  if (!check(wbo[9], 2.8932146955, 1.0e-8, "ALPB Bond order does not match"))
    goto error;

  // COSMO
  //Sensitive to guess, so reset the results
  xtb_delete(res);
  res = xtb_newResults();
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_loadGFN2xTB(env, mol, calc, NULL);
  xtb_setAccuracy(env, calc, 1.0);
  xtb_setElectronicTemp(env, calc, 300.0);
  xtb_setMaxIter(env, calc, 30);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_releaseSolvent(env, calc);
  xtb_setSolvent(env, calc, solvent, NULL, NULL, NULL, cosmo);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_singlepoint(env, mol, calc, res);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_getEnergy(env, res, &energy);
  xtb_getCharges(env, res, q);
  xtb_getDipole(env, res, dipole);
  xtb_getBondOrders(env, res, wbo);
  if (xtb_checkEnvironment(env))
    goto error;

  if (!check(energy, -8.385752224008, 1.0e-9, "COSMO Energy does not match"))
    goto error;
  if (!check(q[5], 0.0597931738, 1.0e-8, "COSMO Charge does not match"))
    goto error;
  if (!check(dipole[2], -0.3709356947, 1.0e-6, "COSMO Dipole does not match"))
    goto error;
  if (!check(wbo[9], 2.8940365143, 1.0e-8, "COSMO Bond order does not match"))
    goto error;


  // CPCMX
#if WITH_CPCMX
  //Sensitive to guess, so reset the results
  xtb_delete(res);
  res = xtb_newResults();
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_releaseSolvent(env, calc);
  xtb_setSolvent(env, calc, solvent, NULL, NULL, NULL, cpcmx);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_cpcmx_calc(env, mol, calc, res);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_getEnergy(env, res, &energy);
  xtb_getSolvationEnergy(env, res, &solv_energy);
  xtb_getCharges(env, res, q);
  xtb_getDipole(env, res, dipole);
  xtb_getBondOrders(env, res, wbo);
  if (xtb_checkEnvironment(env))
    goto error;

  if (!check(energy, -8.383520040765, 1.0e-9, "CPCM-X Energy does not match"))
    goto error;
  if (!check(solv_energy, -0.0010406558065, 1.0e-8, "CPCM-X Solvation Energy does not match"))
    goto error;
  if (!check(q[5], 0.0518386448, 1.0e-8, "CPCM-X Charge does not match"))
    goto error;
  if (!check(dipole[2], -0.2983150104, 1.0e-6, "CPCM-X Dipole does not match"))
    goto error;
  if (!check(wbo[9], 2.8982388543, 1.0e-8, "CPCM-X Bond order does not match"))
    goto error;

#endif
  // Compute Hessian
  //Sensitive to guess, so reset the results
  xtb_delete(res);
  res = xtb_newResults();
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_releaseSolvent(env, calc);
  xtb_hessian(env, mol, calc, res, hess, NULL, NULL, NULL, NULL);
  if (xtb_checkEnvironment(env))
    goto error;

  if (!check(hess[0], 0.4790088649, 1.0e-9, "Hessian[0,0] does not match"))
    goto error;
  if (!check(hess[3], -0.0528761190, 1.0e-9, "Hessian[0,3] does not match"))
    goto error;
  if (!check(hess[3], hess[63], 1.0e-9, "Hessian[0,3] != Hessian[3,0]"))
    goto error;
  if (!check(hess[(9*natsq)-1], 0.3636571159, 1.0e-9, "Hessian[21,21] does not match"))
    goto error;


  xtb_delete(res);
  xtb_delete(calc);
  xtb_delete(mol);
  xtb_delete(env);

  free(q);
  free(wbo);
  free(buffer);
  free(hess);

  tester = !res;
  if (!check(tester, 1, "Results not deleted"))
    goto error;
  tester = !calc;
  if (!check(tester, 1, "Calculator not deleted"))
    goto error;
  tester = !mol;
  if (!check(tester, 1, "Molecule not deleted"))
    goto error;
  tester = !env;
  if (!check(tester, 1, "Environment not deleted"))
    goto error;

  return 0;

error:
  xtb_showEnvironment(env, NULL);
  xtb_delete(res);
  xtb_delete(calc);
  xtb_delete(mol);
  xtb_delete(env);
  return 1;
}

int testSecond() {
  // water tetramer : https://xtb-docs.readthedocs.io/en/latest/pcem.html

  int const natoms = 6;
  int const attyp[6] = {8, 1, 1, 8, 1, 1};
  double const charge = 0.0;
  int const uhf = 0;
  double const coord[3 * 6] = {-2.75237178376284,     2.43247309226225,    -0.01392519847964,
                               -0.93157260886974,     2.79621404458590,    -0.01863384029005,
                               -3.43820531288547,     3.30583608421060,     1.42134539425148,
                               -2.43247309226225,    -2.75237178376284,     0.01392519847964,
                               -2.79621404458590,    -0.93157260886974,     0.01863384029005,
                               -3.30583608421060,    -3.43820531288547,    -1.42134539425148};

  int npc = 6;

  // external point charges
  double pc[6] = {-0.69645733, 0.36031084, 0.33614649,
                  -0.69645733, 0.36031084, 0.33614649};

  // coordinates of point charges
  double pccoord[3 * 6] = {2.75237178376284,    -2.43247309226225,    -0.01392519847964,
                           0.93157260886974,    -2.79621404458590,    -0.01863384029005,
                           3.43820531288547,    -3.30583608421060,     1.42134539425148,
                           2.43247309226225,     2.75237178376284,     0.01392519847964,
                           2.79621404458590,     0.93157260886974,     0.01863384029005,
                           3.30583608421060,     3.43820531288547,    -1.42134539425148};

  int numbers[6] = {8, 1, 1, 8, 1, 1};

  // results will live in here
  double energy = 0.0;
  double grad[3 * 6] = {0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0};

  double pcgrad[3 * 6] = {0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0};

  xtb_TEnvironment env = xtb_newEnvironment();
  xtb_TCalculator calc = xtb_newCalculator();
  xtb_TResults res = xtb_newResults();
  xtb_TMolecule mol =
      xtb_newMolecule(env, &natoms, attyp, coord, &charge, &uhf, NULL, NULL);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_setVerbosity(env, XTB_VERBOSITY_FULL);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_loadGFN2xTB(env, mol, calc, NULL);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_setExternalCharges(env, calc, &npc, numbers, pc, pccoord);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_singlepoint(env, mol, calc, res);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_getEnergy(env, res, &energy);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_getGradient(env, res, grad);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_getPCGradient(env, res, pcgrad);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_releaseExternalCharges(env, calc);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_delete(res);
  xtb_delete(calc);
  xtb_delete(mol);
  xtb_delete(env);

  if (!check(energy, -10.160927754235, 1.0e-9, "Energy does not match"))
    goto error;

  if (!check(grad[0 + 4 * 3], -2.0326750053811E-03, 1.0e-8, "grad[1,5] does not match"))
    goto error;
  if (!check(grad[1 + 1 * 3], 3.3125460324302E-03, 1.0e-8, "grad[2,2] does not match"))
    goto error;
  if (!check(grad[0 + 3 * 3], -1.1929405354823E-03, 1.0e-8, "grad[1,4] does not match"))
    goto error;
  if (!check(grad[2 + 5 * 3], -1.6607683158587E-03, 1.0e-8, "grad[3,6] does not match"))
    goto error;

  if (!check(pcgrad[0 + 4 * 3], -0.000248319582, 1.0e-8, "pcgrad[1,5] does not match"))
    goto error;
  if (!check(pcgrad[1 + 1 * 3], 0.001420844475, 1.0e-8, "pcgrad[2,2] does not match"))
    goto error;
  if (!check(pcgrad[0 + 3 * 3], 0.003746685275, 1.0e-8, "pcgrad[1,4] does not match"))
    goto error;
  if (!check(pcgrad[2 + 5 * 3], 0.000651613444, 1.0e-8, "pcgrad[3,6] does not match"))
    goto error;

  return 0;

error:
  xtb_showEnvironment(env, NULL);
  xtb_delete(res);
  xtb_delete(calc);
  xtb_delete(mol);
  xtb_delete(env);
  return 1;
}

int main(int argc, char **argv) {
  int stat = 0;
  stat += testFirst();
  stat += testSecond();
  return stat > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
