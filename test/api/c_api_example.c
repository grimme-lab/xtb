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
  fprintf(stderr, "FAIL: %s: expected %g, got %g\n", msg, expected, actual);
  return false;
}

int testFirst() {
  // first molecule
  double* q;
  char* buffer;
  double* wbo;

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
  double dipole[3];
  q = (double*) malloc(natoms * sizeof(double));
  wbo = (double*) malloc(natsq * sizeof(double));
  buffer = (char*) malloc(buffersize *sizeof(char));
  char solvent[] = "h2o";

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

  xtb_setSolvent(env, calc, solvent, NULL, NULL, NULL);
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

  if (!check(energy, -8.38393864716134, 1.0e-9, "Energy does not match"))
    goto error;
  if (!check(q[5], 0.06090868805034, 1.0e-8, "Charge does not match"))
    goto error;
  if (!check(dipole[2], -0.35455233974705, 1.0e-6, "Dipole does not match"))
    goto error;
  if (!check(wbo[9], +2.89453979224265, 1.0e-8, "Bond order does not match"))
    goto error;

  xtb_delete(res);
  xtb_delete(calc);
  xtb_delete(mol);
  xtb_delete(env);

  free(q);
  free(wbo);
  free(buffer);

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
