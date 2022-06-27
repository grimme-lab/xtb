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
  // second molecule

  int const natoms = 31;
  int const attyp[31] = {1, 6, 1, 1, 6, 1, 6, 1, 1, 1, 6, 1, 6, 1, 1, 6,
                         1, 1, 6, 1, 6, 1, 1, 1, 6, 1, 1, 6, 1, 8, 1};
  double const charge = 0.0;
  int const uhf = 0;
  double const coord[3 * 31] = {
      -0.442947496, 11.001947210, 23.53843018,  -2.236268461, 11.818985980,
      22.93889444,  -1.841717705, 13.792091510, 22.49830981,  -2.896722482,
      10.768954200, 21.29444454,  -4.133665133, 11.695821690, 25.11713090,
      -4.513283393, 9.713842463,  25.52960654,  -6.643232696, 12.675979330,
      24.06753249,  -6.497152995, 14.672778630, 23.58352390,  -8.148942688,
      12.705927690, 25.47277389,  -7.353363183, 11.592564650, 22.46591770,
      -3.249632431, 12.876986410, 27.66935330,  -4.846956553, 12.734322190,
      28.99201199,  -2.467110848, 15.660343880, 27.33890614,  -4.009123687,
      16.787220770, 26.56744596,  -0.818133300, 15.913359350, 26.13075930,
      -1.719833938, 16.708905520, 29.92782456,  -1.175383185, 18.674510890,
      29.63963875,  -3.270175176, 16.843998070, 31.30603239,  0.393734165,
      15.240980550, 31.24575662,  2.098110773,  15.276181590, 30.05626854,
      1.103721260,  16.300985820, 33.84032285,  -0.510552700, 16.227748650,
      35.11752635,  1.786206150,  18.232492740, 33.62584759,  2.541275947,
      15.177545230, 34.79620160,  -0.431661103, 12.490133300, 31.57724434,
      1.201728308,  11.377478740, 32.22095010,  -1.982252711, 12.263731930,
      32.94292201,  -1.094745099, 11.448025710, 28.96323675,  0.563579412,
      11.458509150, 27.70991388,  -1.947354387, 8.933606299,  29.46637609,
      -0.489290309, 7.918137207,  29.92393411};

  int npc = 32;

  // external point charges
  double pc[32] = {0.431000, -0.664000, 0.173000, 0.020000, 0.020000, 0.020000,
                   0.431000, -0.664000, 0.173000, 0.020000, 0.020000, 0.020000,
                   0.431000, -0.664000, 0.173000, 0.020000, 0.020000, 0.020000,
                   0.431000, -0.664000, 0.173000, 0.020000, 0.020000, 0.020000,
                   0.431000, -0.664000, 0.173000, 0.020000, 0.020000, 0.020000,
                   0.431000, -0.664000};

  // coordinates of point charges
  double pccoord[3 * 32] = {
      -1.696669514,  28.11745897,  55.50995136,   -0.967547429, 28.88443423,
      54.00850230,   -0.950868672, 31.58534217,   53.92332839,  0.439341438,
      32.17158366,   52.52085582,  -0.346867500,  32.42710162,  55.70375073,
      -2.886638164,  32.14280874,  53.49308978,   26.383000520, 21.74765157,
      24.17099786,   27.036716710, 22.30632379,   22.54790876,  26.261114830,
      20.54528062,   20.65041197,  27.011658890,  18.66770536,  21.04304684,
      24.203790340,  20.57032899,  20.55244459,   26.920273440, 21.10742805,
      18.78164663,   25.713072340, 18.66022959,   28.70604561,  26.111998120,
      18.26958272,   26.95615185,  27.664033370,  16.09612865,  26.54280365,
      29.614523860,  16.44665278,  27.10468085,   26.659753370, 14.55538454,
      27.47032282,   27.860595530, 15.78003352,   24.51695451,  12.692343200,
      -11.99272547,  35.30000333,  11.005978380,  -11.84086488, 36.01217832,
      9.538990432,   -10.94052003, 33.92895508,   10.179104340, -11.53047861,
      32.06199633,   9.248569167,  -8.901943552,  33.87774529,  7.688508235,
      -11.83241211,  34.08195518,  -11.936904730, 22.39399025,  52.97048306,
      -10.467911780, 21.31627534,  53.20385904,   -8.951231490, 21.02794412,
      50.98626089,   -8.320084347, 22.92510643,   50.49100234,  -10.224293460,
      20.31250935,   49.53331363,  -7.443136243,  19.66216255,  51.30727310,
      15.373924870,  30.20492464,  53.87317117,   15.097532960, 31.53390939,
      52.63558511};

  int numbers[32] = {7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                     7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};

  // results will live in here
  double pcgrad[3 * 32] = {
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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

  xtb_loadGFN1xTB(env, mol, calc, NULL);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_setExternalCharges(env, calc, &npc, numbers, pc, pccoord);
  if (xtb_checkEnvironment(env))
    goto error;

  xtb_singlepoint(env, mol, calc, res);
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

  if (!check(pcgrad[0], 0.00000755, 1.0e-6, "pcgrad[0] does not match"))
    goto error;
  if (!check(pcgrad[95], 0.00001312, 1.0e-6, "pcgrad[95] does not match"))
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
