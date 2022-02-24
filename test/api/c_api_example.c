#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "xtb.h"

int testFirst() {
   // first molecule
   
   double const thr = 1.0e-10;
   int    const natoms = 7;
   int    const attyp[7] = {6,6,6,1,1,1,1};
   double const charge = 0.0;
   int    const uhf = 0;
   double const coord[3*7] =
      {0.00000000000000, 0.00000000000000,-1.79755622305860,
       0.00000000000000, 0.00000000000000, 0.95338756106749,
       0.00000000000000, 0.00000000000000, 3.22281255790261,
      -0.96412815539807,-1.66991895015711,-2.53624948351102,
      -0.96412815539807, 1.66991895015711,-2.53624948351102,
       1.92825631079613, 0.00000000000000,-2.53624948351102,
       0.00000000000000, 0.00000000000000, 5.23010455462158};

   xtb_TEnvironment env;
   xtb_TMolecule mol;
   xtb_TCalculator calc;
   xtb_TResults res;
   double energy;
   double dipole[3];
   double q[natoms];
   double wbo[natoms*natoms];
   int buffersize = 512;
   char buffer[buffersize];
   char solvent[] = "h2o";

   assert(XTB_API_VERSION == xtb_getAPIVersion());

   env = xtb_newEnvironment();
   calc = xtb_newCalculator();
   res = xtb_newResults();
   mol = xtb_newMolecule(env, &natoms, attyp, coord, NULL, NULL, NULL, NULL);
   if (xtb_checkEnvironment(env)) {
      xtb_showEnvironment(env, NULL);
      return 1;
   }

   xtb_getEnergy(env, res, &energy);
   if (xtb_checkEnvironment(env)) {
      xtb_getError(env, buffer, &buffersize);
      printf("Error message is:\n%s\n", buffer);
   }

   xtb_setVerbosity(env, XTB_VERBOSITY_FULL);
   if (xtb_checkEnvironment(env)) {
      xtb_showEnvironment(env, NULL);
      return 2;
   }

   xtb_loadGFN2xTB(env, mol, calc, NULL);
   xtb_setAccuracy(env, calc, 1.0);
   xtb_setElectronicTemp(env, calc, 300.0);
   xtb_setMaxIter(env, calc, 30);
   if (xtb_checkEnvironment(env)) {
      xtb_showEnvironment(env, NULL);
      return 3;
   }

   xtb_singlepoint(env, mol, calc, res);
   if (xtb_checkEnvironment(env)) {
      xtb_showEnvironment(env, NULL);
      return 4;
   }

   xtb_getEnergy(env, res, &energy);
   xtb_getCharges(env, res, q);
   xtb_getDipole(env, res, dipole);
   xtb_getBondOrders(env, res, wbo);
   if (xtb_checkEnvironment(env)) {
      xtb_showEnvironment(env, NULL);
      return 5;
   }

   assert(fabs(energy + 8.3824793849585) < 1.0e-9);
   assert(fabs(q[5] - 0.05184019996829) < 1.0e-8);
   assert(fabs(dipole[2] + 0.298279305689518) < 1.0e-6);
   assert(fabs(wbo[9] - 2.89823984265213) < 1.0e-8);

   xtb_setSolvent(env, calc, solvent, NULL, NULL, NULL);
   if (xtb_checkEnvironment(env)) {
      xtb_showEnvironment(env, NULL);
      return 6;
   }

   xtb_singlepoint(env, mol, calc, res);
   if (xtb_checkEnvironment(env)) {
      xtb_showEnvironment(env, NULL);
      return 7;
   }

   xtb_getEnergy(env, res, &energy);
   xtb_getCharges(env, res, q);
   xtb_getDipole(env, res, dipole);
   xtb_getBondOrders(env, res, wbo);
   if (xtb_checkEnvironment(env)) {
      xtb_showEnvironment(env, NULL);
      return 8;
   }

   assert(fabs(energy + 8.38393864716134) < 1.0e-9);
   assert(fabs(q[5] - 0.06090868805034) < 1.0e-8);
   assert(fabs(dipole[2] + 0.35455233974705) < 1.0e-6);
   assert(fabs(wbo[9] - 2.89453979224265) < 1.0e-8);

   xtb_delResults(&res);
   xtb_delCalculator(&calc);
   xtb_delMolecule(&mol);
   xtb_delEnvironment(&env);

   assert(!res);
   assert(!calc);
   assert(!mol);
   assert(!env);
}

int testSecond() {
  // second molecule

  int    const natoms = 31;
  int    const attyp[31] = {1, 6, 1, 1, 6, 1, 6, 1, 1, 1, 6, 1, 6, 1, 1, 6, 1, 1, 6, 1, 6, 1, 1, 1, 6, 1, 1, 6, 1, 8, 1};
  double const charge = 0.0;
  int    const uhf = 0;
  double const coord[3*31] = {
   -0.442947496, 11.00194721, 23.53843018,
   -2.236268461, 11.81898598, 22.93889444,
   -1.841717705, 13.79209151, 22.49830981,
   -2.896722482, 10.7689542, 21.29444454,
   -4.133665133, 11.69582169, 25.1171309,
   -4.513283393, 9.713842463, 25.52960654,
   -6.643232696, 12.67597933, 24.06753249,
   -6.497152995, 14.67277863, 23.5835239,
   -8.148942688, 12.70592769, 25.47277389,
   -7.353363183, 11.59256465, 22.4659177,
   -3.249632431, 12.87698641, 27.6693533,
   -4.846956553, 12.73432219, 28.99201199,
   -2.467110848, 15.66034388, 27.33890614,
   -4.009123687, 16.78722077, 26.56744596,
   -0.8181333, 15.91335935, 26.1307593,
   -1.719833938, 16.70890552, 29.92782456,
   -1.175383185, 18.67451089, 29.63963875,
   -3.270175176, 16.84399807, 31.30603239,
   0.393734165, 15.24098055, 31.24575662,
   2.098110773, 15.27618159, 30.05626854,
   1.10372126, 16.30098582, 33.84032285,
   -0.5105527, 16.22774865, 35.11752635,
   1.78620615, 18.23249274, 33.62584759,
   2.541275947, 15.17754523, 34.7962016,
   -0.431661103, 12.4901333, 31.57724434,
   1.201728308, 11.37747874, 32.2209501,
   -1.982252711, 12.26373193, 32.94292201,
   -1.094745099, 11.44802571, 28.96323675,
   0.563579412, 11.45850915, 27.70991388,
   -1.947354387, 8.933606299, 29.46637609,
   -0.489290309, 7.918137207, 29.92393411};

  double energy = 0;
  
  int    npc = 32;

  // external point charges
  double pc[32] = {0.431000,-0.664000, 0.173000, 0.020000, 0.020000, 0.020000, 
    0.431000,-0.664000, 0.173000, 0.020000, 0.020000, 0.020000, 0.431000,-0.664000, 
    0.173000, 0.020000, 0.020000, 0.020000, 0.431000,-0.664000, 0.173000, 0.020000, 
    0.020000, 0.020000, 0.431000,-0.664000, 0.173000, 0.020000, 0.020000, 0.020000, 
    0.431000,-0.664000};

  // coordinates of point charges
  double pccoord[3*32] = {
   -1.696669514, 28.11745897, 55.50995136,
   -0.967547429, 28.88443423, 54.0085023,
   -0.950868672, 31.58534217, 53.92332839,
   0.439341438, 32.17158366, 52.52085582,
   -0.3468675, 32.42710162, 55.70375073,
   -2.886638164, 32.14280874, 53.49308978,
   26.38300052, 21.74765157, 24.17099786,
   27.03671671, 22.30632379, 22.54790876,
   26.26111483, 20.54528062, 20.65041197,
   27.01165889, 18.66770536, 21.04304684,
   24.20379034, 20.57032899, 20.55244459,
   26.92027344, 21.10742805, 18.78164663,
   25.71307234, 18.66022959, 28.70604561,
   26.11199812, 18.26958272, 26.95615185,
   27.66403337, 16.09612865, 26.54280365,
   29.61452386, 16.44665278, 27.10468085,
   26.65975337, 14.55538454, 27.47032282,
   27.86059553, 15.78003352, 24.51695451,
   12.6923432, -11.99272547, 35.30000333,
   11.00597838, -11.84086488, 36.01217832,
   9.538990432, -10.94052003, 33.92895508,
   10.17910434, -11.53047861, 32.06199633,
   9.248569167, -8.901943552, 33.87774529,
   7.688508235, -11.83241211, 34.08195518,
   -11.93690473, 22.39399025, 52.97048306,
   -10.46791178, 21.31627534, 53.20385904,
   -8.95123149, 21.02794412, 50.98626089,
   -8.320084347, 22.92510643, 50.49100234,
   -10.22429346, 20.31250935, 49.53331363,
   -7.443136243, 19.66216255, 51.3072731,
   15.37392487, 30.20492464, 53.87317117,
   15.09753296, 31.53390939, 52.63558511};

  int numbers[32] = {
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
    7, 7};

  // results will live in here
  double pcgrad[3*32] =
      {0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0};

  xtb_TEnvironment env = xtb_newEnvironment();
  xtb_TCalculator calc = xtb_newCalculator();
  xtb_TResults res = xtb_newResults();
  xtb_TMolecule mol = xtb_newMolecule(
      env, &natoms, attyp, coord, &charge, &uhf, NULL, NULL);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return 9;
  }

  xtb_setVerbosity(env, XTB_VERBOSITY_FULL);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return 10;
  }

  xtb_loadGFN1xTB(env, mol, calc, NULL);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return 11;
  }

  xtb_setExternalCharges(env, calc, &npc, numbers, pc, pccoord);
   if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return 12;
   }

  xtb_singlepoint(env, mol, calc, res);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return 13;
  }

  xtb_getPCGradient(env, res, pcgrad);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return 14;
  }

  xtb_releaseExternalCharges(env, calc);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return 15;
  }
  
  xtb_delResults(&res);
  xtb_delCalculator(&calc);
  xtb_delMolecule(&mol);
  xtb_delEnvironment(&env);

  assert(fabs(pcgrad[0] - 0.00000755) < 1.0e-6);
  assert(fabs(pcgrad[95] - 0.00001312) < 1.0e-6);
}

int main(int argc, char **argv) {
   testFirst();
   testSecond();
   return 0;
}
