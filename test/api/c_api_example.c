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
   -0.23439724, 5.82196780, 12.45597529,
   -1.18337988, 6.25432521, 12.13871529,
   -0.97459304, 7.29844555, 11.90556842,
   -1.53287638, 5.69867346, 11.26851166,
   -2.18743690, 6.18914961, 13.29138602,
   -2.38832182, 5.14033352, 13.50965828,
   -3.51544014, 6.70782563, 12.73596360,
   -3.43813825, 7.76448415, 12.47983781,
   -4.31222592, 6.72367359, 13.47958380,
   -3.89122424, 6.13450845, 11.88842729,
   -1.71962790, 6.81419378, 14.64196118,
   -2.56489369, 6.73869928, 15.34188058,
   -1.30553616, 8.28708010, 14.46709643,
   -2.12153254, 8.88339645, 14.05885812,
   -0.43293661, 8.42096985, 13.82777397,
   -0.91009506, 8.84195389, 15.83709025,
   -0.62198472, 9.88210532, 15.68458920,
   -1.73049863, 8.91344164, 16.56640493,
    0.20835472, 8.06516304, 16.53450843,
    1.11027013, 8.08379059, 15.90505974,
    0.58406294, 8.62609252, 17.90749094,
   -0.27017230, 8.58733716, 18.58335654,
    0.94521765, 9.64819987, 17.79399575,
    1.34478256, 8.03159458, 18.41331915,
   -0.22842475, 6.60948035, 16.70992382,
    0.63592593, 6.02069012, 17.05055754,
   -1.04896081, 6.48967415, 17.43260784,
   -0.57931297, 6.05802189, 15.32665341,
    0.29823277, 6.06356948, 14.66342487,
   -1.03049345, 4.72745117, 15.59290274,
   -0.25892075, 4.19008917, 15.83503152};

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
   -0.897837, 14.879088, 29.374541,
   -0.512003, 15.284953, 28.580010,
   -0.503177, 16.714209, 28.534938,
    0.232489, 17.024434, 27.792783,
   -0.183554, 17.159648, 29.477095,
   -1.527540, 17.009207, 28.307266,
   13.961254, 11.508338, 12.790715,
   14.307185, 11.803974, 11.931815,
   13.896755, 10.872072, 10.927705,
   14.293925,  9.878504, 11.135478,
   12.808068, 10.885327, 10.875863,
   14.245566, 11.169547,  9.938799,
   13.606744,  9.874548, 15.190554,
   13.817846,  9.667827, 14.264552,
   14.639146,  8.517687, 14.045818,
   15.671299,  8.703176, 14.343150,
   14.107705,  7.702362, 14.536639,
   14.743162,  8.350417, 12.973787,
    6.716485, -6.346264, 18.679919,
    5.824101, -6.265903, 19.056785,
    5.047806, -5.789462, 17.954393,
    5.386539, -6.101654, 16.966443,
    4.894122, -4.710696, 17.927294,
    4.068575, -6.261430, 18.035357,
   -6.316725, 11.850365, 28.030715,
   -5.539369, 11.280064, 28.154212,
   -4.736778, 11.127486, 26.980712,
   -4.402790, 12.131419, 26.718633,
   -5.410452, 10.748895, 26.211847,
   -3.938730, 10.404747, 27.150584,
    8.135514, 15.983725, 28.508396,
    7.989254, 16.686992, 27.853495};

  int numbers[32] = {
    99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 
    99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 
    99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 
    99, 99};

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

  assert(fabs(pcgrad[0] + 0.00001022) < 1.0e-6);
  assert(fabs(pcgrad[95] - 0.00001225) < 1.0e-6);
}

int main(int argc, char **argv) {
   testFirst();
   testSecond();
   return 0;
}
