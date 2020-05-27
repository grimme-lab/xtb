#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "xtb.h"

int
main (int argc, char **argv)
{
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

   return 0;
}
