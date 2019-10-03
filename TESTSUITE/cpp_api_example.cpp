#include <cmath>
#include <cassert>
#include <cstdio>

#include "xtb.h"

using namespace std;

int
main (int argc, char **argv)
{
   const double thr = 1.0e-10;
   const int    natoms = 7;
   const int    attyp[natoms] = {6,6,6,1,1,1,1};
   const double charge = 0.0;
   const double coord[3*natoms] =
      {0.00000000000000, 0.00000000000000,-1.79755622305860,
       0.00000000000000, 0.00000000000000, 0.95338756106749,
       0.00000000000000, 0.00000000000000, 3.22281255790261,
      -0.96412815539807,-1.66991895015711,-2.53624948351102,
      -0.96412815539807, 1.66991895015711,-2.53624948351102,
       1.92825631079613, 0.00000000000000,-2.53624948351102,
       0.00000000000000, 0.00000000000000, 5.23010455462158};

   const xtb::SCC_options opt {2, 0, 1.0, 300.0, true, false, true, 30, "none"};

   double energy {0.0};
   double dipole[3] {0.0};
   double q[natoms] {0.0};
   double qp[6*natoms] {0.0};
   double wbo[natoms*natoms] {0.0};

   int stat = xtb::GFN2_calculation(&natoms, attyp, &charge, nullptr, coord, &opt, "-",
                                    &energy, nullptr, dipole, q, nullptr, qp, wbo);

   assert(stat == 0);
   assert(fabs(energy + 8.3824793818504) < 1.0e-10);
   assert(fabs(q[5] - 0.05184019996829) < 1.0e-8);
   assert(fabs(dipole[2] + 0.29832384492435) < 1.0e-6);
   assert(fabs(qp[14] - 0.56386138525354) < 1.0e-6);
   assert(fabs(wbo[9] - 2.89823984265213) < 1.0e-8);

   return 0;
}

