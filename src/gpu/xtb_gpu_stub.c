/* This file is part of xtb.
 *
 * xtb is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * xtb is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with xtb.  If not, see <https://www.gnu.org/licenses/>.
 */

/* Stand-ins used when xtb is configured without -Dgpu=true.  Every entry
 * point reports failure, which the Fortran callers read as "use the CPU
 * path", so the rest of the source needs no preprocessor guards. */

#include <string.h>

int xtb_gpu_compiled(void) { return 0; }

int xtb_gpu_open(int dev, int verbose) { (void)dev; (void)verbose; return -1; }

void xtb_gpu_name(char *buf, int buflen) {
   if (buflen > 0) strncpy(buf, "not compiled", (size_t)buflen - 1);
   if (buflen > 0) buf[buflen - 1] = '\0';
}

int xtb_gpu_device(void) { return -1; }
void xtb_gpu_ramp_reset(void) {}
int xtb_gpu_accept_density(void) { return 1; }
void xtb_gpu_report(void) {}

int xtb_gpu_sygvd(int n, double *a, const double *b, double *w, double tol) {
   (void)n; (void)a; (void)b; (void)w; (void)tol; return -1;
}

int xtb_gpu_syevd(int n, double *a, double *w) {
   (void)n; (void)a; (void)w; return -1;
}
