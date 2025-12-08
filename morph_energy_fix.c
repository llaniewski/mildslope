#include <stddef.h>
#include <math.h>

extern size_t n_bord;
extern size_t* bord;
extern double* fix;
extern double* fix_dirs;

void morph_energy_fix(const double *P0, const double *P1, const double *dir_disp, double *res, double *energyb) {
    double energy[2];
    morph_energy_b(P0, P1, res, energy, energyb);
    // $AD II-LOOP
    for (size_t i=0; i<n_bord; i++) {
        size_t idx = bord[i];
        size_t ix = idx*2 + 0;
        size_t iy = idx*2 + 1;
        double vx = fix_dirs[i*2+0];
        double vy = fix_dirs[i*2+1];
        double dispx = P1[ix] - P0[ix];
        double dispy = P1[iy] - P0[iy];
        double disp0 = dispx*vx + dispy*vy;
        double disp1 = -dispx*vy + dispy*vx;
        double res0 =  res[ix]*vx + res[iy]*vy;
        double res1 = -res[ix]*vy + res[iy]*vx;
        res0 = res0 * (1-fix[2*i+0]) + (disp0 - dir_disp[2*i+0]) * fix[2*i+0];
        res1 = res1 * (1-fix[2*i+1]) + (disp1 - dir_disp[2*i+1]) * fix[2*i+1];
        res[ix] = res0*vx - res1*vy;
        res[iy] = res0*vy + res1*vx;
    }
}
