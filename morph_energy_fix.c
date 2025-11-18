#include <stddef.h>
#include <math.h>

extern size_t n_bord;
extern size_t* bord;
extern double* fix;
extern double* fix_dirs;

void morph_energy_fix(const double *P0, const double *P1, const double *Pfix, double *res, double *energy, double *energyb) {
    double energy[2];
    morph_energy_b(P0, P1, res, energy, energyb);
    // $AD II-LOOP
    for (size_t i=0; i<n_bord; i++) {
        size_t idx = bord[i];
        res[idx] = (1-fix[i])*res[i] + fix[i]*(P0[idx] + Pfix[i] - P1[idx]);
    }
}
