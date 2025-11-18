
void morph_energy_fix(const double *P0, const double *P1, const double *Pfix, double *res, double *energy, double *energyb) {
    double energy[2];
    morph_energy_b(P0, P1, res, energy, energyb);
}
