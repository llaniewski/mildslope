#include <stddef.h>
#include <math.h>

extern size_t * triangles;
extern size_t n_triangles;
extern size_t n_points;

void morph_energy(const double *P0, const double* P1, double* energy) {
    energy[0] = 0;
    energy[1] = 0;
    // $AD II-LOOP
    for (size_t i=0;i<n_triangles;i++) {
        size_t i0 = triangles[i*3+0];
        size_t i1 = triangles[i*3+1];
        size_t i2 = triangles[i*3+2];
        double zX00 = -P0[i0*2+0] + P0[i1*2+0];
        double zX10 = -P0[i0*2+1] + P0[i1*2+1];
        double zX01 = -P0[i0*2+0] + P0[i2*2+0];
        double zX11 = -P0[i0*2+1] + P0[i2*2+1];
        double zx00 = -P1[i0*2+0] + P1[i1*2+0];
        double zx10 = -P1[i0*2+1] + P1[i1*2+1];
        double zx01 = -P1[i0*2+0] + P1[i2*2+0];
        double zx11 = -P1[i0*2+1] + P1[i2*2+1];
        double det = -zX01*zX10 + zX11*zX00;
        double Xx00 = ( -zX10*zx01 + zX11*zx00 )/det;
        double Xx10 = ( -zX10*zx11 + zX11*zx10 )/det;
        double Xx01 = ( zX00*zx01 - zX01*zx00 )/det;
        double Xx11 = ( zX00*zx11 - zX01*zx10 )/det;
        double SXX00 = Xx10*Xx10 + Xx00*Xx00;
        double SXX10 = Xx10*Xx11 + Xx00*Xx01;
        double SXX01 = Xx11*Xx10 + Xx01*Xx00;
        double SXX11 = Xx11*Xx11 + Xx01*Xx01;
        double F00 = -1 + Xx10*Xx10 + Xx00*Xx00;
        double F10 = Xx10*Xx11 + Xx00*Xx01;
        double F01 = Xx11*Xx10 + Xx01*Xx00;
        double F11 = -1 + Xx11*Xx11 + Xx01*Xx01;
        double J = Xx00 * Xx11 - Xx10 * Xx01;
        double I1 = F00 + F11;
        double area = -det/2;
        // energy[0] = energy[0] + area * (F00 + F11)*(F00 + F11);
        // energy[1] = energy[1] + area * (F00*F00+2*F01*F10+F11*F11);
        energy[0] = energy[0] + area * (I1 - 3 - 2*log(J));
        energy[1] = energy[1] + area * (J - 1) * (J-1);
    }
}
