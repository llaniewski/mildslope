#include <stddef.h>
#include <math.h>

extern size_t * triangles;
extern size_t n_triangles;
extern size_t n_points;

void morph_energy(const double *P0, const double* P1, double* energy) {
    energy[0] = 0;
    // $AD II-LOOP
    for (size_t i=0;i<n_triangles;i++) {
        size_t i0 = triangles[i*3+0];
        size_t i1 = triangles[i*3+1];
        size_t i2 = triangles[i*3+2];
        double vx0 = P0[i1*2+0] - P0[i0*2+0];
        double vy0 = P0[i1*2+1] - P0[i0*2+1];
        double wx0 = P0[i2*2+0] - P0[i0*2+0];
        double wy0 = P0[i2*2+1] - P0[i0*2+1];
        double vx1 = P1[i1*2+0] - P1[i0*2+0];
        double vy1 = P1[i1*2+1] - P1[i0*2+1];
        double wx1 = P1[i2*2+0] - P1[i0*2+0];
        double wy1 = P1[i2*2+1] - P1[i0*2+1];
        double det = vx0*wy0 - vy0*wx0;
        double detinv = 1/det;
        double mxx = detinv*(-vy1*wx0 + vx1*wy0) - 1;
        double mxy = detinv*( vy1*vx0 - vx1*vy0);
        double myx = detinv*(-wy1*wx0 + wx1*wy0);
        double myy = detinv*( wy1*vx0 - wx1*vy0) - 1;
        double area = -det/2;
        energy[0] = energy[0] + area * (mxx*mxx + myy*myy);
    }
}
