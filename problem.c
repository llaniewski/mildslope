#include <stddef.h>
#include <math.h>

extern size_t * triangles;
extern size_t n_triangles;
extern size_t n_points;
extern size_t * boundary;
extern int * boundary_flag;
extern size_t n_boundary;

void problem(double wave_k, const double *points, const double *depth, const double* x, double* res, double* obj) {
    double M00,M01,M02,M11,M12,M22;
    M00=M11=M22=1.0/6.0;
    M01=M02=M12=1.0/12.0;
    
    for (size_t i=0;i<2*n_points;i++) res[i] = 0;
    obj[0] = 0;
    // $AD II-LOOP
    for (size_t i=0;i<n_triangles;i++) {
        size_t i0 = triangles[i*3+0];
        size_t i1 = triangles[i*3+1];
        size_t i2 = triangles[i*3+2];
        double x0 = points[i0*2+0];
        double y0 = points[i0*2+1];
        double x1 = points[i1*2+0];
        double y1 = points[i1*2+1];
        double x2 = points[i2*2+0];
        double y2 = points[i2*2+1];
        double vx = x1-x0;
        double vy = y1-y0;
        double wx = x2-x0;
        double wy = y2-y0;
        double det = wx*vy-wy*vx;
        double area = -det/2;
        double dp = (depth[i0]+depth[i1]+depth[i2])/3.0;
        //vx*gx1+vy*gy1 = 1
        //wx*gx1+wy*gy1 = 0
        //gx1 = -wy*gy1/wx
        //-wy*vx*gy1/wx+vy*gy1 = 1
        //(-wy*vx/wx+wx*vy/wx)*gy1 = 1
        //(-wy*vx+wx*vy)/wx*gy1 = 1
        //gy1 = wx / (wx*vy-wy*vx)
        //gx1 = -wy / (wx*vy-wy*vx)
        double gx1 = -wy / det;
        double gy1 = wx / det;
        double gx2 = vy / det;
        double gy2 = -vx / det;
        double gx0 = - gx1 - gx2;
        double gy0 = - gy1 - gy2;
        double K00 = gx0*gx0 + gy0*gy0;
        double K01 = gx0*gx1 + gy0*gy1;
        double K02 = gx0*gx2 + gy0*gy2;
        double K11 = gx1*gx1 + gy1*gy1;
        double K12 = gx1*gx2 + gy1*gy2;
        double K22 = gx2*gx2 + gy2*gy2;
        double a = -area*wave_k*wave_k;
        double b = area*dp;
        res[i0*2+0] += a*(M00*x[i0*2+0] + M01*x[i1*2+0] + M02*x[i2*2+0]);
        res[i1*2+0] += a*(M01*x[i0*2+0] + M11*x[i1*2+0] + M12*x[i2*2+0]);
        res[i2*2+0] += a*(M02*x[i0*2+0] + M12*x[i1*2+0] + M22*x[i2*2+0]);
        res[i0*2+0] += b*(K00*x[i0*2+0] + K01*x[i1*2+0] + K02*x[i2*2+0]);
        res[i1*2+0] += b*(K01*x[i0*2+0] + K11*x[i1*2+0] + K12*x[i2*2+0]);
        res[i2*2+0] += b*(K02*x[i0*2+0] + K12*x[i1*2+0] + K22*x[i2*2+0]);
        res[i0*2+1] += a*(M00*x[i0*2+1] + M01*x[i1*2+1] + M02*x[i2*2+1]);
        res[i1*2+1] += a*(M01*x[i0*2+1] + M11*x[i1*2+1] + M12*x[i2*2+1]);
        res[i2*2+1] += a*(M02*x[i0*2+1] + M12*x[i1*2+1] + M22*x[i2*2+1]);
        res[i0*2+1] += b*(K00*x[i0*2+1] + K01*x[i1*2+1] + K02*x[i2*2+1]);
        res[i1*2+1] += b*(K01*x[i0*2+1] + K11*x[i1*2+1] + K12*x[i2*2+1]);
        res[i2*2+1] += b*(K02*x[i0*2+1] + K12*x[i1*2+1] + K22*x[i2*2+1]);
    }
    // $AD II-LOOP
    for (size_t i=0;i<n_boundary;i++) {
        size_t i0 = boundary[i*2+0];
        size_t i1 = boundary[i*2+1];
        double x0 = points[i0*2+0];
        double y0 = points[i0*2+1];
        double x1 = points[i1*2+0];
        double y1 = points[i1*2+1];
        double vx = x1-x0;
        double vy = y1-y0;
        double len = sqrt(vx*vx+vy*vy);
        double b = len;
        double robin = wave_k;
        double EM00, EM01, EM11;
        EM00 = EM11 = 1.0/3.0;
        EM01 = 1.0/6.0;
        double val0 = 0;
        double val1 = 0;
        double o0, o1;
        if (boundary_flag[i] == 2) val1 = 1;
        if (boundary_flag[i] > 1) {
            res[i0*2+0] +=  b*robin*(EM00 * x[i0*2+1] + EM01 * x[i1*2+1] - val0);
            res[i1*2+0] +=  b*robin*(EM01 * x[i0*2+1] + EM11 * x[i1*2+1] - val0);
            res[i0*2+1] += -b*robin*(EM00 * x[i0*2+0] + EM01 * x[i1*2+0] - val1);
            res[i1*2+1] += -b*robin*(EM01 * x[i0*2+0] + EM11 * x[i1*2+0] - val1);
            o0 = (x[i0*2+0] - val1)*(x[i0*2+0] - val1) + (x[i0*2+1] - val0)*(x[i0*2+1] - val0);
            o1 = (x[i1*2+0] - val1)*(x[i1*2+0] - val1) + (x[i1*2+1] - val0)*(x[i1*2+1] - val0);
            obj[boundary_flag[i]-2] += 0.5*b*(o0 + o1);
        }
    }
}
