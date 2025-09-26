#include "global.h"
#include "vtu_write.h"
#include "solve.h"

#include <math.h>

std::vector<double> P;
size_t nP;
std::vector<size_t> T;
size_t nT;
std::vector<size_t> B1;
std::vector<size_t> B2;
size_t nB1,nB2;

void MMult(const double* x, double* Mx) {
    double M00,M01,M02,M11,M12,M22;
    M00=M11=M22=1.0/6.0;
    M01=M02=M12=1.0/12.0;
    for (size_t i=0;i<2*nP;i++) Mx[i] = 0;
    for (size_t i=0;i<nT;i++) {
        size_t i0 = T[i*3+0];
        size_t i1 = T[i*3+1];
        size_t i2 = T[i*3+2];
        double x0 = P[i0*2+0];
        double y0 = P[i0*2+1];
        double x1 = P[i1*2+0];
        double y1 = P[i1*2+1];
        double x2 = P[i2*2+0];
        double y2 = P[i2*2+1];
        double vx = x1-x0;
        double vy = y1-y0;
        double wx = x2-x0;
        double wy = y2-y0;
        double det = wx*vy-wy*vx;
        double area = fabs(det)/2;
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
        double a = area;
        double k = 1;
        double b = area;
        Mx[i0*2+0] += a*(M00*x[i0*2+0] + M01*x[i1*2+0] + M02*x[i2*2+0]);
        Mx[i1*2+0] += a*(M01*x[i0*2+0] + M11*x[i1*2+0] + M12*x[i2*2+0]);
        Mx[i2*2+0] += a*(M02*x[i0*2+0] + M12*x[i1*2+0] + M22*x[i2*2+0]);
        Mx[i0*2+0] += b*(K00*x[i0*2+0] + K01*x[i1*2+0] + K02*x[i2*2+0]);
        Mx[i1*2+0] += b*(K01*x[i0*2+0] + K11*x[i1*2+0] + K12*x[i2*2+0]);
        Mx[i2*2+0] += b*(K02*x[i0*2+0] + K12*x[i1*2+0] + K22*x[i2*2+0]);
        Mx[i0*2+1] += a*(M00*x[i0*2+1] + M01*x[i1*2+1] + M02*x[i2*2+1]);
        Mx[i1*2+1] += a*(M01*x[i0*2+1] + M11*x[i1*2+1] + M12*x[i2*2+1]);
        Mx[i2*2+1] += a*(M02*x[i0*2+1] + M12*x[i1*2+1] + M22*x[i2*2+1]);
        Mx[i0*2+1] += b*(K00*x[i0*2+1] + K01*x[i1*2+1] + K02*x[i2*2+1]);
        Mx[i1*2+1] += b*(K01*x[i0*2+1] + K11*x[i1*2+1] + K12*x[i2*2+1]);
        Mx[i2*2+1] += b*(K02*x[i0*2+1] + K12*x[i1*2+1] + K22*x[i2*2+1]);
    }
    for (size_t i=0;i<nB1;i++) {
        size_t idx = B1[i];
        Mx[idx*2+0] = x[idx*2+0];
        Mx[idx*2+1] = x[idx*2+1];
    }
    for (size_t i=0;i<nB2;i++) {
        size_t idx = B2[i];
        Mx[idx*2+0] = x[idx*2+0];
        Mx[idx*2+1] = x[idx*2+1];
    }
}



int main() {
    {
        FILE* f = fopen("mesh/mesh2_points.txt","rb");
        nP = 0;
        while (! feof(f)) {
            double x,y;
            fscanf(f, "%lf %lf",&x,&y);
            P.push_back(x);
            P.push_back(y);
            nP++;
        }
        fclose(f);
    }
    printf("Points %ld\n", P.size());
    {
        double maxx=P[0], minx=P[0];
        for (size_t i=0;i<nP;i++) {
            if (P[i*2+0] > maxx) maxx = P[i*2+0];
            if (P[i*2+0] < minx) minx = P[i*2+0];
        }
        printf("X: [%lg, %lg]\n", minx, maxx);
        for (size_t i=0;i<nP;i++) {
            if (fabs(P[i*2+0] - minx) < 1e-6) B1.push_back(i);
            if (fabs(P[i*2+0] - maxx) < 1e-6) B2.push_back(i);
        }
    }
    nB1 = B1.size();
    nB2 = B2.size();
    printf("Borders: %ld, %ld\n", B1.size(), B2.size());
    {
        FILE* f = fopen("mesh/mesh2_triangles.txt","rb");
        nT = 0;
        while (! feof(f)) {
            size_t i1,i2,i3;
            fscanf(f, "%ld %ld %ld",&i1,&i2,&i3);
            T.push_back(i1);
            T.push_back(i2);
            T.push_back(i3);
            nT++;
        }
        fclose(f);
    }
    printf("Triangles: %ld\n", T.size());


    vec x;
    vec rhs;
    rhs.resize(2*nP);
    x.resize(2*nP);
    for (const size_t i : B1) rhs[i*2+0] = 1;
    for (const size_t i : B2) rhs[i*2+1] = 1;
    std::function<void(const vec&, vec&)> mult = [](const vec& x, vec& Mx){
        MMult(x.data(), Mx.data());
    };
    Solve(mult, rhs, x, 10000);
    
    write_vtu("out.vtu", P, T, std::vector{std::make_tuple(std::string("Eta"),2,&x)});
    return 0;
}
