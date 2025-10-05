#include "global.h"
#include "vtu_write.h"
#include "solve.h"

#include <math.h>


#include <Eigen/Sparse>
#include <Eigen/KLUSupport>
#include <vector>
#include <iostream>
#include <set>

std::vector<double> P;
size_t nP;
std::vector<size_t> T;
size_t nT;
std::vector<size_t> B;
std::vector<int> B_flag;
size_t nB;

double k = 2;

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
        double a = -area*k*k;
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
    for (size_t i=0;i<nB;i++) {
        size_t i0 = B[i*2+0];
        size_t i1 = B[i*2+1];
        double x0 = P[i0*2+0];
        double y0 = P[i0*2+1];
        double x1 = P[i1*2+0];
        double y1 = P[i1*2+1];
        double vx = x1-x0;
        double vy = y1-y0;
        double len = sqrt(vx*vx+vy*vy);
        double b = len;
        double robin = k;
        double EM00, EM01, EM11;
        EM00 = EM11 = 1.0/3.0;
        EM01 = 1.0/6.0;
        if (B_flag[i] == 1) {
            // Mx[i0*2+0] =  x[i0*2+0];
            // Mx[i0*2+1] =  x[i0*2+1];
            // Mx[i1*2+0] =  x[i1*2+0];
            // Mx[i1*2+1] =  x[i1*2+1];
            Mx[i0*2+0] +=  b*robin*(EM00 * x[i0*2+1] + EM01 * x[i1*2+1]);
            Mx[i1*2+0] +=  b*robin*(EM01 * x[i0*2+1] + EM11 * x[i1*2+1]);
            Mx[i0*2+1] += -b*robin*(EM00 * x[i0*2+0] + EM01 * x[i1*2+0]);
            Mx[i1*2+1] += -b*robin*(EM01 * x[i0*2+0] + EM11 * x[i1*2+0]);
        } else if (B_flag[i] == 2) {
            Mx[i0*2+0] +=  b*robin*(EM00 * x[i0*2+1] + EM01 * x[i1*2+1]);
            Mx[i1*2+0] +=  b*robin*(EM01 * x[i0*2+1] + EM11 * x[i1*2+1]);
            Mx[i0*2+1] += -b*robin*(EM00 * x[i0*2+0] + EM01 * x[i1*2+0]);
            Mx[i1*2+1] += -b*robin*(EM01 * x[i0*2+0] + EM11 * x[i1*2+0]);
        }
    }
}



int main() {
    std::string mesh = "mesh/mesh2";
    {
        FILE* f = fopen((mesh+"_points.txt").c_str(),"rb");
        nP = 0;
        while (! feof(f)) {
            double x,y;
            fscanf(f, "%lf %lf",&x,&y);
            if (feof(f)) break;
            P.push_back(x);
            P.push_back(y);
            nP++;
        }
        fclose(f);
    }
    printf("Points %ld\n", nP);
    {
        FILE* f = fopen((mesh+"_triangles.txt").c_str(),"rb");
        nT = 0;
        while (! feof(f)) {
            size_t i1,i2,i3;
            fscanf(f, "%ld %ld %ld",&i1,&i2,&i3);
            if (feof(f)) break;
            T.push_back(i1);
            T.push_back(i2);
            T.push_back(i3);
            nT++;
        }
        fclose(f);
    }
    printf("Triangles: %ld\n", nT);

    {
        double maxx=P[0], minx=P[0];
        for (size_t i=0;i<nP;i++) {
            if (P[i*2+0] > maxx) maxx = P[i*2+0];
            if (P[i*2+0] < minx) minx = P[i*2+0];
        }
        printf("X: [%lg, %lg]\n", minx, maxx);
        auto fun = [](double val, std::vector<size_t>& B, int flag){
            for (size_t i=0;i<nT;i++) {
                int count=0;
                for (int j=0; j<3; j++) {
                    size_t idx = T[i*3+j];
                    if (fabs(P[idx*2+0] - val) < 1e-6) {
                        //printf("Added: %ld %ld\n",i, idx);
                        B.push_back(idx);
                        count++;
                    }
                }
                if (count == 1) {
                    B.pop_back();
                } else if (count == 2) {
                    B_flag.push_back(flag);
                } else if (count > 2) {
                    fprintf(stderr, "Whole triangle on edge\n");
                    exit(2);
                }
            }
            if (B.size() % 2 != 0) {
                fprintf(stderr, "border edges wrong\n");
                exit(2);
            }
        };
        fun(minx, B, 1);
        fun(maxx, B, 2);
    }
    nB = B.size()/2;
    printf("Borders: %ld\n", nB);


    vec x;
    vec rhs;
    rhs.resize(2*nP);
    x.resize(2*nP);
    for (size_t i=0; i<x.size(); i++) x[i] = 0;
    for (size_t i=0; i<nB; i++) {
        if (B_flag[i] == 1) {
            size_t i0 = B[i*2+0];
            size_t i1 = B[i*2+1];
            double x0 = P[i0*2+0];
            double y0 = P[i0*2+1];
            double x1 = P[i1*2+0];
            double y1 = P[i1*2+1];
            double vx = x1-x0;
            double vy = y1-y0;
            double len = sqrt(vx*vx+vy*vy);
            double b = 2*len;
            rhs[i0*2+0] += b;
            rhs[i0*2+1] += 0;
            rhs[i1*2+0] += b;
            rhs[i1*2+1] += 0;
        }
    }
    std::function<void(const vec&, vec&)> mult = [](const vec& x, vec& Mx){
        MMult(x.data(), Mx.data());
    };


    typedef Eigen::SparseMatrix<double> SpMat;  // declares a column-major sparse matrix type of double
    typedef Eigen::Triplet<double> Trip;

    const size_t DOFperP = 2;
    const size_t DOF = nP*DOFperP;

    const size_t W = 40;
    const size_t NA = -1;
    int maxk = 0;
    Eigen::Matrix<size_t, Eigen::Dynamic, W> ref_j(DOF,W);
    Eigen::Matrix<double, Eigen::Dynamic, W> ref_x(DOF,W);
    {
        Eigen::Matrix<bool, Eigen::Dynamic, W> ref_b(DOF,W);
        printf("construct graph\n");
        std::vector<std::set<size_t>> graph(DOF);
        for (size_t i=0;i<nT;i++) {
            size_t Ti[3];
            for (char j=0;j<3;j++) Ti[j] = T[i*3+j];
            for (char j=0;j<3;j++)
                for (char k=0;k<3;k++)
                    for (char d1=0;d1<DOFperP;d1++)
                        for (char d2=0;d2<DOFperP;d2++)
                            graph[Ti[j]*DOFperP+d1].insert(Ti[k]*DOFperP+d2);
        }
        printf("done\n");
        for (size_t i=0;i<DOF;i++)
            for (int k=0;k<W;k++)
                ref_j(i,k) = NA;
        for (size_t i=0;i<DOF;i++) {
            //printf("i:%ld\n",i);
            int k;
            for (k=0;k<W;k++) {
                if (!ref_b(i,k)) {
                    ref_x(i,k) = 1.0;
                    //printf("i:%ld k:%ld\n",i,k);
                    for (size_t j : graph[i]) {
                        //printf("i:%ld k:%ld j:%ld\n",i,k,j);
                        if (ref_j(j,k) == NA) {
                            ref_j(j,k) = i;
                            for (size_t p : graph[j]) {
                                ref_b(p,k) = true;
                            }
                        } else {
                            printf("wrong\n");
                            exit(2);
                        }
                    }
                    break;
                }
            }
            if (k >= maxk) maxk = k + 1;
            //printf("i:%ld\n",i);
        }   
        printf("maxk: %d\n",maxk);
        if (maxk > W) {
            printf("Not enought colors\n");
            exit(3);
        }
    }


    assert(x.size() == DOF);

    std::vector<Trip> coef;
    vec Mx; Mx.resize(x.size());
    // for (size_t i=0; i<x.size(); i++) {
    //     x[i] = 1;
    //     MMult(x.data(), Mx.data());
    //     x[i] = 0;
    //     for (size_t j=0; j<x.size(); j++){
    //         if (fabs(Mx[j]) > 1e-6) {
    //             coef.push_back(Trip(j,i,Mx[j]));
    //         }
    //     }
    //     printf("mult %ld -> %ld\n", i, coef.size());
    // }
    // coef.clear();

    for (int k=0; k<maxk; k++) {
        MMult(ref_x.col(k).data(), Mx.data());
        for (size_t j=0; j<DOF; j++){
            if (fabs(Mx[j]) > 1e-6) {
                coef.push_back(Trip(j,ref_j(j,k),Mx[j]));
            }
        }
       printf("mult %ld -> %ld\n", k, coef.size());
    }


    //    fclose(outf);
    printf("SpMat\n");
    SpMat A(x.size(), x.size());
    printf("setFromTriplets\n");
    A.setFromTriplets(coef.begin(), coef.end());
 
    printf("KLU\n");
    Eigen::KLU<SpMat> solver;  // performs a Cholesky factorization of A
    printf("compute\n");
    solver.compute(A);
    assert(solver.info() == Eigen::Success);
    printf("solve\n");
    Eigen::VectorXd b(rhs.size());
    for (size_t i=0; i<rhs.size(); i++) b[i] = rhs[i];
    Eigen::VectorXd ret = solver.solve(b);         // use the factorization to solve for the given right hand side
    for (size_t i=0; i<ret.size(); i++) x[i] = ret[i];
    Eigen::VectorXd ret_m = A * ret;
    vec Ax; Ax.resize(ret_m.size());
    for (size_t i=0; i<ret_m.size(); i++) Ax[i] = ret_m[i];

    write_vtu("out.vtu", P, T, {
        std::make_tuple(std::string("Eta"),2,&x),
        std::make_tuple(std::string("RHS"),2,&rhs),
        std::make_tuple(std::string("Ax"),2,&Ax)
    });
    return 0;
}
