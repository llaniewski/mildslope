#include "global.h"
#include "vtu_write.h"
#include "solve.h"

#include <math.h>

#include <Eigen/Sparse>
#include <Eigen/KLUSupport>
#include <vector>
#include <iostream>
#include <set>


size_t * triangles;
size_t n_triangles;
size_t n_points;
size_t * boundary;
int * boundary_flag;
size_t n_boundary;
extern "C" {
    void problem(double wave_k, const double *points, const double* x, double* res, double* obj);
    void problem_d(double wave_k, const double *points, const double *pointsd, const double *x, const double *xd, double *res, double *resd, double *obj, double *objd);
    //void problem_d(double k, const double *points, const double *x, const double *xd, double *res, double *resd, double *obj);
}

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Trip;

int main() {

    std::vector<double> P;
    size_t nP;
    std::vector<size_t> T;
    size_t nT;
    std::vector<size_t> B;
    std::vector<int> B_flag;
    size_t nB;

    std::string mesh = "mesh/empty2";

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
        auto fun = [&nT,&T,&nP,&P,&nB,&B,&B_flag](double val, int flag){
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
        fun(minx, 1);
        fun(maxx, 2);
    }
    nB = B.size()/2;
    printf("Borders: %ld\n", nB);



    const size_t DOFperP = 2;
    const size_t DOF = nP*DOFperP;


    // Graph coloring
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


    // Pass mesh data as globals to C
    triangles = T.data();
    n_triangles = nT;
    n_points = nP;
    boundary = B.data();
    boundary_flag = B_flag.data();
    n_boundary = nB;

    // problem coefficient
    double wave_k = 4.0;

    Eigen::VectorXd x(DOF);
    Eigen::VectorXd res(DOF);
    Eigen::VectorXd obj(1);

    problem(wave_k, P.data(), x.data(), res.data(), obj.data());
    double resL2 = res.norm();
    printf("Residual: %lg\n", resL2);
    
    Eigen::VectorXd Pd(DOF);
    Eigen::VectorXd res_tmp(DOF);
    Eigen::VectorXd obj_tmp(1);
    Eigen::VectorXd objd_tmp(1);

    Eigen::VectorXd Mx(DOF);

    printf("Gathering Jacobian...\n");
    std::vector<Trip> coef;
    for (int k=0; k<maxk; k++) {
        problem_d(wave_k, P.data(), Pd.data(), x.data(), ref_x.col(k).data(), res_tmp.data(), Mx.data(), obj_tmp.data(), objd_tmp.data());
        for (size_t j=0; j<DOF; j++){
            if (fabs(Mx[j]) > 1e-6) {
                coef.push_back(Trip(j,ref_j(j,k),Mx[j]));
            }
        }
        printf("mult %ld -> %ld\n", k, coef.size());
    }
    SpMat A(DOF,DOF);

    printf("Constructing Jacobian from triplets\n");
    A.setFromTriplets(coef.begin(), coef.end());
 
    printf("Solving linear problem");
    Eigen::KLU<SpMat> solver;  // performs a Cholesky factorization of A
    printf(" [compute]");
    solver.compute(A); assert(solver.info() == Eigen::Success);
    printf(" [solve]");
    Eigen::VectorXd ret = solver.solve(res); assert(solver.info() == Eigen::Success);
    printf(" [done]\n");

    for (size_t i=0; i<ret.size(); i++) x[i] -= ret[i];

    problem(wave_k, P.data(), x.data(), res.data(), obj.data());

    write_vtu("out.vtu", P, T, {
        std::make_tuple(std::string("Eta"),2,std::span(x.data(),x.size())),
        std::make_tuple(std::string("res"),2,std::span(res.data(),res.size()))
    });
    return 0;
}
