#include "global.h"
#include "vtu_write.h"
#include "solve.h"

#include <math.h>

#include <Eigen/Sparse>
#include <Eigen/KLUSupport>
#include <Eigen/LU>
#include <vector>
#include <iostream>
#include <set>
#include <fstream>

size_t * triangles;
size_t n_triangles;
size_t n_points;
size_t * boundary;
int * boundary_flag;
size_t n_boundary;
extern "C" {
    void problem(double wave_k, const double *points, const double* x, double* res, double* obj);
    void problem_d(double wave_k, const double *points, const double *pointsd, const double *x, const double *xd, double *res, double *resd, double *obj);
    void problem_dP(double wave_k, const double *points, const double *pointsd, const double *x, double *res, double *resd, double *obj);
    void problem_dX(double wave_k, const double *points, const double *x, const double *xd, double *res, double *resd, double *obj);
    void problem_b(double wave_k, const double *points, double *pointsb, const double *x, double *xb, double *res, double *obj, const double *objb);
}

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Trip;

template <class T>
std::span<T> to_span(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec) { return std::span(vec.data(), vec.size()); }

int main() {

    std::vector<double> P;
    size_t nP;
    std::vector<size_t> T;
    size_t nT;
    std::vector<size_t> B;
    std::vector<int> B_flag;
    size_t nB;

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

    const size_t NPAR_SIDE = 10;
    const size_t NPAR = 2*NPAR_SIDE;
    Eigen::Matrix<double, Eigen::Dynamic, NPAR> par(nP,NPAR);
    if (true) {
        const auto& fun = [](double x1, double y1, double x2, double y2) {
            double ret = 1;
            if (x2 < 3) ret = 0;
            else if (x2 < x1) ret = (x2 - 3)/(x1 - 3);
            else if (x2 < 7) ret = (x2 - 7)/(x1 - 7);
            else ret = 0;
            if (y1 == 0.0) ret = ret * (1-y2) / (1-y1);
            else if (y1 == 1.0) ret = ret * (0 - y2) / (0-y1);
            else exit(3);
            double r = sqrt((5-x2)*(5-x2) + (0.5-y2)*(0.5-y2));
            if (r < 0.3) ret = 0;
            else if (r < 0.5) ret = ret * (r-0.3)/(0.5-0.3);
            else ret = ret;
            return ret;
        };
        Eigen::Matrix<double, NPAR, 2> par_points(NPAR,2);
        for (int i = 0; i<NPAR_SIDE; i++) {
            par_points(i,0) = 3.0+(7.0-3.0)*(i+1)/(NPAR_SIDE+1);
            par_points(i,1) = 0.0;
            par_points(i+NPAR_SIDE,0) = 3.0+(7.0-3.0)*(i+1)/(NPAR_SIDE+1);
            par_points(i+NPAR_SIDE,1) = 1.0;
        }
        std::cout << par_points << '\n';
        Eigen::Matrix<double, NPAR, NPAR> par_cross(NPAR, NPAR);
        for (int i = 0; i<NPAR; i++) {
            for (int j = 0; j<NPAR; j++) {
                par_cross(i,j) = fun(par_points(j,0),par_points(j,1),par_points(i,0),par_points(i,1));
            }
        }
        // std::cout << par_cross << '\n';
        // {
        //     std::ofstream file("mat.txt");
        //     file << par_cross << '\n';
        // }
        // {
        //     std::ofstream file("inv.txt");
        //     file << par_cross.inverse() << '\n';
        // }
        Eigen::Matrix<double, Eigen::Dynamic, NPAR> par_a(nP, NPAR);
        for (int i = 0; i<nP; i++) {
            for (int j = 0; j<NPAR; j++) {
                par_a(i,j) = fun(par_points(j,0),par_points(j,1),P[2*i+0],P[2*i+1]);
            }
        }
        
        par = par_a * par_cross.inverse();

        // Eigen::Matrix<double, Eigen::Dynamic, NPAR> par_test = par_cross * par_cross.inverse();
        // std::cout << par_test << '\n';
    }

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
    printf("obj:%lg\n", obj[0]);

    double resL2 = res.norm();
    printf("Residual: %lg\n", resL2);
    
    Eigen::VectorXd Pd(DOF);
    Eigen::VectorXd res_tmp(DOF);
    Eigen::VectorXd obj_tmp(1);
    Eigen::VectorXd objd_tmp(1);

    Eigen::VectorXd Mx(DOF);

    SpMat A(DOF,DOF);
    {
        printf("Gathering Jacobian...\n");
        std::vector<Trip> coef;
        for (int k=0; k<maxk; k++) {
            problem_dX(wave_k, P.data(), x.data(), ref_x.col(k).data(), res_tmp.data(), Mx.data(), obj_tmp.data());
            for (size_t j=0; j<DOF; j++){
                if (fabs(Mx[j]) > 1e-6) {
                    coef.push_back(Trip(j,ref_j(j,k),Mx[j]));
                }
            }
            printf("mult %ld -> %ld\n", k, coef.size());
        }
        printf("Constructing Jacobian from triplets\n");
        A.setFromTriplets(coef.begin(), coef.end());
    }

    {
        printf("Solving linear problem");
        Eigen::KLU<SpMat> solver;  // performs a Cholesky factorization of A
        printf(" [compute]");
        solver.compute(A); assert(solver.info() == Eigen::Success);
        printf(" [solve]");
        Eigen::VectorXd ret = solver.solve(res); assert(solver.info() == Eigen::Success);
        printf(" [done]\n");
        for (size_t i=0; i<ret.size(); i++) x[i] -= ret[i];
    }

    

    problem(wave_k, P.data(), x.data(), res.data(), obj.data());
    printf("obj:%lg\n", obj[0]);
    // ADJOINT

    Eigen::VectorXd Pb(DOF);
    Eigen::VectorXd xb(DOF);
    Eigen::VectorXd objb(1);
    objb[0] = 1;
    problem_b(wave_k, P.data(), Pb.data(), x.data(), xb.data(), res_tmp.data(), obj_tmp.data(), objb.data());

    Eigen::VectorXd x_adj;
    {   
        printf("Solving adjoint problem");
        Eigen::KLU<SpMat> solver;  // performs a Cholesky factorization of A
        printf(" [compute]");
        solver.compute(A.transpose()); assert(solver.info() == Eigen::Success);
        printf(" [solve]");
        x_adj = solver.solve(xb); assert(solver.info() == Eigen::Success);
        printf(" [done]\n");
    }
    
    SpMat dRdP(DOF,DOF);
    {
        printf("Gathering dRdP...\n");
        std::vector<Trip> coef;
        for (int k=0; k<maxk; k++) {
            problem_dP(wave_k, P.data(), ref_x.col(k).data(), x.data(), res_tmp.data(), Mx.data(), obj_tmp.data());
            for (size_t j=0; j<DOF; j++){
                if (fabs(Mx[j]) > 1e-6) {
                    coef.push_back(Trip(j,ref_j(j,k),Mx[j]));
                }
            }
            printf("mult %ld -> %ld\n", k, coef.size());
        }
        printf("Constructing B from triplets\n");
        dRdP.setFromTriplets(coef.begin(), coef.end());
    }
    
    Eigen::VectorXd p_adj = Pb - dRdP.transpose() * x_adj;

    Eigen::MatrixXd grad = p_adj.reshaped(2,nP) * par;
    std::cout << grad << "\n";
    
    write_vtu("out.vtu", P, T, {
        std::make_tuple(std::string("Eta"), 2, to_span(x)),
        std::make_tuple(std::string("Eta_adj"), 2, to_span(x_adj)),
        std::make_tuple(std::string("Pb"), 2, to_span(Pb)),
        std::make_tuple(std::string("grad"), 2, to_span(p_adj)),
        std::make_tuple(std::string("res"), 2, to_span(res)),
        std::make_tuple(std::string("par"), 1, std::span(par.col(0).data(),par.col(0).size()))
    });
    return 0;
}
