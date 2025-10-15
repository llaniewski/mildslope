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
#include <nlopt.h>

const double pi = 4.0 * atan(1.0);

size_t * triangles;
size_t n_triangles;
size_t n_points;
size_t * boundary;
int * boundary_flag;
size_t n_boundary;
extern "C" {
    void problem(double wave_k, const double *points, const double* x, double* res, double* obj);
    void problem_d(double wave_k, const double *points, const double *x, const double *xd, double *res, double *resd, double *obj);
    void problem_bP(double wave_k, const double *points, double *pointsb, const double *x, double *res, double *resb, double *obj, double *objb);
    void problem_bX(double wave_k, const double *points, const double *x, double *xb, double *res, double *obj, double *objb);
}

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Trip;

template <class T>
std::span<T> to_span(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec) { return std::span(vec.data(), vec.size()); }

int main() {

    size_t nP;
    std::vector<size_t> T;
    size_t nT;
    std::vector<size_t> B;
    std::vector<int> B_flag;
    size_t nB;

    std::string mesh = "mesh/mesh2";
    Eigen::Matrix<double, 2, Eigen::Dynamic> P;
    {
        std::vector<double> Pv;
        FILE* f = fopen((mesh+"_points.txt").c_str(),"rb");
        nP = 0;
        while (! feof(f)) {
            double x,y;
            fscanf(f, "%lf %lf",&x,&y);
            if (feof(f)) break;
            Pv.push_back(x);
            Pv.push_back(y);
            nP++;
        }
        fclose(f);
        P.resize(2, nP);
        for (size_t i=0;i<nP;i++) {
            P(0,i) = Pv[2*i+0];
            P(1,i) = Pv[2*i+1];
        }
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
        double maxx=P(0,0), minx=P(0,0);
        for (size_t i=0;i<nP;i++) {
            if (P(0,i) > maxx) maxx = P(0,i);
            if (P(0,i) < minx) minx = P(0,i);
        }
        printf("X: [%lg, %lg]\n", minx, maxx);
        auto fun = [&nT,&T,&nP,&P,&nB,&B,&B_flag](double val, int flag){
            for (size_t i=0;i<nT;i++) {
                int count=0;
                for (int j=0; j<3; j++) {
                    size_t idx = T[i*3+j];
                    if (fabs(P(0,idx) - val) < 1e-6) {
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

    const size_t NPAR_SIDE = 60;
    const size_t NPAR_PER_NODE = 1;
    const size_t NPAR = NPAR_SIDE*NPAR_PER_NODE;
    Eigen::Matrix<double, Eigen::Dynamic, NPAR> par(DOF,NPAR); par.setZero();
    if (true) {
        std::vector<double> nodes_x;
        for (size_t i=0; i<NPAR_SIDE+2; i++) nodes_x.push_back(3.0+(7.0-3.0)*(i)/(NPAR_SIDE+1));
        const auto& fun = [&nodes_x](size_t i, int side, double x, double y) {
            double x0 = nodes_x[i+0];
            double x1 = nodes_x[i+1];
            double x2 = nodes_x[i+2];
            double ret = 1;
            if (x < x0) ret = 0;
            else if (x < x1) ret = (x-x0)/(x1-x0);
            else if (x < x2) ret = (x-x2)/(x1-x2);
            else ret = 0;
            double r = sqrt((5-x)*(5-x) + (0.5-y)*(0.5-y));
            if (r < 0.3) ret = 0;
            else if (r < 0.5) ret = ret * (r-0.3)/(0.5-0.3);
            else ret = ret;
            if (side) ret = ret * (y-1.0) / (0.0-1.0);
            else ret = ret * (y-0.0) / (1.0-0.0);
            return ret;
        };
        for (int j = 0; j<NPAR_SIDE; j++) {
            for (int i = 0; i<nP; i++) {
                par(i*2+1,j*NPAR_PER_NODE+0) = - fun(j,true,P(0,i),P(1,i)) + fun(j,false,P(0,i),P(1,i));
                // par(i*2+0,j*NPAR_PER_NODE+0) = fun(j,true,P(0,i),P(1,i));
                // par(i*2+1,j*NPAR_PER_NODE+1) = -fun(j,true,P(0,i),P(1,i));
                // par(i*2+0,j*NPAR_PER_NODE+2) = fun(j,false,P(0,i),P(1,i));
                // par(i*2+1,j*NPAR_PER_NODE+3) = fun(j,false,P(0,i),P(1,i));
            }
        }
    }

    {
        std::vector<std::tuple<std::string, int, std::span<double> > >fields;
        for (int j = 0; j<NPAR; j++) {
            char buf[1024];
            sprintf(buf, "par_%02d", j);
            fields.push_back(std::make_tuple(
                std::string(buf),
                2,
                std::span(par.col(j).data(),par.col(j).size())
            ));
        }
        write_vtu("output/par.vtu", std::span(P.data(), P.size()), T, fields);
    }

    // Graph coloring
    const size_t W = 40;
    const size_t NA = -1;
    int maxk = 0;
    Eigen::Matrix<size_t, Eigen::Dynamic, W> ref_j(DOF,W); ref_j.setZero();
    Eigen::Matrix<double, Eigen::Dynamic, W> ref_x(DOF,W); ref_x.setZero();
    {
        Eigen::Matrix<bool, Eigen::Dynamic, W> ref_b(DOF,W); ref_b.setZero();
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
    //double wave_k = 4.0;
    std::vector<std::pair<double, double>> k_integral;
    const size_t KINT = 100;
    const double wave_k_min = 0.01 * pi;
    const double wave_k_max = 1 * pi;
    {
        double k_dist = (wave_k_max-wave_k_min)/(KINT - 1);
        for (size_t i=0; i<KINT; i++) {
            double weight = k_dist;
            if ((i == 0) || (i == KINT-1)) weight = weight/2;
            double wave_k = wave_k_min + i*k_dist;
            k_integral.push_back(std::make_pair(weight, wave_k));
        }
    }
    Eigen::Matrix<double, 2, Eigen::Dynamic> P0 = P;

    int iter = 0;
    const auto& objective = [&](const double *pr_, double* grad_, bool export_all=false) -> double {
        Eigen::Map< const Eigen::VectorXd > pr(pr_, NPAR);
        P = P0 + (par * pr).reshaped(2,nP);
        double total_obj = 0;
        Eigen::Map< Eigen::VectorXd >total_grad(grad_, NPAR);
        if (grad_ != NULL) total_grad.setZero();
        std::vector<std::pair<double, double> > objs; 
        for(size_t m = 0; m < KINT; m++) {
            double weight = k_integral[m].first;
            double wave_k = k_integral[m].second;
            Eigen::VectorXd x(DOF); x.setZero();
            Eigen::VectorXd res(DOF);
            Eigen::VectorXd obj(1);

            problem(wave_k, P.data(), x.data(), res.data(), obj.data());
            //printf("obj:%lg\n", obj[0]);

            {
                double resL2 = res.norm();
                printf("Residual (before): %lg\n", resL2);
            }
            
            Eigen::VectorXd Pd(DOF);
            Eigen::VectorXd res_tmp(DOF);
            Eigen::VectorXd obj_tmp(1);

            Eigen::VectorXd Mx(DOF);

            SpMat A(DOF,DOF);
            {
                printf("Gathering Jacobian");
                std::vector<Trip> coef;
                printf(" [mult]");
                for (size_t k=0; k<maxk; k++) {
                    problem_d(wave_k, P.data(), x.data(), ref_x.col(k).data(), res_tmp.data(), Mx.data(), obj_tmp.data());
                    for (size_t j=0; j<DOF; j++){
                        if (fabs(Mx[j]) > 1e-6) {
                            coef.push_back(Trip(j,ref_j(j,k),Mx[j]));
                        }
                    }
                    //printf("mult %ld -> %ld\n", k, coef.size());
                }
                printf(" [sparse]");
                A.setFromTriplets(coef.begin(), coef.end());
                printf(" [done]\n");
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
            {
                double resL2 = res.norm();
                printf("Residual (after): %lg\n", resL2);
            }
            printf("obj:%lg\n", obj[0]);
            total_obj += obj[0]*weight;
            objs.push_back(std::make_pair(wave_k, obj[0]));

            // ADJOINT
            if (grad_ != NULL) {
                Eigen::VectorXd Pb(DOF); Pb.setZero();
                Eigen::VectorXd xb(DOF); xb.setZero();
                Eigen::VectorXd objb(1); objb.setZero();
                objb[0] = 1;
                problem_bX(wave_k, P.data(), x.data(), xb.data(), res_tmp.data(), obj_tmp.data(), objb.data());

                Eigen::VectorXd resb;
                {   
                    printf("Solving adjoint problem");
                    Eigen::KLU<SpMat> solver;  // performs a Cholesky factorization of A
                    printf(" [compute]");
                    solver.compute(A.transpose()); assert(solver.info() == Eigen::Success);
                    printf(" [solve]");
                    resb = solver.solve(-xb); assert(solver.info() == Eigen::Success);
                    printf(" [done]\n");
                }
                
                problem_bP(wave_k, P.data(), Pb.data(), x.data(), res_tmp.data(), resb.data(), obj_tmp.data(), objb.data());

                Eigen::VectorXd grad = Pb.transpose() * par;
                total_grad += grad * weight;
            }
            if (export_all || (m == KINT-1)) {
                char buf[1024];
                sprintf(buf, "output/res_%lg_%04d.vtu", wave_k, iter);
                write_vtu(buf, std::span(P.data(), P.size()), T, {
                    std::make_tuple(std::string("Eta"), 2, to_span(x))
                });
            }
        }
        {
            char buf[1024];
            sprintf(buf, "output/res_%04d.csv", iter);
            FILE* f = fopen(buf, "w");
            fprintf(f, "wave_k,obj\n");
            for (size_t i=0;i<objs.size();i++) {
                fprintf(f, "%.15lg, %.15lg\n", objs[i].first, objs[i].second);
            }
            fclose(f);
        }
        iter++;
        return total_obj;
    };
    
    if (false) { // FD test
        Eigen::VectorXd pr(NPAR); pr.setZero();
        Eigen::VectorXd gr(NPAR);
        pr(2) += 0.05;
        pr(5) += -0.05;
        double val = objective(pr.data(), gr.data());
        double h = 1e-4;
        for (int i=0; i<NPAR; i++) {
            pr(i) += h;
            double val1 = objective(pr.data(), NULL);
            pr(i) -= 2*h;
            double val2 = objective(pr.data(), NULL);
            pr(i) += h;
            double fd = (val1-val2)/(2*h);
            printf("grad: %d %lg -- %lg => %lg\n", i, fd, gr(i), fd - gr(i));
        }
        return 0;
    }

    using obj_type = decltype(&objective);
    nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, NPAR);
    nlopt_result opt_res;
    opt_res = nlopt_set_min_objective(opt, 
        [](unsigned n, const double* x, double* grad, void* f_data) -> double {
            obj_type fun = (obj_type) f_data;
            return (*fun)(x, grad);
        }, (void*) &objective);
    Eigen::VectorXd lower(NPAR);
    Eigen::VectorXd upper(NPAR);
    for (size_t i=0;i<NPAR_SIDE; i++) {
        lower(i*NPAR_PER_NODE+0) = -0.2;
        upper(i*NPAR_PER_NODE+0) =  1.0;
        // lower(i*NPAR_PER_NODE+0) = -0.01;
    //     upper(i*NPAR_PER_NODE+0) = 0.01;
    //     lower(i*NPAR_PER_NODE+1) = -0.2;
    //     upper(i*NPAR_PER_NODE+1) = 1;
    //     lower(i*NPAR_PER_NODE+2) = -0.01;
    //     upper(i*NPAR_PER_NODE+2) = 0.01;
    //     lower(i*NPAR_PER_NODE+3) = -0.2;
    //     upper(i*NPAR_PER_NODE+3) = 1;
    }
    opt_res = nlopt_set_lower_bounds(opt, lower.data());
    opt_res = nlopt_set_upper_bounds(opt, upper.data());
    opt_res = nlopt_set_maxeval(opt, 500);

    Eigen::VectorXd pr(NPAR); pr.setZero();
    double obj;
    opt_res = nlopt_optimize(opt, pr.data(), &obj);
    std::cout << pr << "\n";
    printf("Objective: %lg\n", obj);
    objective(pr.data(),NULL,true);

    // for (int k=0; k<60; k++) {
    //     double val = objective(pr.data(), gr.data());
    //     pr += -1e-3*gr;
    // }
    return 0;
}
