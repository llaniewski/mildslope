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
size_t n_bord;
size_t* bord;
double* fix;
double* fix_dirs;

extern "C" {
    void problem(double wave_k, const double *points, const double *depth, const double* x, double* res, double* obj);
    void problem_d(double wave_k, const double *points, const double *depth, const double *x, const double *xd, double *res, double *resd, double *obj);
    void problem_bP(double wave_k, const double *points, double *pointsb, const double *depth, const double *depthb, const double *x, double *res, double *resb, double *obj, double *objb);
    void problem_bX(double wave_k, const double *points, const double *depth, const double *x, double *xb, double *res, double *obj, double *objb);
    void morph_energy_fix(const double *P0, const double *P1, const double *Pfix, double *res, double *energy, double *energyb);
    void morph_energy_fix_d(const double *P0, const double *P1, const double *P1d, const double *Pfix, double *res, double *resd, double *energyb);
    void morph_energy_fix_b(const double *P0, const double *P1, const double *Pfix, double *Pfixb, double *res, double *resb, double *energyb);
}

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Trip;

// template <class T>
// std::span<T> to_span(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec) { return std::span(vec.data(), vec.size()); }
template <typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
std::span<Scalar> to_span(Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>& vec) { return std::span(vec.data(), vec.size()); }
template <typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
std::span<Scalar> to_span(Eigen::Array<Scalar, RowsAtCompileTime, ColsAtCompileTime>& vec) { return std::span(vec.data(), vec.size()); }

FILE* fopen_safe(std::string fn, char* mode) {
    FILE* f = fopen(fn.c_str(), mode);
    if (f == NULL) {
        fprintf(stderr, "Cannot open file: %s\n", fn.c_str());
        exit(2);
    }
    return f;
}

class mesh {
public:
    size_t nP;
    size_t nT;
    size_t nB;
    size_t nAttr;
    Eigen::Matrix<double, 2, Eigen::Dynamic> P;
    Eigen::Array<size_t, 3, Eigen::Dynamic> T;
    Eigen::Array<size_t, 2, Eigen::Dynamic> B;
    Eigen::Array<int, Eigen::Dynamic, 1> B_flag;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Attr;
    mesh() {
        nP = 0;
        nB = 0;
        nT = 0;
        nAttr = 0;        
    }
    int read_mesh(const std::string& mesh_name) {
        FILE* f;
        int dim, bord, el_size, ignore;
        f = fopen_safe(mesh_name+".node","rb");
        fscanf(f, "%ld %d %ld %d", &nP, &dim, &nAttr, &bord);
        assert(dim == 2);
        P.resize(dim, nP);
        Attr.resize(nP, nAttr);
        for (size_t i=0;i<nP;i++) {
            size_t idx;
            fscanf(f, "%ld", &idx);
            assert(idx == i+1);
            for (int j=0;j<dim;j++) fscanf(f,"%lf", &P(j,i));
            for (size_t j=0;j<nAttr;j++) fscanf(f,"%lf", &Attr(i,j));
            for (int j=0;j<bord;j++) fscanf(f,"%d", &ignore);
        }
        fclose(f);
        printf("Points %ld\n", nP);
        f = fopen_safe(mesh_name+".ele","rb");
        fscanf(f, "%ld %d %d", &nT, &el_size, &bord);
        assert(el_size == 3);
        T.resize(el_size, nT);
        for (size_t i=0;i<nT;i++) {
            size_t idx;
            fscanf(f, "%ld", &idx);
            assert(idx == i+1);
            for (int j=0;j<el_size;j++) { size_t tmp; fscanf(f,"%ld", &tmp); T(j,i) = tmp-1; }
            for (int j=0;j<bord;j++) fscanf(f,"%d", &ignore);
        }
        fclose(f);
        printf("Triangles: %ld\n", nT);
        f = fopen_safe(mesh_name+".poly","rb");
        {
            size_t nP_, nAttr_;
            int dim_, bord_;
            fscanf(f, "%ld %d %ld %d", &nP_, &dim_, &nAttr_, &bord_);
            assert(nP_ == 0);
        }
        const int edge_size = 2;
        fscanf(f, "%ld %d", &nB, &bord);        
        assert(bord == 1);
        B.resize(edge_size, nB);
        B_flag.resize(nB, bord);
        for (size_t i=0;i<nB;i++) {
            size_t idx;
            fscanf(f, "%ld", &idx);
            assert(idx == i+1);
            for (int j=0;j<edge_size;j++) { size_t tmp; fscanf(f,"%ld", &tmp); B(j,i) = tmp-1; }
            for (int j=0;j<bord;j++) fscanf(f,"%d", &B_flag(i,j));
        }
        fclose(f);
        printf("Border: %ld\n", nB);
        return 0;
    }
};

int main(int argc, char **argv) {
    std::vector<std::string> args;
    for (size_t i = 0;i<argc;i++) args.push_back(argv[i]);

    if (args.size() != 2) {
        fprintf(stderr, "Usage: ./main case\n");
        exit(2);
    }
    std::string mesh_name = args[1];

    mesh m;

    m.read_mesh(mesh_name);

    const size_t DOFperP = 2;
    const size_t DOF = m.nP*DOFperP;

    Eigen::VectorXd Pfix(DOF);


    // size_t NPAR_SIDE = nAttr;
    // const size_t NPAR_PER_NODE = 1;
    // const size_t NPAR_SHAPE = NPAR_SIDE*NPAR_PER_NODE;
    const size_t NPAR_SHAPE = 0;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> par_shape(DOF,NPAR_SHAPE); par_shape.setZero();
    // {
    //     for (int j = 0; j<NPAR_SIDE; j++) {
    //         for (int i = 0; i<m.nP; i++) {
    //             par_shape(i*2+1,j*NPAR_PER_NODE+0) = Attr(i,j);
    //         }
    //     }
    // }

    // {
    //     std::vector<std::tuple<std::string, int, std::span<double> > >fields;
    //     for (int j = 0; j<NPAR_SHAPE; j++) {
    //         char buf[1024];
    //         sprintf(buf, "par_%02d", j);
    //         fields.push_back(std::make_tuple(
    //             std::string(buf),
    //             2,
    //             std::span(par_shape.col(j).data(),par_shape.col(j).size())
    //         ));
    //     }
    //     write_vtu("output/par.vtu", to_span(m.P), to_span(m.T), fields);
    // }

    size_t NPAR_DEPTH=0;
    SpMat par_depth;
    {
        size_t npar=0;
        std::vector<Trip> coef;
        for (size_t i=0;i<m.nP;i++) {
            double x = m.P(0,i);
            if ((3<x) && (x<7)) {
                coef.push_back(Trip(i,npar,1));
                npar++;
            }
        }
        NPAR_DEPTH = npar;
        printf("depth parameters: %ld (coef.size: %ld)\n", NPAR_DEPTH, coef.size());
        par_depth.resize(m.nP,NPAR_DEPTH);
        par_depth.setFromTriplets(coef.begin(), coef.end());
    }
    size_t NPAR = NPAR_SHAPE + NPAR_DEPTH;

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
        for (size_t i=0;i<m.nT;i++) {
            size_t Ti[3];
            for (char j=0;j<3;j++) Ti[j] = m.T(j,i);
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
    triangles = m.T.data();
    n_triangles = m.nT;
    n_points = m.nP;
    boundary = m.B.data();
    boundary_flag = m.B_flag.data();
    n_boundary = m.nB;

    std::vector<bool> P_bord(m.nP);
    std::map<size_t, std::vector< Eigen::Vector2d > > bvec;
    for (size_t i=0;i<P_bord.size();i++) P_bord[i] = false;
    for (size_t i=0;i<m.nB;i++) {
        size_t i0 = m.B(0,i);
        size_t i1 = m.B(1,i);
        Eigen::Vector2d p0 = m.P.col(i0);
        Eigen::Vector2d p1 = m.P.col(i1);
        Eigen::Vector2d v = p1 - p0;
        //v.normalize();
        double tmp = v(1); v(1) = -v(0); v(0) = tmp;
        bvec[i0].push_back(v);
        bvec[i1].push_back(v);
        P_bord[i0] = true;
        P_bord[i1] = true;
    }

    size_t nBP = 0;
    for (size_t i=0;i<m.nP;i++) if (P_bord[i]) nBP++;
    Eigen::Array<size_t, Eigen::Dynamic, 1> border_indexes(nBP,1);
    Eigen::Matrix<double, 4, Eigen::Dynamic> border_directions(4,nBP);
    Eigen::Matrix<double, 2, Eigen::Dynamic> border_coef(2,nBP);

    {
        size_t j = 0;
        for (size_t i=0;i<m.nP;i++) if (P_bord[i]) {
            border_indexes(j) = i;
            std::vector< Eigen::Vector2d > vecs = bvec[i];
            assert(vecs.size() == 2);
            Eigen::Vector2d v0 = vecs[0];
            Eigen::Vector2d v1 = vecs[1];
            Eigen::Vector2d w = v0 + v1;
            w.normalize();
            Eigen::Vector2d v = w;
            double tmp = v(1); v(1) = -v(0); v(0) = tmp;
            v = v0; w = v1;
            border_directions(0,j) = w(0);
            border_directions(1,j) = w(1);
            border_directions(2,j) = v(0);
            border_directions(3,j) = v(1);
            border_coef(0,j) = 1;
            border_coef(1,j) = 1;
            j++;
        }
        assert(j == nBP);
    }

    {
        std::vector<std::tuple<std::string, int, std::span<double> > >fields;
        Eigen::Matrix<double, 3, Eigen::Dynamic> v1(3,m.nP); v1.setZero();
        Eigen::Matrix<double, 3, Eigen::Dynamic> v2(3,m.nP); v2.setZero();
        for (size_t i=0;i<nBP;i++)  {
            size_t j = border_indexes[i];
            v1(0,j) = border_directions(0,i);
            v1(1,j) = border_directions(1,i);
            v2(0,j) = border_directions(2,i);
            v2(1,j) = border_directions(3,i);
        }
        fields.push_back(std::make_tuple(
            "v1",
            3,
            std::span(v1.data(),v1.size())
        ));
        fields.push_back(std::make_tuple(
            "v2",
            3,
            std::span(v2.data(),v2.size())
        ));
        write_vtu("output/bord.vtu", to_span(m.P), to_span(m.T), fields);
    }

    return 0;

    SpMat K(DOF,DOF);
    {
        printf("Gathering morphing energy Hessian at 0");
        std::vector<Trip> coef;
        printf(" [mult]");
        Eigen::VectorXd P1b_tmp(DOF);
        Eigen::VectorXd energy_tmp(2);
        Eigen::VectorXd energy_weights(2);
        double YoungMod = 1.0;
        double Poiss = 0.4;
        // energy_weights(0) = YoungMod*Poiss/((1+Poiss)*(1-2*Poiss)); //(tr(E))^2
        // energy_weights(1) = 2*YoungMod/(2*(1+Poiss)); // tr(E^2)
        energy_weights(0) = Poiss; //(tr(E))^2
        energy_weights(1) = (1-2*Poiss); // tr(E^2)
        Eigen::VectorXd Mx(DOF);
        for (size_t k=0; k<maxk; k++) {
            Mx.setZero();
            morph_energy_fix_d(m.P.data(), m.P.data(), ref_x.col(k).data(), Pfix.data(), P1b_tmp.data(), Mx.data(), energy_weights.data());
            for (size_t j=0; j<DOF; j++){
                if (fabs(Mx[j]) > 1e-6) {
                    coef.push_back(Trip(j,ref_j(j,k),Mx[j]));
                }
            }
        }
        printf(" [sparse]");
        K.setFromTriplets(coef.begin(), coef.end());
        printf(" [done]\n");
    }

    {
        for (size_t i=0; i<DOF; i++) {
            for (int j = 0; j<NPAR_SHAPE; j++) {
                if (!P_bord[i]){
                    par_shape(i,j) = 0;
                }                
            }
        }
        printf("Solving linear problem");
        Eigen::KLU<SpMat> solver;  // performs a Cholesky factorization of A
        printf(" [compute]");
        solver.compute(K); assert(solver.info() == Eigen::Success);
        printf(" [solve]");
        for (int j = 0; j<NPAR_SHAPE; j++) {
            par_shape.col(j) = solver.solve(par_shape.col(j));
            assert(solver.info() == Eigen::Success);
        }
        printf(" [done]\n");
    }
    

    {
        std::vector<std::tuple<std::string, int, std::span<double> > >fields;
        for (int j = 0; j<NPAR_SHAPE; j++) {
            char buf[1024];
            sprintf(buf, "par_%02d", j);
            fields.push_back(std::make_tuple(
                std::string(buf),
                2,
                std::span(par_shape.col(j).data(),par_shape.col(j).size())
            ));
        }
        write_vtu("output/par_morph.vtu", to_span(m.P), to_span(m.T), fields);
    }


//    return 0;

    // problem coefficient
    //double wave_k = 4.0;
    const int NOBJ = 3;
    const size_t KINT = 100;
    Eigen::VectorXd integral_k(KINT);
    Eigen::Matrix<double, NOBJ, Eigen::Dynamic> integral_weights(NOBJ,KINT);
    integral_weights.setZero();
    const double wave_k_min = 0.01 * pi;
    const double wave_k_max = 1 * pi;
    {
        double k_dist = (wave_k_max-wave_k_min)/(KINT - 1);
        for (size_t i=0; i<KINT; i++) {
            double weight = k_dist;
            if ((i == 0) || (i == KINT-1)) weight = weight/2;
            double wave_k = wave_k_min + i*k_dist;
            integral_k(i) = wave_k;
            integral_weights(0,i) = weight;
        }
    }
    Eigen::Matrix<double, 2, Eigen::Dynamic> P0 = m.P;
    Eigen::VectorXd D0(m.nP);
    for (size_t i=0; i<D0.size(); i++) D0(i) = 1;

    int iter = 0;
    const auto& objective = [&](const double *pr_, double* grad_, bool export_all=false) -> double {
        Eigen::Map< const Eigen::VectorXd > pr_shape(pr_, NPAR_SHAPE);
        Eigen::Map< const Eigen::VectorXd > pr_depth(pr_ + NPAR_SHAPE, NPAR_DEPTH);
        m.P = P0 + (par_shape * pr_shape).reshaped(2,m.nP);
        //P = P0;
        Eigen::VectorXd D = D0 + par_depth * pr_depth;
        double total_obj = 0;
        Eigen::Map< Eigen::VectorXd >total_grad_shape(grad_, NPAR_SHAPE);
        Eigen::Map< Eigen::VectorXd >total_grad_depth(grad_ + NPAR_SHAPE, NPAR_DEPTH);
        if (grad_ != NULL) {
            total_grad_shape.setZero();
            total_grad_depth.setZero();
        }
        std::vector<std::pair<double, double> > objs; 
        for(size_t kidx = 0; kidx < KINT; kidx++) {
            double wave_k = integral_k(kidx);
            Eigen::VectorXd weights = integral_weights.col(kidx);
            Eigen::VectorXd x(DOF); x.setZero();
            Eigen::VectorXd res(DOF);
            Eigen::VectorXd obj(NOBJ);

            problem(wave_k, m.P.data(), D.data(), x.data(), res.data(), obj.data());
            //printf("obj:%lg\n", obj[0]);

            {
                double resL2 = res.norm();
                printf("Residual (before): %lg\n", resL2);
            }
            
            Eigen::VectorXd Pd(DOF);
            Eigen::VectorXd res_tmp(DOF);
            Eigen::VectorXd obj_tmp(NOBJ);

            Eigen::VectorXd Mx(DOF);

            SpMat A(DOF,DOF);
            {
                printf("Gathering Jacobian");
                std::vector<Trip> coef;
                printf(" [mult]");
                for (size_t k=0; k<maxk; k++) {
                    problem_d(wave_k, m.P.data(), D.data(), x.data(), ref_x.col(k).data(), res_tmp.data(), Mx.data(), obj_tmp.data());
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

            

            problem(wave_k, m.P.data(), D.data(), x.data(), res.data(), obj.data());
            {
                double resL2 = res.norm();
                printf("Residual (after): %lg\n", resL2);
            }
            printf("obj:%lg\n", obj[0]);
            total_obj += weights.dot(obj);
            objs.push_back(std::make_pair(wave_k, obj[0]));

            // ADJOINT
            if (grad_ != NULL) {
                Eigen::VectorXd Pb(DOF); Pb.setZero();
                Eigen::VectorXd Db(m.nP); Db.setZero();
                Eigen::VectorXd xb(DOF); xb.setZero();
                Eigen::VectorXd objb(NOBJ); objb.setZero();
                objb = weights;
                problem_bX(wave_k, m.P.data(), D.data(), x.data(), xb.data(), res_tmp.data(), obj_tmp.data(), objb.data());

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
                
                problem_bP(wave_k, m.P.data(), Pb.data(), D.data(), Db.data(), x.data(), res_tmp.data(), resb.data(), obj_tmp.data(), objb.data());

                Eigen::VectorXd grad_shape = Pb.transpose() * par_shape;
                total_grad_shape += grad_shape;
                Eigen::VectorXd grad_depth = Db.transpose() * par_depth;
                total_grad_depth += grad_depth;
            }
            if (export_all || (kidx == KINT-1)) {
                char buf[1024];
                sprintf(buf, "output/res_%lg_%04d.vtu", wave_k, iter);
                write_vtu(buf, to_span(m.P), to_span(m.T), {
                    std::make_tuple(std::string("Eta"), 2, to_span(x)),
                    std::make_tuple(std::string("depth"), 1, to_span(D))
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
        // {
        //     char buf[1024];
        //     sprintf(buf, "output/res_%04d.points", iter);
        //     FILE* f = fopen(buf, "w");
        //     for (size_t i=0;i<nP;i++) {
        //         fprintf(f, "%.15lg %.15lg\n", m.P(0,i), m.P(1,i));
        //     }
        //     fclose(f);
        // }
        // {
        //     char buf[1024];
        //     sprintf(buf, "output/mesh_%04d.triangles", iter);
        //     FILE* f = fopen(buf, "w");
        //     for (size_t i=0;i<m.nT;i++) {
        //         fprintf(f, "%ld %ld %ld\n", T[3*i+0], T[3*i+1], T[3*i+2]);
        //     }
        //     fclose(f);
        // }
        iter++;
        return total_obj;
    };
    
    if (false) { // FD test
        Eigen::VectorXd pr(NPAR); pr.setZero();
        Eigen::VectorXd gr(NPAR);
        pr(2) += 0.05;
        pr(5) += -0.05;
        double val = objective(pr.data(), gr.data(), true);
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
    Eigen::VectorXd lower(NPAR); lower.setZero();
    Eigen::VectorXd upper(NPAR); upper.setZero();
    for (size_t i=0;i<NPAR_SHAPE; i++) {
        lower(i) = -0.2;
        upper(i) =  1.0;
    }
    for (size_t i=0;i<NPAR_DEPTH; i++) {
        lower(i + NPAR_SHAPE) = -0.5;
        upper(i + NPAR_SHAPE) =  10;
    }

    opt_res = nlopt_set_lower_bounds(opt, lower.data());
    opt_res = nlopt_set_upper_bounds(opt, upper.data());
    opt_res = nlopt_set_maxeval(opt, 500);

    Eigen::VectorXd pr(NPAR); pr.setZero();
    double obj = 0;
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
