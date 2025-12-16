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
#include <thread>
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
    void morph_energy_fix(const double *P0, const double *P1, const double *dir_disp, double *res, double *energyb);
    void morph_energy_fix_d(const double *P0, const double *P1, const double *P1d, const double *dir_disp, double *res, double *resd, double *energyb);
    void morph_energy_fix_b(const double *P0, const double *P1, const double *dir_disp, double *dir_dispb, double *res, double *resb, double *energyb);
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
        f = fopen_safe(mesh_name+".node","r");
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
        f = fopen_safe(mesh_name+".ele","r");
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
        f = fopen_safe(mesh_name+".poly","r");
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
    int write_mesh(const std::string& mesh_name) {
        FILE* f;
        int dim=2, el_size=3;
        f = fopen_safe(mesh_name+".node","w");
        fprintf(f, "%ld %d %ld %d\n", nP, dim, nAttr, 0/*bord*/);
        for (size_t i=0;i<nP;i++) {
            fprintf(f, "%ld", (long int) i+1);
            for (int j=0;j<dim;j++) fprintf(f," %lf", P(j,i));
            for (size_t j=0;j<nAttr;j++) fprintf(f," %lf", Attr(i,j));
            fprintf(f,"\n");
        }
        fclose(f);
        f = fopen_safe(mesh_name+".ele","w");
        fprintf(f, "%ld %d %d\n", nT, el_size, 0/*bord*/);
        for (size_t i=0;i<nT;i++) {
            size_t idx;
            fprintf(f, "%ld", (long int) i+1);
            for (int j=0;j<el_size;j++) fprintf(f," %ld", T(j,i)+1);
            fprintf(f,"\n");
        }
        fclose(f);
        f = fopen_safe(mesh_name+".poly","w");
        int bord = 1;
        fprintf(f, "%ld %d %ld %d\n", 0/*n points*/, dim, nAttr, bord);
        const int edge_size = 2;
        fprintf(f, "%ld %d\n", nB, 1/*bord*/);        
        for (size_t i=0;i<nB;i++) {
            fprintf(f, "%ld", (long int) i+1);
            for (int j=0;j<edge_size;j++) fprintf(f," %ld", B(j,i)+1);
            for (int j=0;j<bord;j++) fprintf(f," %d", B_flag(i,j));
            fprintf(f,"\n");
        }
        fprintf(f,"0\n");
        fprintf(f,"# Generated by mildslope\n");
        fclose(f);
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
    double mesh_area = 0.05;
    double mesh_area_coef = 0.5;
    double mesh_area_limit = 0.0002;
    //mesh_area = mesh_area_limit;
    {
        char buf[1024];
        sprintf(buf,"triangle/triangle -a%.6lf -q30 -j -p %s.poly", mesh_area, mesh_name.c_str());
        system(buf);
    }

    double ftol_abs = 1e-8;
    int mesh_idx = 1;
    int iter = 0;

start:
    {
        char buf[1024];
        sprintf(buf,"%s.%d", mesh_name.c_str(), mesh_idx);
        m.read_mesh(buf);
    }

    const size_t DOFperP = 2;
    const size_t DOF = m.nP*DOFperP;

    size_t NPAR_DEPTH=0;
    SpMat par_depth;
    size_t NPAR_SHAPE = 0;
    SpMat par_shape;

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

    SpMat border_mat;
    SpMat edge_mass_mat;
    size_t nBP = 0;
    Eigen::Array<size_t, Eigen::Dynamic, 1> border_indexes;
    Eigen::Matrix<double, 2, Eigen::Dynamic> border_directions;
    Eigen::Matrix<double, 2, Eigen::Dynamic > border_fix;
    Eigen::Matrix<double, Eigen::Dynamic, 2> par_shape_limits;
    {
        std::set< std::array< size_t, 2 > > bpairs;
        for (size_t i=0;i<m.nB;i++) bpairs.emplace(std::array<size_t, 2>{m.B(0,i), m.B(1,i)});            
        std::map<size_t, std::vector< Eigen::Vector2d > > bvec;
        std::vector<bool> P_bord(m.nP);
        for (size_t i=0;i<P_bord.size();i++) P_bord[i] = false;
        edge_mass_mat.resize(2*m.nP,2*m.nP);
        {
            std::vector<Trip> coef;
            auto fun = [&](size_t i0, size_t i1, size_t i2) {
                if (bpairs.count({i0,i1}) + bpairs.count({i1,i0}) == 0) return;
                P_bord[i0] = true;
                P_bord[i1] = true;
                Eigen::Vector2d p0 = m.P.col(i0);
                Eigen::Vector2d p1 = m.P.col(i1);
                Eigen::Vector2d p2 = m.P.col(i2);
                Eigen::Vector2d v = p1 - p0;
                double len = v.norm();
                double K_in_mat = 0.4;
                for (int k=0; k<2; k++) {
                    coef.push_back(Trip(i0*2+k,i0*2+k,len/3 + K_in_mat/(2*len)));
                    coef.push_back(Trip(i0*2+k,i1*2+k,len/6 - K_in_mat/(2*len)));
                    coef.push_back(Trip(i1*2+k,i0*2+k,len/6 - K_in_mat/(2*len)));
                    coef.push_back(Trip(i1*2+k,i1*2+k,len/3 + K_in_mat/(2*len)));
                }
                //v.normalize();
                double tmp = v(1); v(1) = -v(0); v(0) = tmp;
                Eigen::Vector2d w = p2 - p0;
                if (v.dot(w) > 0.0) v = -v;
                bvec[i0].push_back(v);
                bvec[i1].push_back(v);
                
            };
            for (size_t i=0;i<m.nT;i++) {
                size_t i0 = m.T(0,i);
                size_t i1 = m.T(1,i);
                size_t i2 = m.T(2,i);
                fun(i0,i1,i2);
                fun(i1,i2,i0);
                fun(i2,i0,i1);
            }
            
            edge_mass_mat.setFromTriplets(coef.begin(), coef.end());
        }
        for (size_t i=0;i<m.nP;i++) if (P_bord[i]) nBP++;
        border_indexes.resize(nBP,1);
        border_directions.resize(2,nBP);
        border_fix.resize(2,nBP);
        border_mat.resize(2*m.nP,2*nBP);
        std::vector< std::tuple< size_t, double, double, double > > shape_par;
        {
            std::vector<Trip> coef;
            size_t j = 0;
            for (size_t i=0;i<m.nP;i++) if (P_bord[i]) {
                border_indexes(j) = i;
                std::vector< Eigen::Vector2d > vecs = bvec[i];
                if (vecs.size() != 2) {
                    printf("Number of normal vectors for %ld is %ld\n", i, vecs.size());
                }
                assert(vecs.size() == 2);
                Eigen::Vector2d v0 = vecs[0];
                Eigen::Vector2d v1 = vecs[1];
                Eigen::Vector2d w0 = v0 + v1;
                w0.normalize();
                v0.normalize();
                v1.normalize();
                double skal = v0.dot(v1);
                coef.push_back(Trip(i*2+0,j*2+0, w0(0)));
                coef.push_back(Trip(i*2+1,j*2+0, w0(1)));
                coef.push_back(Trip(i*2+0,j*2+1,-w0(1)));
                coef.push_back(Trip(i*2+1,j*2+1, w0(0)));

                border_directions(0,j) = w0(0);
                border_directions(1,j) = w0(1);
                border_fix(0,j) = 1.0;
                border_fix(1,j) = 1.0;

                double outer, inner, sides, scale=1;
                if (m.nAttr > 0) outer = m.Attr(i,0); else outer = 1;
                if (m.nAttr > 1) inner = m.Attr(i,1); else inner = -outer;
                if (m.nAttr > 2) sides = m.Attr(i,2); else sides = 0.2*(outer-inner);
                assert(outer - inner > -1e-6);
                if (outer - inner > 1e-6) {
                    shape_par.push_back(std::make_tuple(2*j+0,scale,inner,outer));
                    if (fabs(skal) < 0.5) {
                        shape_par.push_back(std::make_tuple(2*j+1,scale,-sides,sides));
                    } else {
                        //border_fix(1,j) = 0.0;
                    }
                }
                j++;
            }
            assert(j == nBP);
            border_mat.setFromTriplets(coef.begin(), coef.end());
        }
        NPAR_SHAPE = shape_par.size();
        par_shape.resize(nBP*2, NPAR_SHAPE);
        par_shape_limits.resize(NPAR_SHAPE,2);
        {
            std::vector<Trip> coef;
            for (size_t i=0; i<NPAR_SHAPE; i++) {
                size_t idx = std::get<0>(shape_par[i]);
                double scale = std::get<1>(shape_par[i]);
                double lower = std::get<2>(shape_par[i]);
                double upper = std::get<3>(shape_par[i]);
                coef.push_back(Trip(idx,i,scale));
                par_shape_limits(i,0) = lower/scale;
                par_shape_limits(i,1) = upper/scale;
            }
            par_shape.setFromTriplets(coef.begin(), coef.end());
        }
    }

    {
        std::vector<std::tuple<std::string, int, std::span<double> > >fields;
        Eigen::Matrix<double, 3, Eigen::Dynamic> v1(3,m.nP); v1.setZero();
        Eigen::Matrix<double, Eigen::Dynamic, 1> fix0(m.nP); fix0.setZero();
        Eigen::Matrix<double, Eigen::Dynamic, 1> fix1(m.nP); fix1.setZero();
        for (size_t i=0;i<nBP;i++)  {
            size_t j = border_indexes[i];
            v1(0,j) = border_directions(0,i);
            v1(1,j) = border_directions(1,i);
            fix0(j) = border_fix(0,i);
            fix1(j) = border_fix(1,i);
        }
        fields.push_back(std::make_tuple(
            "v1",
            3,
            std::span(v1.data(),v1.size())
        ));
        fields.push_back(std::make_tuple("fix0", 1, to_span(fix0)));
        fields.push_back(std::make_tuple("fix1", 1, to_span(fix1)));
        write_vtu("output/bord.vtu", to_span(m.P), to_span(m.T), fields);
    }

    
    SpMat par_shape_mass = par_shape.transpose() * border_mat.transpose() * edge_mass_mat * border_mat * par_shape;

    Eigen::KLU<SpMat> par_shape_mass_inv;
    par_shape_mass_inv.compute(par_shape_mass);
    assert(par_shape_mass_inv.info() == Eigen::Success);
    
    size_t NPAR = NPAR_SHAPE + NPAR_DEPTH;

    n_bord = nBP;
    bord = border_indexes.data();
    fix_dirs = border_directions.data();
    fix = border_fix.data();

    double YoungMod = 1.0;
    double Poiss = 0.4;
    Eigen::VectorXd energy_weights(2);
    // energy_weights(0) = YoungMod*Poiss/((1+Poiss)*(1-2*Poiss)); //(tr(E))^2
    // energy_weights(1) = 2*YoungMod/(2*(1+Poiss)); // tr(E^2)
    energy_weights(0) = Poiss; //(tr(E))^2
    energy_weights(1) = (1-2*Poiss); // tr(E^2)


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
            // if (i*2 < KINT) {
            //     //integral_weights(1,i) = weight;
            // } else {
            //     integral_weights(1,i) = weight;
            // }
        }
    }

    Eigen::VectorXd D0(m.nP);
    Eigen::Matrix<double, 2, Eigen::Dynamic> P1(2,m.nP);
    for (size_t i=0; i<D0.size(); i++) D0(i) = 1;

    //const int iter_morph_ramp = 100;
    const int iter_morph_max = 15;
    const int iter_problem_max = 20;
    
    const auto& objective = [&](const double *pr_, double* grad_, bool export_all=false) -> double {
        Eigen::Map< const Eigen::VectorXd > pr_shape(pr_, NPAR_SHAPE);
        Eigen::Map< const Eigen::VectorXd > pr_depth(pr_ + NPAR_SHAPE, NPAR_DEPTH);
        Eigen::VectorXd D;
        Eigen::VectorXd dir_disp;

        D = D0 + par_depth * pr_depth;
        dir_disp = par_shape * pr_shape;

        P1 = m.P;
        Eigen::VectorXd res_morph(DOF);
        bool do_exit = false;
        SpMat K(DOF,DOF);
        for(int iter_morph = 0; iter_morph < iter_morph_max; iter_morph++) {
            res_morph.setZero();
            morph_energy_fix(m.P.data(), P1.data(), dir_disp.data(), res_morph.data(), energy_weights.data());
            if (export_all) {
                std::vector<std::tuple<std::string, int, std::span<double> > >fields;
                fields.push_back(std::make_tuple("res", 2, to_span(res_morph)));
                char buf[1024];
                sprintf(buf, "output/morph_%04d_%04d.vtu", iter, iter_morph);
                write_vtu(buf, to_span(P1), to_span(m.T), fields);
            }

            {
                double resL2 = res_morph.norm();
                printf("%8d %3d Residual (morph): %lg\n", iter, iter_morph, resL2);
                if (resL2 < 1e-8) do_exit = true;
            }
            if (grad_ == NULL && do_exit) break;
            
            {
                printf("Gathering morphing energy Hessian at 0");
                std::vector<Trip> coef;
                printf(" [mult]");
                Eigen::VectorXd res_tmp(DOF);
                Eigen::VectorXd energy_tmp(2);
                
                Eigen::VectorXd Mx(DOF);
                for (size_t k=0; k<maxk; k++) {
                    Mx.setZero();
                    //                (    P0    ,    P1    ,           P1d      ,     dir_disp   ,       res     ,   resd   ,  energyb             );
                    morph_energy_fix_d(m.P.data(), P1.data(), ref_x.col(k).data(), dir_disp.data(), res_tmp.data(), Mx.data(), energy_weights.data());
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
            if (grad_ != NULL && do_exit) break;
            {
                printf("Solving linear problem");
                Eigen::KLU<SpMat> solver;  // performs a Cholesky factorization of A
                printf(" [compute]");
                solver.compute(K); assert(solver.info() == Eigen::Success);
                printf(" [solve]");
                Eigen::VectorXd ret = solver.solve(res_morph);
                printf(" [done]\n");
                // double coef = 1;
                // if (iter_morph < iter_morph_ramp) {
                //     coef = (iter_morph+1.0) / iter_morph_ramp;
                // }

                P1 = P1 - ret.reshaped(2, m.nP);// * coef;
            }
        }

        double total_obj = 0;
        Eigen::Map< Eigen::VectorXd >total_grad_shape(grad_, NPAR_SHAPE);
        Eigen::Map< Eigen::VectorXd >total_grad_depth(grad_ + NPAR_SHAPE, NPAR_DEPTH);
        if (grad_ != NULL) {
            total_grad_shape.setZero();
            total_grad_depth.setZero();
        }
        std::mutex outputs_mutex;
        std::mutex print_mutex;
        std::vector<bool> kdone(KINT);
        for (size_t i=0;i<kdone.size();i++) kdone[i] = false;
        Eigen::Array<double, 2+NOBJ, Eigen::Dynamic> objs(2+NOBJ, KINT); objs.setZero();
        Eigen::VectorXd P_grad(DOF); P_grad.setZero();
        Eigen::VectorXd D_grad(m.nP); D_grad.setZero();
        const auto& solve_problem = [&](size_t kidx, bool print=false) {
            double wave_k = integral_k(kidx);
            Eigen::VectorXd weights = integral_weights.col(kidx);
            Eigen::VectorXd x(DOF); x.setZero();
            Eigen::VectorXd res(DOF);
            Eigen::VectorXd obj(NOBJ);
            Eigen::VectorXd res_tmp(DOF);
            Eigen::VectorXd obj_tmp(NOBJ);
            SpMat A(DOF,DOF);
            bool do_exit = false;
            //for(int iter_problem = 0; iter_problem < iter_problem_max; iter_problem++) {
            {   int iter_problem = 0;
                problem(wave_k, P1.data(), D.data(), x.data(), res.data(), obj.data());

                {
                    double resL2 = res.norm();
                    if (print) printf("%8d %6lg %3d Residual (problem): %lg\n", iter, wave_k, iter_problem, resL2);
                    if (resL2 < 1e-8) do_exit = true;
                }
            //    if (grad_ == NULL && do_exit) break;
      
                {
                    Eigen::VectorXd Mx(DOF);
                    if (print) printf("Gathering Jacobian");
                    std::vector<Trip> coef;
                    if (print) printf(" [mult]");
                    for (size_t k=0; k<maxk; k++) {
                        problem_d(wave_k, P1.data(), D.data(), x.data(), ref_x.col(k).data(), res_tmp.data(), Mx.data(), obj_tmp.data());
                        for (size_t j=0; j<DOF; j++){
                            if (fabs(Mx[j]) > 1e-6) {
                                coef.push_back(Trip(j,ref_j(j,k),Mx[j]));
                            }
                        }
                        //printf("mult %ld -> %ld\n", k, coef.size());
                    }
                    if (print) printf(" [sparse]");
                    A.setFromTriplets(coef.begin(), coef.end());
                    if (print) printf(" [done]\n");
                }
            //    if (grad_ != NULL && do_exit) break;
                {
                    if (print) printf("Solving linear problem");
                    Eigen::KLU<SpMat> solver;  // performs a Cholesky factorization of A
                    if (print) printf(" [compute]");
                    solver.compute(A); assert(solver.info() == Eigen::Success);
                    if (print) printf(" [solve]");
                    Eigen::VectorXd ret = solver.solve(res); assert(solver.info() == Eigen::Success);
                    if (print) printf(" [done]\n");
                    for (size_t i=0; i<ret.size(); i++) x[i] -= ret[i];
                }
            }

            problem(wave_k, P1.data(), D.data(), x.data(), res.data(), obj.data());
            {
                double resL2 = res.norm();
                if (print) printf("Residual (after): %lg\n", resL2);
            }
            double wobj = weights.dot(obj);
            if (print) {
                printf("obj:");
                for (int k=0;k<NOBJ;k++) printf(" %lg", obj[k]);
                printf(" -> %lg\n", wobj);
            }
            objs(0,kidx) = wave_k;
            objs(1,kidx) = wobj;
            for (int k=0;k<NOBJ;k++) objs(2+k,kidx) = obj[k];
            Eigen::VectorXd Pb(DOF); Pb.setZero();
            // ADJOINT
            if (grad_ != NULL) {
                
                Eigen::VectorXd Db(m.nP); Db.setZero();
                Eigen::VectorXd xb(DOF); xb.setZero();
                Eigen::VectorXd objb(NOBJ); objb.setZero();
                objb = weights;
                problem_bX(wave_k, P1.data(), D.data(), x.data(), xb.data(), res_tmp.data(), obj_tmp.data(), objb.data());

                Eigen::VectorXd resb;
                {   
                    if (print) printf("Solving adjoint problem");
                    Eigen::KLU<SpMat> solver;  // performs a Cholesky factorization of A
                    if (print) printf(" [compute]");
                    solver.compute(A.transpose()); assert(solver.info() == Eigen::Success);
                    if (print) printf(" [solve]");
                    resb = solver.solve(-xb); assert(solver.info() == Eigen::Success);
                    if (print) printf(" [done]\n");
                }
                
                problem_bP(wave_k, P1.data(), Pb.data(), D.data(), Db.data(), x.data(), res_tmp.data(), resb.data(), obj_tmp.data(), objb.data());
                {
                    std::lock_guard<std::mutex> lock(outputs_mutex);        
                    P_grad += Pb;
                    D_grad += Db;
                }
            }
            if (!print) {
                std::lock_guard<std::mutex> lock(print_mutex);        
                kdone[kidx] = true;
                printf("[");
                for (size_t i=0;i<KINT;i++) if (kdone[i]) printf("X"); else printf(" ");
                printf("]\r"); fflush(stdout);
            }
            if (export_all || (kidx == KINT-1)) {
                char buf[1024];
                sprintf(buf, "output/res_%lg_%04d.vtu", wave_k, iter);
                write_vtu(buf, to_span(P1), to_span(m.T), {
                    std::make_tuple(std::string("Eta"), 2, to_span(x)),
                    std::make_tuple(std::string("Pb"), 2, to_span(Pb)),
                    std::make_tuple(std::string("depth"), 1, to_span(D))
                });
            }
        };
        //for(size_t kidx = 0; kidx < KINT; kidx++) solve_problem(kidx);
        {
            std::atomic<size_t> akidx = 0;
            std::vector<std::jthread> thr;
            for (int i=0;i<8;i++) {
                thr.push_back(std::jthread([&](){
                    for(size_t kidx = (akidx++); kidx < KINT; kidx = (akidx++)) solve_problem(kidx);
                }));
            }
            printf("\n");
        }
        if (grad_ != NULL) {
            Eigen::VectorXd res_morphb;
            {   
                printf("Solving adjoint morph");
                Eigen::KLU<SpMat> solver;  // performs a Cholesky factorization of A
                printf(" [compute]");
                solver.compute(K.transpose()); assert(solver.info() == Eigen::Success);
                printf(" [solve]");
                res_morphb = solver.solve(-P_grad); assert(solver.info() == Eigen::Success);
                printf(" [done]\n");
            }

            Eigen::VectorXd dir_dispb(nBP*2); dir_dispb.setZero();
            morph_energy_fix_b(m.P.data(), P1.data(), dir_disp.data(), dir_dispb.data(), res_morph.data(), res_morphb.data(), energy_weights.data());

            Eigen::VectorXd grad_shape = dir_dispb.transpose() * par_shape;
            total_grad_shape += grad_shape;
            Eigen::VectorXd grad_depth = D_grad.transpose() * par_depth;
            total_grad_depth += grad_depth;
        }
        total_obj = 0;
        for (size_t j=0;j<objs.cols();j++) total_obj += objs(1,j);
        {
            char buf[1024];
            sprintf(buf, "output/res_%04d.csv", iter);
            FILE* f = fopen(buf, "w");
            fprintf(f, "idx,wave_k,objw");
            for (size_t i=0;i<NOBJ;i++) {
                fprintf(f, ",obj%d", (int) i);
            }
            fprintf(f, "\n");
            for (size_t j=0;j<objs.cols();j++) {
                fprintf(f, "%ld", (long int) j);
                for (size_t i=0;i<objs.rows();i++) {
                    fprintf(f, ",%.15lg", objs(i,j));
                }
                fprintf(f, "\n");
            }
            fclose(f);
        }
        // {
        //     char buf[1024];
        //     sprintf(buf, "output/res_%04d.points", iter);
        //     FILE* f = fopen(buf, "w");
        //     for (size_t i=0;i<nP;i++) {
        //         fprintf(f, "%.15lg %.15lg\n", P1(0,i), P1(1,i));
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
        pr(2) += 0.1;
        pr(5) += -0.1;
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

    const auto& precond = [&](const double *grad_, double* direction_) {
        //printf("Called precond!\n"); exit(0);
        Eigen::Map< const Eigen::VectorXd > grad_shape(grad_, NPAR_SHAPE);
        Eigen::Map< const Eigen::VectorXd > grad_depth(grad_ + NPAR_SHAPE, NPAR_DEPTH);
        Eigen::Map< Eigen::VectorXd > direction_shape(direction_, NPAR_SHAPE);
        Eigen::Map< Eigen::VectorXd > direction_depth(direction_ + NPAR_SHAPE, NPAR_DEPTH);
        Eigen::VectorXd ret = par_shape_mass_inv.solve(grad_shape);
        ret = ret * (grad_shape.norm()/ret.norm());
        direction_shape = ret;
        direction_depth = grad_depth;
    };

    using obj_type = decltype(objective);
    using pre_type = decltype(precond);
    struct obj_data_t {
        obj_type& obj;
        pre_type& pre;
    } obj_data{objective, precond};

    nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, NPAR);
    nlopt_result opt_res;
    opt_res = nlopt_set_precond_min_objective(
        opt,
        [](unsigned n, const double* x, double* grad, void* f_data) -> double { return static_cast<obj_data_t*>(f_data)->obj(x, grad); },
        [](unsigned n, const double *x, const double *v, double *vpre, void *f_data) { return static_cast<obj_data_t*>(f_data)->pre(v, vpre); },
        (void*) &obj_data
    );

    Eigen::VectorXd lower(NPAR); lower.setZero();
    Eigen::VectorXd upper(NPAR); upper.setZero();

    for (size_t i=0;i<NPAR_SHAPE; i++) {
        lower(i) =  par_shape_limits(i,0);
        upper(i) =  par_shape_limits(i,1);
    }
    for (size_t i=0;i<NPAR_DEPTH; i++) {
        lower(i + NPAR_SHAPE) = -0.5;
        upper(i + NPAR_SHAPE) =  10;
    }

    opt_res = nlopt_set_lower_bounds(opt, lower.data());
    opt_res = nlopt_set_upper_bounds(opt, upper.data());
    opt_res = nlopt_set_maxeval(opt, 200);
    opt_res = nlopt_set_ftol_abs(opt, ftol_abs);
    Eigen::VectorXd pr(NPAR); pr.setZero();
    double obj = 0;
    opt_res = nlopt_optimize(opt, pr.data(), &obj);
    std::cout << pr << "\n";
    printf("Objective: %lg\n", obj);
    objective(pr.data(),NULL,false);
    m.P = P1;
    mesh_idx++;
    {
        char buf[1024];
        sprintf(buf,"%s.%d", mesh_name.c_str(), mesh_idx);
        m.write_mesh(buf);
    }
    if (mesh_area < mesh_area_limit) return 0;
    mesh_area *= mesh_area_coef;
    ftol_abs *= 0.1;
    if (ftol_abs < 1e-8) ftol_abs = 1e-8;
    {
        char buf[1024];
        sprintf(buf,"triangle/triangle -a%.6lf -q30 -j -p %s.%d", mesh_area, mesh_name.c_str(), mesh_idx);
        system(buf);
    }
    mesh_idx++;
    goto start;
    // for (int k=0; k<60; k++) {
    //     double val = objective(pr.data(), gr.data());
    //     pr += -1e-3*gr;
    // }
    return 0;
}
