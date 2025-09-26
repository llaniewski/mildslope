#include "solve.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double skal(const vec& a, const vec& b) {
    size_t n = a.size();
    double sum = 0;
    for (int i=0;i<n;i++) sum += a[i]*b[i];
    return sum;
}


int Solve(std::function<void(const vec&, vec&)> mult, const vec& b, vec& x, int max_iter, double eps) {
    size_t n = b.size();
    static vec r; r.resize(n);
    static vec p; p.resize(n);
    static vec Ap; Ap.resize(n);
    static vec q; q.resize(n);
    static vec Aq; Aq.resize(n);
    double res;
    for (int iter = 0; iter < max_iter; iter++) {
        mult(x, r);
        for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
        res = sqrt(skal(r, r));
        //res_draw(res);
        printf("res=%lg (%d)\n", res, iter);
        if (res < eps) {
            printf("linear %5d iterations converged (%lg)\n", iter, res);
            return iter;
        }
        for (int i = 0; i < n; i++) {
            p[i] = r[i];
        }
        if (iter > 0) {
            mult(p, Ap);
            mult(q, Aq);
            double beta = skal(Aq, Ap) / skal(Aq, Aq);
            for (int i = 0; i < n; i++) p[i] = p[i] - q[i] * beta;
        }
        mult(p, Ap);
        //printf("|p|^2=%lg\n",skal(p,p));
        //printf("|Ap|^2=%lg\n",skal(Ap,Ap));
        double alpha = skal(Ap, r) / skal(Ap, Ap);
        //printf("alpha=%lg\n",alpha);
        for (int i = 0; i < n; i++) x[i] = x[i] + p[i]*alpha;
        for (int i = 0; i < n; i++) q[i] = p[i];
    }
    printf("linear unconverged final residual=%lg\n", res);
    return max_iter;
}


int SolveCG(std::function<void(const vec&, vec&)> mult, const vec& b, vec& x, int max_iter, double eps) {
    size_t n = b.size();
    static vec r; r.resize(n);
    static vec p; p.resize(n);
    static vec Ap; Ap.resize(n);
    static vec q; q.resize(n);
    static vec Aq; Aq.resize(n);
    double res;
    for (int iter = 0; iter < max_iter; iter++) {
        mult(x, r);
        for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
        res = sqrt(skal(r, r));
        //res_draw(res);
        printf("res=%lg (%d)\n", res, iter);
        if (res < eps) {
            printf("linear %5d iterations converged (%lg)\n", iter, res);
            return iter;
        }
        for (int i = 0; i < n; i++) {
            p[i] = r[i];
        }
        if (iter > 0) {
            mult(p, Ap);
            mult(q, Aq);
            double beta = skal(Aq, p) / skal(Aq, q);
            for (int i = 0; i < n; i++) p[i] = p[i] - q[i] * beta;
        }
        mult(p, Ap);
        //printf("|p|^2=%lg\n",skal(p,p));
        //printf("|Ap|^2=%lg\n",skal(Ap,Ap));
        double alpha = skal(p, r) / skal(Ap, p);
        //printf("alpha=%lg\n",alpha);
        for (int i = 0; i < n; i++) x[i] = x[i] + p[i]*alpha;
        for (int i = 0; i < n; i++) q[i] = p[i];
    }
    printf("linear unconverged final residual=%lg\n", res);
    return max_iter;
}
