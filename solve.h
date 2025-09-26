#ifndef SOLVE_H
#define SOLVE_H

#include "global.h"
#include <functional>

double skal(const vec& a, const vec& b);
int Solve(std::function<void(const vec&, vec&)> mult, const vec& b, vec& x, int max_iter=2000, double eps=1e-6);

#endif