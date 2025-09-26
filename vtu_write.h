#ifndef VTU_WRITE_H
#define VTU_WRITE_H

#include "global.h"

#include <string>
#include <tuple>

void write_vtu(char* filename, std::vector<double>& points, std::vector<size_t>& triangles, const std::vector<std::tuple<std::string, int, std::vector<double>*>> fields);

#endif
