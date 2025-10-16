#ifndef VTU_WRITE_H
#define VTU_WRITE_H

#include "global.h"

#include <string>
#include <tuple>
#include <span>

void write_vtu(
    char* filename,
    std::span<double> points,
    std::span<size_t> triangles,
    const std::vector<std::tuple<std::string, int, std::span<double>>> point_fields={},
    const std::vector<std::tuple<std::string, int, std::span<double>>> cell_fields={}
);

#endif
