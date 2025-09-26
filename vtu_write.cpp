#include "vtu_write.h"

#include <stdio.h>
#include <stdlib.h>

void write_vtu(char* filename, vec& points, std::vector<size_t>& triangles, std::vector<std::tuple<std::string, int, vec*>> fields) {
    FILE* f = fopen(filename, "w");
    if (f == NULL) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(-1);
    }
    size_t npoints = points.size()/2;
    size_t ncells = triangles.size()/3;
    if ((points.size() % 2 != 0) || (triangles.size() % 3 != 0)) {
        fprintf(stderr, "Wrong number of points or triangles\n");
        exit(-1);
    }
    fprintf(f, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    fprintf(f, "  <UnstructuredGrid>\n");
    fprintf(f, "    <Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\">\n", npoints, ncells);
    fprintf(f, "      <PointData>\n");
    for (const auto& field : fields) {
        const std::string name = std::get<0>(field);
        const int comp = std::get<1>(field);
        const vec& vals = *std::get<2>(field);
        printf("Writing %s with %d components (%ld values)\n", name.c_str(), comp, vals.size());
        if (vals.size() != npoints*comp) {
            fprintf(stderr, "Wrong number of elements in field %s\n", name.c_str());
            exit(-1);
        }
        fprintf(f, "        <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\">\n", name.c_str(), comp);
        for (size_t i=0; i<vals.size(); i++) {
            fprintf(f, " %.15lg", vals[i]);
            if (i % 6 == 0) fprintf(f, "\n");
        }
        fprintf(f, "\n");
        fprintf(f, "        </DataArray>\n");
    }
    fprintf(f, "      </PointData>\n");
    fprintf(f, "      <CellData>\n");
    fprintf(f, "      </CellData>\n");
    fprintf(f, "      <Points>\n");
    fprintf(f, "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">");
    for (size_t i=0; i<points.size(); i++) {
        if (i % 6 == 0) fprintf(f, "\n");
        fprintf(f, " %.15lg", points[i]);
        if (i % 2 == 1) fprintf(f, " 0");
    }
    fprintf(f, "\n");
    fprintf(f, "        </DataArray>\n");
    fprintf(f, "      </Points>\n");
    fprintf(f, "      <Cells>\n");
    fprintf(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for (size_t i=0; i<triangles.size(); i++) {
        fprintf(f, " %ld", triangles[i]);
        if (i % 6 == 0) fprintf(f, "\n");
    }
    fprintf(f, "\n");
    fprintf(f, "        </DataArray>\n");
    fprintf(f, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">");
    for (size_t i=0; i<triangles.size(); i++) {
        if (i % 6 == 0) fprintf(f, "\n");
        fprintf(f, " %ld", (i+1)*3);
    }
    fprintf(f, "\n");
    fprintf(f, "        </DataArray>\n");
    fprintf(f, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">");
    for (size_t i=0; i<triangles.size(); i++) {
        if (i % 6 == 0) fprintf(f, "\n");
        fprintf(f, " %d", 5);
    }
    fprintf(f, "\n");
    fprintf(f, "        </DataArray>\n");
    fprintf(f, "      </Cells>\n");
    fprintf(f, "    </Piece>\n");
    fprintf(f, "  </UnstructuredGrid>\n");
    fprintf(f, "</VTKFile>\n");
    fprintf(f, "\n");
    fclose(f);
}
