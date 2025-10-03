CXXFLAGS+=-O3 -I/usr/include/eigen3/ -I/usr/include/suitesparse/

main: main.o vtu_write.o solve.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lklu


