CXXFLAGS=-O3

main: main.o vtu_write.o solve.o
	$(CXX) $(CXXFLAGS) -o $@ $^


