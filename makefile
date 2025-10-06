CXXFLAGS+=-O3 -I/usr/include/eigen3/ -I/usr/include/suitesparse/ -std=c++20
CFLAGS+=-O3
TAPENADE=$(HOME)/tapenade/tapenade_3.16/bin/tapenade

main: main.o vtu_write.o solve.o problem.o problem_d.o problem_b.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lklu

problem_d.c : problem.c
	$(TAPENADE) -d -head 'problem(res,obj)/(x,points)' $<

problem_b.c : problem.c
	$(TAPENADE) -b -head 'problem(obj)/(x,points)' $<
