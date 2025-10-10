CXXFLAGS+=-I/usr/include/eigen3/ -I/usr/include/suitesparse/ -std=c++20
CXXFLAGS+=-O3
CFLAGS+=-O3
# CXXFLAGS+=-g
# CFLAGS+=-g
TAPENADE=$(HOME)/tapenade/tapenade_3.16/bin/tapenade

main: main.o vtu_write.o solve.o problem.o problem_d.o problem_b.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lklu -lnlopt

problem_d.c : problem.c
	$(TAPENADE) -fixinterface -d -head 'problem[X](res)/(x)' -head 'problem[P](res)/(points)' $<

problem_b.c : problem.c
	$(TAPENADE) -fixinterface -b -head 'problem(obj)/(x,points)' $<
	sed -e '/adStack/s|^|//|' -i $@
