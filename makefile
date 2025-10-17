CXXFLAGS+=-I/usr/include/eigen3/ -I/usr/include/suitesparse/ -std=c++20
CXXFLAGS+=-O3 -fopenmp
CFLAGS+=-O3 -fopenmp
# CXXFLAGS+=-g
# CFLAGS+=-g
TAPENADE=$(HOME)/tapenade/tapenade_3.16/bin/tapenade

main: main.o vtu_write.o solve.o problem.o problem_d.o problem_b.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lklu -lnlopt

problem_d.c : problem.c
	$(TAPENADE) -fixinterface -d -head 'problem(res)/(x)' $<

problem_b.c : problem.c
	$(TAPENADE) -fixinterface -b -head 'problem[X](obj)/(x)' -head 'problem[P](res,obj)/(points,depth)' $<
	sed -e '/adStack/s|^|//|' -i $@
	./ADmod.R -f $@ -o $@

morph_energy_b.c : morph_energy.c
	$(TAPENADE) -fixinterface -b -head 'morph_energy(energy)/(P1)' $<
	sed -e '/adStack/s|^|//|' -i $@
	./ADmod.R -f $@ -o $@

morph_energy_b_d.c : morph_energy_b.c
	$(TAPENADE) -fixinterface -d -head 'morph_energy_b(P1b)/(P1)' $<
