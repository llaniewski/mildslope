CXXFLAGS+=-I/usr/include/eigen3/ -I/usr/include/suitesparse/ -std=c++20
CXXFLAGS+=-O3 -fopenmp
CFLAGS+=-O3 -fopenmp
CXXFLAGS+=-Wno-unused-result -Wno-write-strings
# CXXFLAGS+=-g
# CFLAGS+=-g
TAPENADE=$(HOME)/tapenade/tapenade_3.16/bin/tapenade

main: main.o vtu_write.o solve.o problem.o problem_d.o problem_b.o morph_energy_fix_g_d.o morph_energy_fix_g_b.o
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

morph_energy_fix_g.c : morph_energy_b.c morph_energy_fix.c
	cat $^ >$@

morph_energy_fix_g_d.c : morph_energy_fix_g.c 
	$(TAPENADE) -fixinterface -d -head 'morph_energy_fix(res)/(P1)' -head 'morph_energy_fix[F](res)/(Pfix)' $<

morph_energy_fix_g_b.c : morph_energy_fix_g.c
	$(TAPENADE) -fixinterface -b -head 'morph_energy_fix(res)/(Pfix)' -copyname "_nodiff2" $<
	sed -e '/adStack/s|^|//|' -i $@
	./ADmod.R -f $@ -o $@