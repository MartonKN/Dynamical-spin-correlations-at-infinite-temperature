PNAME=holeMotion
CC=g++
CFLAGS=-I.
LIBS=
DEPS=MC_main.hpp MC_random.hpp MC_coord.hpp MC_Cycle.hpp MC_PermutationCycles.hpp MC_Permutation.hpp MC_holePropagation.hpp MC_holePropagationQState.hpp MC_Monte_Carlo_simulations.hpp MC_holePropagationSpinStateAvg.hpp MC_pathData.hpp MC_singleHoleQState.hpp MC_singleHoleQStateBin.hpp MC_singleHoleQStateContainer.hpp MC_pathContainer.hpp MC_clusterRunQ.hpp MC_holePropagationInfiniteSpin.hpp MC_holePropagationSpinStateAvgFwdBwd.hpp MC_holePropagationInfiniteSpinFwdBwd.hpp MC_holePropagationTraceFwdBwd.hpp
OBJ=MC_main.o MC_random.o MC_Cycle.o MC_PermutationCycles.o MC_Permutation.o MC_holePropagation.o MC_holePropagationQState.o MC_Monte_Carlo_simulations.o MC_holePropagationSpinStateAvg.o MC_pathData.o MC_singleHoleQState.o MC_singleHoleQStateBin.o MC_singleHoleQStateContainer.o MC_pathContainer.o MC_holePropagationInfiniteSpin.o MC_holePropagationSpinStateAvgFwdBwd.o MC_holePropagationInfiniteSpinFwdBwd.o MC_holePropagationTraceFwdBwd.o

all: $(PNAME)

%.o: %.cpp $(DEPS)
	$(CC) -O3 -std=c++0x -c -o $@ $< $(CFLAGS) $(LIBS)

$(PNAME): $(OBJ)
	$(CC) -O3 -std=c++0x -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm  *.o *~
