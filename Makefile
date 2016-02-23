CXXFLAGS=-O3 

bodyforce: main_bodyforce.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o
	g++ -o run_bodyforce main_bodyforce.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o

init: main_init.o initialize_lattice_arrays.o streamCollComputeInit.o domain_noSlipWalls.o square.o force.o write_vtk.o
	g++ -o run_init main_init.o initialize_lattice_arrays.o streamCollComputeInit.o domain_noSlipWalls.o square.o force.o write_vtk.o 	

prog: main_prog.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o
	g++ -o run_prog main_prog.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o

random: main_random.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o
	g++ -o run_random main_random.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o 

analysis: main_analysis.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk_mean.o
	g++ -o run_analysis main_analysis.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk_mean.o

pressure: main_pressure.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o domain_InletOutlet.o square.o force.o write_vtk.o
	g++ -o run_pressure main_pressure.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o domain_InletOutlet.o square.o force.o write_vtk.o

benchmark: main_benchmark.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o
	g++ -o run_benchmark main_benchmark.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o 

# --- MAIN FUNCTIONS

main_bodyforce.o: main_bodyforce.cpp
	g++ -o main_bodyforce.o -c main_bodyforce.cpp $(CXXFLAGS)
main_init.o: main_init.cpp
	g++ -o main_init.o -c main_init.cpp $(CXXFLAGS)
main_analysis.o: main_bodyforce.cpp
	g++ -o main_analysis.o -c main_analysis.cpp $(CXXFLAGS)
main_pressure.o: main_pressure.cpp
	g++ -o main_pressure.o -c main_pressure.cpp $(CXXFLAGS)
main_prog.o: main_prog.cpp
	g++ -o main_prog.o -c main_prog.cpp $(CXXFLAGS)
main_random.o: main_random.cpp
	g++ -o main_random.o -c main_random.cpp $(CXXFLAGS)
main_benchmark.o: main_benchmark.cpp
	g++ -o main_benchmark.o -c main_benchmark.cpp $(CXXFLAGS)

# --- CORE FUNCTIONS ---

initialize_lattice_arrays.o: initialize_lattice_arrays.cpp
	g++ -o initialize_lattice_arrays.o -c initialize_lattice_arrays.cpp $(CXXFLAGS)
streamCollCompute.o: streamCollCompute.cpp
	g++ -o streamCollCompute.o -c streamCollCompute.cpp $(CXXFLAGS)
streamCollComputeInit.o: streamCollComputeInit.cpp
	g++ -o streamCollComputeInit.o -c streamCollComputeInit.cpp $(CXXFLAGS)
domain_noSlipWalls.o: domain_noSlipWalls.cpp
	g++ -o domain_noSlipWalls.o -c domain_noSlipWalls.cpp $(CXXFLAGS)
domain_InletOutlet.o: domain_InletOutlet.cpp
	g++ -o domain_InletOutlet.o -c domain_InletOutlet.cpp $(CXXFLAGS)
square.o: square.cpp
	g++ -o square.o -c square.cpp $(CXXFLAGS)
force.o: force.cpp
	g++ -o force.o -c force.cpp $(CXXFLAGS)
write_vtk.o: write_vtk.cpp
	g++ -o write_vtk.o -c write_vtk.cpp $(CXXFLAGS)
write_vtk_mean.o: write_vtk_mean.cpp
	g++ -o write_vtk_mean.o -c write_vtk_mean.cpp $(CXXFLAGS)
clean:
	rm -rf *.o
mrproper: clean
	rm -rf bodyforce; rm -rf pressure; rm -rf analysis; rm -rf prog;

