CXXFLAGS=-O3

bodyforce: main_bodyforce.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o
	g++ -o run_bodyforce main_bodyforce.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o

analysis: main_analysis.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk_mean.o
	g++ -o run_analysis main_analysis.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk_mean.o


#pressure:
	#echo 'ERROR : pressure driven flow unavailable at the moment'
pressure: main_pressure.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o domain_InletOutlet.o square.o force.o write_vtk.o
	g++ -o run_pressure main_pressure.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o domain_InletOutlet.o square.o force.o write_vtk.o

main_bodyforce.o: main_bodyforce.cpp
	g++ -o main_bodyforce.o -c main_bodyforce.cpp $(CXXFLAGS)
main_analysis.o: main_bodyforce.cpp
	g++ -o main_analysis.o -c main_analysis.cpp $(CXXFLAGS)
main_pressure.o: main_pressure.cpp
	g++ -o main_pressure.o -c main_pressure.cpp $(CXXFLAGS)
initialize_lattice_arrays.o: initialize_lattice_arrays.cpp
	g++ -o initialize_lattice_arrays.o -c initialize_lattice_arrays.cpp $(CXXFLAGS)
streamCollCompute.o: streamCollCompute.cpp
	g++ -o streamCollCompute.o -c streamCollCompute.cpp $(CXXFLAGS)
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
	rm -rf bodyforce; rm -rf pressure; rm -rf analysis;

