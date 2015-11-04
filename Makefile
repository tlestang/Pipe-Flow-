CXXFLAGS=-O3

bodyforce: main_bodyforce.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o
	g++ -o run_bodyforce main_bodyforce.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o

pressure:
	echo 'ERROR : pressure driven flow unavailable at the moment'
# pressure: main_pressure.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o
#	g++ -o run_pressure main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o

main_bodyforce.o: main_bodyforce.cpp
	g++ -o main_bodyforce.o -c main_bodyforce.cpp $(CXXFLAGS)
main_pressure.o: main_pressure.cpp
	g++ -o main_pressure.o -c main_pressure.cpp $(CXXFLAGS)
initialize_lattice_arrays.o: initialize_lattice_arrays.cpp
	g++ -o initialize_lattice_arrays.o -c initialize_lattice_arrays.cpp $(CXXFLAGS)
streamCollCompute.o: streamCollCompute.cpp
	g++ -o streamCollCompute.o -c streamCollCompute.cpp $(CXXFLAGS)
domain_noSlipWalls.o: domain_noSlipWalls.cpp
	g++ -o domain_noSlipWalls.o -c domain_noSlipWalls.cpp $(CXXFLAGS)
square.o: square.cpp
	g++ -o square.o -c square.cpp $(CXXFLAGS)
force.o: force.cpp
	g++ -o force.o -c force.cpp $(CXXFLAGS)
write_vtk.o: write_vtk.cpp
	g++ -o write_vtk.o -c write_vtk.cpp $(CXXFLAGS)
clean:
	rm -rf *.o
mrproper: clean
	rm -rf bodyforce
