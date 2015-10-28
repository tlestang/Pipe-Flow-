all: run

run: main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o
	g++ -o run main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o

main.o: main.cpp
	g++ -o main.o -c main.cpp
initialize_lattice_arrays.o: initialize_lattice_arrays.cpp
	g++ -o initialize_lattice_arrays.o -c initialize_lattice_arrays.cpp
streamCollCompute.o: streamCollCompute.cpp
	g++ -o streamCollCompute.o -c streamCollCompute.cpp
domain_noSlipWalls.o: domain_noSlipWalls.cpp
	g++ -o domain_noSlipWalls.o -c domain_noSlipWalls.cpp
square.o: square.cpp
	g++ -o square.o -c square.cpp
force.o: force.cpp
	g++ -o force.o -c force.cpp
write_vtk.o: write_vtk.cpp
	g++ -o write_vtk.o -c write_vtk.cpp
clean:
	rm -rf *.o
mrproper: clean
	rm -rf pois
