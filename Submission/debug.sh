#!/bin/sh
echo "generating input file:\n"
python3 generate_input.py 10
echo "generated input file already ..."
echo "\ndebugging mode:"
echo "\ncleaning...\n"
make clean
echo "\ncompiling...\n"
make all
echo "\nsequential testing...\n"
./seq_tests
echo "\nmpi testings, np 4...\n"
mpirun -np 4 ./mpi_tests
echo "\ntesting mode, running testing cases... (not finished yet)\n"