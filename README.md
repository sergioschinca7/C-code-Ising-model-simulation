# C-code-Ising-model-simulation
Ising Model C code
This code is an Ising model with Metropolis algorithm. There are five independet simulations to inmprove decorrelation, you can set number of steps, size and temperature.

There is an output "Ising_statisc_L.txt" with mean magnetization, mean magnetization square and same for energy. There is quantity of acepted steps with name "acept".

To run this code, you clone this repository and compile these into an executable:

gcc -Wall -O3 -o Ising.exe Ising.c -lm

Once you have the executable, it can be run from the command line:

./Ising.e

NOTE: if you change a parameter in Ising.c you must recompile the code.
