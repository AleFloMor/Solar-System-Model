# Solar-System-Model

I have provided all necessary files to run the simulation.
How to run the simulation:

1. Open your terminal (use Anaconda Prompt on Windows).
2. Activate Anaconda's base environment (if necessary): conda activate base
3. Navigate to the directory containing the provided simulation files: cd "path/to/Solar System Code"
4. Run the simulation using the command: python main.py parameters.txt


The following is for if you want to edit your parameter file:

The user must provide a text file with the parameters in the terminal (e.g. parameters.txt) in the following format:

numsteps = #
stepsize = #
input_file = '...'.txt
output_file = '...'.xyz

Ensure that:
- The keys i.e.numsteps, stepsize etc. are spelt as above
- There is an '=' sign between the key and the parameter
The code will pick up if there are any problems with the user's parameter input.
The step size is in days.
The input file  must describe each particle correctly in the format required by Particle3D as given in "solar_system.txt".
This code is for a Solar System simulation, please ensure there is a body for the Sun.
