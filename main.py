"""

Student Number: s2155825
NOTE: Pycharm was used to write, run and test the code

"""

import sys
import numpy as np
import matplotlib.pyplot as pyplot
from particle3D import Particle3D
import scipy as scipy

def initial_conditions(particles, positions, potential, n, output, particle_indexes, relative_distance, orbit_indexes, energies):

    """

    Initialises the trajectory plotting array, the energy array, and subtracts the system's centre-of-mass velocity.

    """

    # Calculates the initial centre-of-mass velocity of the system
    CoM_velocity = Particle3D.com_velocity(particles)
    output.write(f"{n} \n")
    output.write(f"0 \n")

    for i in range(n):
        # Subtracts the centre-of-mass velocity from each particle in the system
        particles[i].velocity -= CoM_velocity

        # indexes which body each particle is orbiting around
        if particles[i].label.lower() == 'moon':
            orbit_indexes[i] = particle_indexes['earth']
        else:
            orbit_indexes[i] = particle_indexes['sun']

        # Stores the initial position of the particles in the trajectories array
        positions[i, 0] = particles[i].position

        # Writes the initial positions of the particles to file
        output.write(f"{(particles[i].__str__())} \n")

        # Stores the initial distance between each particle and the body they are orbiting
        relative_distance[i, 0] = np.linalg.norm(positions[i, 0] - positions[orbit_indexes[i], 0])

    # Stores the initial energy of the system
    energies[0] = total_energy(particles, potential)

def total_energy(particles, pe):

    """

    Sums the kinetic and potential energy of the system to return the total energy.

    """

    # Calculates the kinetic energy of the system
    ke = Particle3D.total_kinetic_energy(particles)

    # Calculates the total energy of the system
    te = ke + pe

    return te

def parameter_validation(file_name):

    """

    Checks the parameters taken from the text file provided by ensuring that:
    - The parameter names are correct
    - The files and values are of the correct type
    - There are no missing parameters

    Outputs the error type to the user if there are any

    """

    # Sets a dictionary for the data/file type expected for each parameter
    expected_parameters = {'numsteps': int, 'stepsize': float, 'input_file': '.txt', 'output_file': '.xyz'}

    # Creates an empty dictionary for the valid parameters
    valid_parameters = {}

    # Iterates for each line in the file
    for line in open(file_name, 'r'):

        # If the line is empty, the function goes to the next line
        if len(line.split()) == 0:
            continue

        # Checks that there is an equals sign
        if '=' not in line:
            print("Error: There must be an equals sign = separating the key and the value of the parameter")
            return None

        # Splits the line into a list containing the string before and after the equals sign
        line = line.split('=')

        # Strips the blank spaces and stores the key and value of the parameter
        # If there is are capital letters in the key, the key will still match if it is spelt correctly
        key = line[0].strip().lower()
        value = line[1].strip()

        # If the key is not in the dictionary for the expected parameters, prints which key is missing
        if key not in expected_parameters:
            print(f"Error: Unknown Parameter {key}")
            return None

        # If the key is a match, it checks that the value attached is of the correct data/file type
        try:
            if key == 'numsteps' and type(int(value)) == expected_parameters['numsteps']:
                valid_parameters[key] = int(value)
            elif key == 'stepsize' and type(float(value)) == expected_parameters['stepsize']:
                valid_parameters[key] = float(value)
            elif key == 'input_file' and value[-4:] == expected_parameters['input_file']:
                valid_parameters[key] = str(value)
            elif key == 'output_file' and value[-4:] == expected_parameters['output_file']:
                valid_parameters[key] = str(value)

        # If value is of the wrong type, and there is no value error, the error must have come from a file
        # the key and expected file type is printed for the user
            else:
                print(f"Error: Invalid filetype for {key}. '{expected_parameters[key]}' expected ")
                return None

        # If there is a value error, the key and expected value type is printed for the user
        except ValueError:
            if key == 'stepsize':
                print(f"Error: Invalid datatype for {key}, float or integer is expected")
            else:
                print(f"Error: Invalid datatype for {key}, '{expected_parameters[key].__name__}' expected")
            return None

    # Checks for missing parameters
    missing_parameters = set(expected_parameters.keys()) - set(valid_parameters.keys())

    # If there are missing parameters, it prints a list of these parameters
    if missing_parameters:
        print(f"Error: Missing parameters: {', '.join(missing_parameters)}.")
        return None

    return valid_parameters

def generate_system(input_file):

    """

    Creates a system of particles from the text file provided by the user

    """
    # Creates blank list to store the individual particles
    particles = []

    # Loop to go through each line of the system file
    for line in open(input_file):

        # Stores the particle object and appends it to the list of particles
        p = Particle3D.read_line(line.strip())
        particles.append(p)

    return particles

def index(particles, particle_indexes, n):

    """

    Creates a dictionary for the name of each body in the system

    """

    for i in range(n):
        particle_indexes[particles[i].label.lower()] = i

def main():

    # Inputs the parameter file name
    file_name = str(sys.argv[1])

    # Validates the parameters
    valid_parameters = parameter_validation(file_name)

    # Ends the program if any of the parameters are not valid
    if valid_parameters is None:
        sys.exit()

    # Stores the validated file names, and parameters
    output_file = valid_parameters['output_file']
    input_file = valid_parameters['input_file']
    numstep = valid_parameters['numsteps']
    dt = valid_parameters['stepsize']
    numdays = numstep * dt
    print(f"Time step = {dt} days")
    print(f"Simulation duration = {numdays} days")

    # Opens the output (XYZ) file to write the trajectories of the particles throughout the simulation
    output = open(output_file, 'w')

    # Sets up the system used to run the simulation
    particles = generate_system(input_file)

    # Uses a converted value of G such that the simulation numbers are reasonably close to 1
    # Earth's mass, AU, and days, are used to convert G
    G = 8.887692593e-10

    # Initialises the time
    time = 0.0

    # Stores the number of particles in the system
    n = len(particles)

    # Sets up the array used to plot the trajectories of the particles
    # Stores the position of each particle for every timestep
    positions = np.zeros((n, numstep, 3))

    # Stores the current time for every timestep
    times = np.zeros(numstep)
    energies = np.zeros(numstep)

    # Stores the relative distance between each particle and the body they are orbiting
    relative_distance = np.zeros((n, numstep))

    # function that creates a dictionary for the name of each body in the system
    particle_indexes = {}
    index(particles, particle_indexes, n)
    # Stores the index of the body each particle is orbiting
    orbit_indexes = list(range(n))

    # Initialises all variables for the system at time = 0
    # Potential is calculated to find the initial total energy
    # Force is calculated to initiate the verlet algorithm
    forces, potential = Particle3D.compute_forces_potential(particles, Particle3D.compute_separations(particles), G)

    initial_conditions(particles, positions, potential, n, output, particle_indexes, relative_distance, orbit_indexes, energies)

    # Begins the Verlet algorithm with the initial positions, velocities, and forces established.

    for i in range(1, numstep):
        # Updates the time of the system by adding a timestep
        time += dt

        # Stores the time of the system
        times[i] = time

        # Writes the number of particles and the iteration to file
        output.write(f"{n} \n")
        output.write(f"{i} \n")

        for j in range(n):
            # Updates the position of the particles for the new time - the 1st step of the Verlet algorithm
            Particle3D.update_position_verlet(particles[j], dt, forces[j])

            # Writes the particle label and updated position to .xyz file in string format used in Particle3D file
            output.write(f"{(particles[j].__str__())} \n")

            # Stores the updated position in the trajectories array
            positions[j, i] = particles[j].position

        # Calculates the new forces acting on the particles - the 2nd step of the Verlet algorithm
        new_forces, potential = Particle3D.compute_forces_potential(particles, Particle3D.compute_separations(particles), G)

        for j in range(n):
            relative_distance[j, i] = np.linalg.norm(positions[j, i] - positions[orbit_indexes[j], i])
            # Updates the velocities of the particles - the 3rd step of the Verlet algorithm
            Particle3D.update_velocity_verlet(particles[j], dt, forces[j], new_forces[j])

        # Stores the new forces as the initial forces such that
        # the algorithm can begin the next iteration of the algorithm
        forces = new_forces
        # Stores the total energy of the system
        energies[i] = total_energy(particles, potential)

    # Finds the orbiting properties of each body
    # Calculates the energy deviation of the system in the simulation
    energy_deviation = abs((max(energies) - min(energies))/energies[0])
    print(f"Energy deviation = {energy_deviation} \n")

    for i in range(n):
        # Only calculates the orbiting properties of the planets and moon
        # Finds the x displacement between each particle and the body they are orbiting
        relative_position = positions[i, :, 0] - positions[orbit_indexes[i], :, 0]

        # Uses library to find when the peaks in the relative position occured
        peaks_array = (scipy.signal.find_peaks(relative_position))
        peaks = peaks_array[0]

        if particles[i].label.lower() != 'sun':
            print(particles[i].label)
            # If the number of peaks observed is less than 2, the period cannot be calculated with this method
            if len(peaks) < 2:
                print(f"Not enough time has elapsed to find a period \n")
            # Otherwise the orbital period, apsides, and semi-major axes are computed
            else:

                # Calculates the orbital period by summing the difference between each peak
                # This is multiplied by dt to ensure that the timescale of the period is correct
                # This is then divided by 1 less than the number of peaks (i.e. the number of orbits)
                # to give the orbital period of the body
                period = dt * sum(peaks[-j] - peaks[-j - 1] for j in range(1, len(peaks))) / (len(peaks[1:]))
                print(f"Orbital period around the {particles[orbit_indexes[i]].label} = {period} days")

                # Prints the correct language for the apsides
                # Given by the minimum and maximum distance between each body and the body they are orbiting
                if particles[i].label.lower() == 'moon':
                    print(f"Perigee = {min(relative_distance[i])} AU")
                    print(f"Apogee = {max(relative_distance[i])} AU")
                else:
                    print(f"Perihelion = {min(relative_distance[i])} AU")
                    print(f"Aphelion = {max(relative_distance[i])} AU")
                # Computes and prints the semi-major axis of the body
                semi_major_axis = (min(relative_distance[i]) + max(relative_distance[i])) / 2
                print(f"Semi-major axis = {semi_major_axis} AU\n")

    # Plots the Energy of the simulation over time
    pyplot.title('Energy vs time ')
    pyplot.xlabel('Time (Days)')
    pyplot.ylabel('Total energy')
    pyplot.plot(times, energies)
    pyplot.show()

    # Plots the particle trajectories in the x and y directions for bodies orbiting the Sun
    pyplot.title('Solar system trajectories about the Sun')
    pyplot.xlabel('X position (AU)')
    pyplot.ylabel('Y position (AU)')
    pyplot.plot(positions[particle_indexes['sun'], :, 0], positions[particle_indexes['sun'], :, 1], label='Sun')
    for i in range(n):
        particle_label = particles[i].label.lower()
        if particle_label != 'sun' and particle_label != 'moon':
            relative_x_positions = positions[particle_indexes[particle_label], :, 0] - positions[particle_indexes['sun'], :, 0]
            relative_y_positions = positions[particle_indexes[particle_label], :, 1] - positions[particle_indexes['sun'], :, 1]
            pyplot.plot(relative_x_positions, relative_y_positions, label=particles[i].label)
    pyplot.legend()
    pyplot.show()

main()
