"""

Particle3D, a class to describe point particles in 3D space
An instance describes a particle in Euclidean 3D space:
velocity and position are [3] arrays
Number: s2155825

"""
import numpy as np


class Particle3D(object):

    """

     Class to describe point-particles in 3D space
     Attributes

    ----------

    label: name of the particle
     mass: mass of the particle
     position: position of the particle
     velocity: velocity of the particle
     Methods

    -------
     __init__
     __str__
     kinetic_energy: computes the kinetic energy
     momentum: computes the linear momentum
     update_position_1st: updates the position to 1st order
     update_position_2nd: updates the position to 2nd order
     update_velocity: updates the velocity

    Static Methods

    --------------

    read_file: initializes a P3D instance from a file handle
     total_kinetic_energy: computes total K.E. of a list of particles

    com_velocity: computes centre-of-mass velocity of a list of particles

    """

    def __init__(self, label, mass, position, velocity):

        """

        Initialises a particle in 3D space.

        label: str

        name of the particle

        mass: float

        mass of the particle

        position: [3] float array

        position vector

        velocity: [3] float array

        velocity vector

        """

        self.label = str(label)

        self.mass = float(mass)

        self.position = np.array(position)

        self.velocity = np.array(velocity)

    def __str__(self):

        """

        Return an XYZ-format string. The format is

        label x y z

        """

        f = self.label + " " + str(self.position[0]) + " " + str(self.position[1]) + " " + str(self.position[2])

        return f

    def kinetic_energy(self):

        """

        Returns the kinetic energy of a Particle3D instance by following the formula: KE = 0.5 * m * v^2

        """

        ke = 0.5 * self.mass * np.linalg.norm(self.velocity) ** 2

        return ke

    def momentum(self):

        """

        follows formula for linear momentum as: momentum = mass * velocity


        """

        p = self.mass * self.velocity

        return p

    def update_position_1st(self, dt):

        """

        follows first position updating algorithm as: new position = old position + velocity * timestep

        """

        self.position += self.velocity * dt

        return self.position

    def update_position_verlet(self, dt, force):

        """

        follows second position updating algorithm as: new position = old position + velocity * timestep + timestep^2 * force/ 2 * mass

        """

        self.position += self.velocity * dt + (dt ** 2) * (force / (2 * self.mass))

        return self.position

    def update_velocity_verlet(self, dt, force1, force2):
        """

        follows velocity updating verlet algorithm as: new velocity = old velocity + timestep * (force_1 + force_2)/(2 * mass)

        """

        self.velocity += dt * ((force1 + force2) / (2 * self.mass))

        return self.velocity

    @staticmethod
    def read_line(line):

        """

        Splits the line of text into a list and then associates the relevant data to the correct label. finally it runs the particle

        through the class to initialise the object

        """

        elements = line.split()

        label = str(elements[0])

        mass = float(elements[1])

        position = np.array([float(elements[2]), float(elements[3]), float(elements[4])])

        velocity = np.array([float(elements[5]), float(elements[6]), float(elements[7])])

        p = Particle3D(label, mass, position, velocity)

        return p

    @staticmethod
    def total_kinetic_energy(particles):

        """

        Creates a loop over the particles calling the kinetic energy function and summing them in this function over the loop

        """

        tke = 0

        for i in range(len(particles)):

            tke += Particle3D.kinetic_energy(particles[i])



        return tke

    @staticmethod
    def com_velocity(particles):

        """

        Computes the CoM velocity of a list of P3D's by summing the momentum of all the particles and dividing the result by the sum of the
        masses

        """

        total_momentum = np.zeros(3)

        total_mass = 0

        for i in range(len(particles)):

            total_momentum += Particle3D.momentum(particles[i])

            total_mass += particles[i].mass

        Com_V = total_momentum / total_mass

        return Com_V

    @staticmethod
    def compute_separations(particles):

        """

        Computes the separation between each particle by iterating through every combination of particles.

        """

        # Finds the number of particles in the system in order to iterate through the particles the correct number of times
        n = len(particles)

        # Creates a 3D numpy array filled with zeros to store the vector separating each particle
        separations = np.zeros((n, n, 3))

        # Iterates through every possible combination of particles
        x = 0
        for j in range(n - 1):
            x += 1
            for i in range(x, n):
                # Calculates the vector separating particle i and particle j
                separation = particles[j].position - particles[i].position

                # Stores the separation into the appropriate index in the array
                separations[i, j] = separation

                # Stores the negative of the separation into the other index associated with the separation of the two particles
                separations[j, i] = -separation

        return separations

    @staticmethod
    def compute_forces_potential(particles, separations, G):

        """

        Computes the gravitational force acting on each particle by iterating through every combination of particles and calculating the resulting force

        """

        # finds the number of particles in the system in order to iterate through the particles the correct number of times
        n = len(particles)

        # creates a 2D numpy array filled with zeros to store the force acting on each particle
        forces = np.zeros((n, 3))

        # creates a 3D numpy array filled with zeros to store the force acting between each particle
        all_forces = np.zeros((n, n, 3))

        # creates a variable to store the sum of the potential energy of the system
        potential = 0

        # iterates through every possible combination of particles
        x = 0
        for j in range(n):
            x += 1
            for i in range(x, n):

                # Computes the distance between the two particles by utilising the separations variable created in the compute_separations function
                distance = np.linalg.norm(separations[i, j])

                # Computes the force acting on particle i due to particle j using the formula for gravitational force:
                # F_ij = (G * m_i * m_j * r_ij)/ r^3
                force = (G * particles[i].mass * particles[j].mass * separations[i, j]) / (distance ** 3)

                # Stores the force acting on particle i due to particle j into the appropriate index in the array
                all_forces[i, j] = force

                # Utilising Newton's third law, we can store the negative of the calculated force into the other index,
                # which is associated with the force between the two particles
                all_forces[j, i] = - force

                # Computes the gravitational potential energy between the two particles using the formula for gravitational potential energy:
                # U_ij = -G * M_i * M_j / r
                # The condition for the computation of the forces means the potential energy is only calculated once
                potential -= float((G * particles[i].mass * particles[j].mass) / distance)

            # Once all the forces acting on a particle has been calculated, the total force acting on the particle is computed,
            # by summing the (x,y,z) components of all the forces acting on the particle
            forces[j] = np.sum(all_forces[j], axis=0)

        return forces, potential

