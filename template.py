import numpy as np

class Ray(object):

    """
    contains data about the ray.

    I am not sure if to separate variables and function to x, y, z,
    or to leave it in a 'vector' mode. ( k=[k_x, k_y, k_z] )
    it will be easier to switch to 2d when it is in vector mode,
    but the code will be less 'readable'.
    """

    def __init__(self, k_x, k_y, k_z, omega, x, y, z):
        # initialize ray data.
        self.k = [k_x, k_y, k_z]
        self.omega = omega
        self.r = [x, y, z]
        self.Cg = [0,0,0]

    def calc_Cg(self):
        """
        :return: Cg_x, Cg_Y, Cg_z
        """
        pass

    def calc_omega(self):
        """
        :return: omega
        """
        pass

    def calc_k(self):
        """
        :return: k_x, k_y, k_z
        """
        pass

    def calc_r(self):
        """
        :return: x, y, z
        """
        pass


class BackFlow(object):
    """
    back flow data
    """
    def __init__(self, N, f):
        # initialize flow data. (we start at steady flow)
        self.N = N
        self.f = f

    def V(self):
        pass


class RayTracing(object):
    # simulates a coupled pendulums system, including pendulums and springs
    def __init__(self):
        # initialize pendulums system data. accepts input file name or path (string)
        self.ray = Ray(0,0,0,0,0,0,0)
        self.back_flow = BackFlow(0,0)
        self.t = 0
        self.dt = 0.01
        self.t_final = 1
        N_measure = 10
        self.dt_measure = int((self.t_final-self.t)/self.dt/N_measure)*self.dt

    def set_system(self):
        # set the system to initial time
        pass

    def jump(self):
        """
        jump the system one time leap.
        """
        pass

    def output_line(self):
        """
        write data to output file
        """
        output_line = '\n{0:6.6f}\t{1:6.6f}  {2:6.6f}  {3:6.6f}\t{4:6.6f}  {5:6.6f}  {6:6.6f}\t{7:6.6f}\t{8:6.6f}' \
                      '  {9:6.6f}  {10:6.6f}'.format(self.t, *self.ray.r, *self.ray.k, self.ray.omega, *self.ray.Cg)
        return output_line

    def open_data_file(self):
        file_output = open('Data.out', 'w')
        file_output.write('\ntime\tx\ty\tz\tk_x\tk_y\tk_z\tomega\tCg_x\tCg_y\tCg_z\n')
        return file_output

    def main(self):
        """
        the main function runs the algorithm in a loop.
        it writes data to the output file N times.
        """
        # open output file
        file_output = self.open_data_file()
        # set system for simulation and output initial data
        self.set_system()
        file_output.write(self.output_line())

        # start moving time forward using algorithm
        while self.t <= self.t_final:
            self.jump()
            if self.t % self.dt_measure == 0:
                file_output.write(self.output_line())
            self.t += self.dt # before or after the step?
        # close output files
        file_output.close()


ray_system = RayTracing()
ray_system.main()
