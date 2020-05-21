import numpy as np
from scipy import integrate


class BackFlow(object):
    """
    back flow data
    """

    def __init__(self):
        pass

    # initialize flow data. (we start at steady flow)

    @staticmethod
    def N(x, z):
        pass

    @staticmethod
    def dN_dx(x, z):
        pass

    @staticmethod
    def dN_dz(x, z):
        pass

    @staticmethod
    def f(x, z):
        pass

    @staticmethod
    def Vx(x, z):
        pass

    @staticmethod
    def dVx_dx(x, z):
        pass

    @staticmethod
    def dVx_dz(x, z):
        pass

    @staticmethod
    def d2Vx_dz2(x, z):
        pass

    @staticmethod
    def d2Vx_dx2(x, z):
        pass

    @staticmethod
    def d2Vx_dxdz(x, z):
        pass


class Ray(BackFlow):
    """
    contains data about the ray.

    I am not sure if to separate variables and function to x, y, z,
    or to leave it in a 'vector' mode. ( k=[k_x, k_y, k_z] )
    it will be easier to switch to 2d when it is in vector mode,
    but the code will be less 'readable'.
    """

    def __init__(self):
        BackFlow.__init__()

    def calc_Cg_x(self, x, z, kx, kz):
        return self.N(x, z) ** 2 * kx / self.f(x, z) / kz ** 2

    def calc_Cg_z(self, x, z, kx, kz):
        return -self.N(x, z) ** 2 * kx ** 2 / self.f(x, z) / kz ** 3 + kx / kz ** 2 * self.dVx_dz(x, z)

    def ode_x(self, x, z):
        return BackFlow.Vx(x, z) + self.calc_Cg_x

    def ode_z(self):
        return self.calc_Cg_z

    def ode_kx(self, x, z, kx, kz):
        return -0.5 * self.d2Vx_dx2(x, z) \
               + kx / kz * self.d2Vx_dxdz(x, z) \
               - kx ** 2 / self.f / kz ** 2 * self.N * self.dN_dx(x, z) \
               - kx * self.dVx_dx(x, z)  # only if doppler shift k*V is relevant

    def ode_kz(self, x, z, kx, kz):
        return -0.5 * self.d2Vx_dxdz(x, z) \
               + kx / kz * self.d2Vx_dz2(x, z) \
               - kx ** 2 / self.f / kz ** 2 * self.N * self.dN_dz(x, z) \
               - kz * self.dVx_dz(x, z)  # only if doppler shift k*V is relevant


class RayTracing(object):
    # simulates a coupled pendulums system, including pendulums and springs
    def __init__(self, x0, z0, kx0, kz0):
        # initialize pendulums system data. accepts input file name or path (string)
        # initialize ray data.
        self.x0, self.z0, self.kx0, self.kz0 = x0, z0, kx0, kz0
        self.ray = Ray()
        self.t = 0
        self.t_final = 1
        N_measure = 11
        self.dt = int((self.t_final - self.t) / N_measure)

    def set_system(self):
        """
        set the system to initial time.
        this function is still not needed. maybe later on
        """
        pass

    def progressing_eq(self):
        """
        this equation returns differential equations, Runge-Kutta use.
        """
        dx_dt = 0
        dz_dt = 0
        dkx_dt = 0
        dkz_dt = 0
        return dx_dt, dz_dt, dkx_dt, dkz_dt

    def jump(self, t0, t_bound, y0):
        """
        jump the system one time leap.
        """
        y = integrate.RK45(self.progressing_eq, t0, y0, t_bound)
        return y

    def output_line(self, x, z, kx, kz):
        """
        write data to output file
        """
        output_line = '\n{0:6.6f}\t{1:6.6f}  {2:6.6f}\t{4:6.6f}  {5:6.6f}\t{6:6.6f}\t{7:6.6f}  {8:6.6f}  {9:6.6f}' \
            .format(self.t, x, z, kx, kz, self.ray.omega, *self.ray.Cg)
        return output_line

    def open_data_file(self):
        file_output = open('Data.out', 'w')
        file_output.write('\ntime\tx\tz\tk_x\tk_z\tomega\tCg_x\tCg_y\tCg_z\n')
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
        x, z, kx, kz = self.x0, self.z0, self.kx0, self.kz0
        file_output.write(self.output_line(x, z, kx, kz))

        # start moving time forward using algorithm
        while self.t <= self.t_final:
            x, z, kx, kz = self.jump(self.t, self.t + self.dt, x, z, kx, kz)
            self.t += self.dt  # before or after the step?
            file_output.write(self.output_line(x, z, kx, kz))

        # close output files
        file_output.close()


ray_system = RayTracing()
ray_system.main()
