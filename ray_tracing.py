import numpy as np
from scipy import integrate


class BarotropicBackFlow(object):
    """
    back flow data
    """

    def __init__(self):
        self.f = 1e-4  # coriolis frequency
        self.N = 0.01

    def f_eff(self, x):
        return self.f + 0.5 * self.dVx_dx(x)

    @staticmethod
    def Vx(x):
        return -0.2 * np.exp(-1e-8 * x ** 2)

    @staticmethod
    def dVx_dx(x):
        return 4e-9 * x * np.exp(-1e-8 * x ** 2)

    @staticmethod
    def d2Vx_dx2(x):
        return (-8e-17 * x ** 2 + 4e-9) * np.exp(-1e-8 * x ** 2)


class Ray(object):
    """
    contains data about the ray.

    I am not sure if to separate variables and function to x, y, z,
    or to leave it in a 'vector' mode. ( k=[k_x, k_y, k_z] )
    it will be easier to switch to 2d when it is in vector mode,
    but the code will be less 'readable'.
    """

    def __init__(self, x0, z0, kx0, kz0):
        option = 1
        if option == 1:
            self.BF = BarotropicBackFlow()
            self.Barotropic = True

        self.x, self.z, self.kx, self.kz = x0, z0, kx0, kz0

    def calc_Cg_x(self):
        if self.Barotropic:
            return self.BF.N ** 2 * self.kx / self.BF.f_eff(self.x) / self.kz ** 2
        else:
            return self.BF.N ** 2 * self.kx / self.BF.f_eff(self.x, self.z) / self.kz ** 2

    def calc_Cg_z(self):
        if self.Barotropic:
            return -self.BF.N ** 2 * self.kx ** 2 / self.BF.f_eff(self.x) / self.kz ** 3
        else:
            return -self.BF.N ** 2 * self.kx ** 2 / self.BF.f_eff(self.x, self.z) / self.kz ** 3 \
                   + self.kx / self.kz ** 2 * self.BF.dV
            x_dz(self.x0, self.z0)

    def ode_x(self):
        if self.Barotropic:
            return self.BF.Vx(self.x) + self.calc_Cg_x()
        else:
            return self.BF.Vx(self.x) + self.calc_Cg_x()

    def ode_z(self):
        if self.Barotropic:
            return self.calc_Cg_z()
        else:
            return self.calc_Cg_z()

    def ode_kx(self):
        if self.Barotropic:
            return -0.5 * self.BF.d2Vx_dx2(self.x) - self.kx * self.BF.dVx_dx(self.x)  # only if doppler shift k*V is relevant
        else:
            return -0.5 * self.BF.d2Vx_dx2(self.x, self.z) \
                   + self.kx / self.kz * self.BF.d2Vx_dxdz(self.x, self.z) \
                   - self.kx ** 2 / self.BF.f_eff(self.x, self.z) / self.kz ** 2 * self.BF.N(self.x, self.z) \
                   * self.BF.dN_dx(self.x, self.z) \
                   - self.kx * self.BF.dVx_dx(self.x, self.z)  # only if doppler shift k*V is relevant

    def ode_kz(self, x, kx, kz):
        if self.Barotropic:
            return 0
        else:
            return -0.5 * self.BF.d2Vx_dxdz(self.x, self.z) \
                   + self.kx / kz * self.BF.d2Vx_dz2(self.x, self.z) \
                   - self.kx ** 2 / self.BF.f / kz ** 2 * self.BF.N(self.x, self.z) \
                   * self.BF.dN_dz(self.x, self.z) \
                   - kz * self.BF.dVx_dz(self.x, self.z)  # only if doppler shift k*V is relevant


class RayTracing(Ray):
    # simulates a coupled pendulums system, including pendulums and springs
    def __init__(self, x0, z0, kx0, kz0):
        Ray.__init__(self, x0, z0, kx0, kz0)
        # initialize pendulums system data. accepts input file name or path (string)
        # initialize ray data.
        self.t = 0
        self.t_final = 1
        N_measure = 11
        self.dt = int((self.t_final - self.t) / N_measure)

    def jump(self, t0, t_bound):
        """
        jump the system one time leap.
        """
        self.x, self.z, self.kx, self.kz = integrate.RK45((self.ode_x, self.ode_z, self.ode_kx, self.ode_kz),
                                                          t0, (self.x, self.z, self.kx, self.kz), t_bound)

    def output_line(self):
        """
        write data to output file
        """
        output_line = '\n{0:6.6f}\t{1:6.6f}  {2:6.6f}\t{3:6.6f}  {4:6.6f}'\
            .format(self.t, self.x, self.z, self.kx, self.kz)
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
        file_output.write(self.output_line())

        # start moving time forward using algorithm
        while self.t <= self.t_final:
            self.jump(self.t, self.t + self.dt)
            self.t += self.dt  # before or after the step?
            file_output.write(self.output_line())

        # close output files
        file_output.close()


ray_system = RayTracing(0, 0, 0, 0)
ray_system.main()
