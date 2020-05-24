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

    def __init__(self):
        option = 1
        if option == 1:
            self.BF = BarotropicBackFlow()
            self.Barotropic = True

    def calc_Cg_x(self, x, kx, kz):
        if self.Barotropic:
            return self.BF.N ** 2 * kx / self.BF.f_eff(x) / kz ** 2
        else:
            return self.BF.N ** 2 * kx / self.BF.f_eff(x) / kz ** 2

    def calc_Cg_z(self, x, kx, kz):
        if self.Barotropic:
            return -self.BF.N ** 2 * kx ** 2 / self.BF.f_eff(x) / kz ** 3
        else:
            return -self.BF.N ** 2 * kx ** 2 / self.BF.f_eff(x) / kz ** 3 + kx / kz ** 2 * self.BF.dVx_dz(x)

    def ode_x(self, x, kx, kz):
        if self.Barotropic:
            return self.BF.Vx(x) + self.calc_Cg_x(x, kx, kz)
        else:
            return self.BF.Vx(x) + self.calc_Cg_x(x, kx, kz)

    def ode_z(self, x, kx, kz):
        if self.Barotropic:
            return self.calc_Cg_z(x, kx, kz)
        else:
            return self.calc_Cg_z(x, kx, kz)

    def ode_kx(self, x, kx, kz):
        if self.Barotropic:
            return -0.5 * self.BF.d2Vx_dx2(x) - kx * self.BF.dVx_dx(x)  # only if doppler shift k*V is relevant
        else:
            return -0.5 * self.BF.d2Vx_dx2(x) \
                   + kx / kz * self.BF.d2Vx_dxdz(x) \
                   - kx ** 2 / self.BF.f(x) / kz ** 2 * self.BF.N(x) * self.BF.dN_dx(x) \
                   - kx * self.BF.dVx_dx(x)  # only if doppler shift k*V is relevant

    def ode_kz(self, x, kx, kz):
        if self.Barotropic:
            return 0
        else:
            return -0.5 * self.BF.d2Vx_dxdz(x) \
                   + kx / kz * self.BF.d2Vx_dz2(x) \
                   - kx ** 2 / self.BF.f / kz ** 2 * self.BF.N(x) * self.BF.dN_dz(x) \
                   - kz * self.BF.dVx_dz(x)  # only if doppler shift k*V is relevant


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
