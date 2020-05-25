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

    def __init__(self, option):
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


class RayTracing(Ray):
    # simulates a coupled pendulums system, including pendulums and springs
    def __init__(self, x0, z0, kx0, kz0, t0=0, t_final=1, N_eval=11, option=1):
        # initialize pendulums system data. accepts input file name or path (string)
        # initialize ray data.
        Ray.__init__(self, option=option)
        self.x0, self.z0, self.kx0, self.kz0 = x0, z0, kx0, kz0
        self.t0 = t0
        self.t_final = t_final
        self.t_eval = np.linspace(t0, t_final, N_eval)

    def progressing_eq(self, t, var0):
        """
        this equation returns differential equations, Runge-Kutta use.
        """
        if self.Barotropic:
            return self.ode_x(var0[0], var0[2], var0[3]), self.ode_z(var0[0], var0[2], var0[3]), \
                   self.ode_kx(var0[0], var0[2], var0[3]), self.ode_kz(var0[0], var0[2], var0[3])
        else:
            return self.ode_x(*var0), self.ode_z(*var0), self.ode_kx(*var0), self.ode_kz(*var0)

    def output_line(self, t, x, z, kx, kz):
        """
        write data to output file
        """
        output_line = '\n{0:6.6f}\t{1:6.6f}  {2:6.6f}\t{3:6.6f}  {4:6.6f}'.format(t, x, z, kx, kz)
        return output_line

    def open_data_file(self):
        file_output = open('Data.out', 'w')
        file_output.write('\ntime\tx\tz\tk_x\tk_z')
        return file_output

    def main(self):
        """
        the main function runs the algorithm in a loop.
        it writes data to the output file N times.
        """
        # open output file

        sol = integrate.solve_ivp(fun=self.progressing_eq, t_span=(self.t0, self.t_final),
                                  y0=(self.x0, self.z0, self.kx0, self.kz0), method='RK45', t_eval=self.t_eval)

        file_output = self.open_data_file()
        for i in range(len(self.t_eval)):
            file_output.write(self.output_line(sol.t[i], *sol.y[:, i]))
        file_output.close()


ray_system = RayTracing(0, 0, 1, 1)
ray_system.main()
