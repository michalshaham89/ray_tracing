import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

"""
this code runs a simulation for ray tracing, according to the data in Kunze 1985.
fits midlatitude open ocean, southward geostrophic jet.
simulation parameters:
2 option for background flows: barotropic and baroclinic.
angle between wave and flow. 
starting coordinate.
(set the parameters in the 'main' part)
"""


class BarotropicBackFlow(object):
    """
    back flow data
    """

    def __init__(self):
        self.f = 1e-4  # coriolis frequency
        self.N = 0.01

    def f_eff(self, x):
        return self.f + 0.5 * self.dV_dx(x)

    @staticmethod
    def V(x):
        return -0.2 * np.exp(-1e-8 * x ** 2)

    @staticmethod
    def dV_dx(x):
        return 4e-9 * x * np.exp(-1e-8 * x ** 2)

    @staticmethod
    def d2V_dx2(x):
        return (-8e-17 * x ** 2 + 4e-9) * np.exp(-1e-8 * x ** 2)


class BaroclinicBackFlow(object):
    """
    back flow data
    ** neglecting dN/dx (N is not dependant on x)
    """

    def __init__(self):
        self.f = 1e-4  # coriolis frequency

    def f_eff(self, x, z):
        return self.f + 0.5 * self.dV_dx(x, z)

    @staticmethod
    def N(z):
        return 0.01 * np.exp(0.002 * z)

    @staticmethod
    def dN_dz(z):
        return 2e-5 * np.exp(0.002 * z)

    @staticmethod
    def V(x, z):
        return -0.2 * np.exp(-1e-8 * x ** 2) * np.exp(0.005 * z)

    @staticmethod
    def dV_dx(x, z):
        return 4e-9 * x * np.exp(-1e-8 * x ** 2) * np.exp(0.005 * z)

    @staticmethod
    def d2V_dx2(x, z):
        return (-8e-17 * x ** 2 + 4e-9) * np.exp(-1e-8 * x ** 2) * np.exp(0.005 * z)

    @staticmethod
    def dV_dz(x, z):
        return -0.001 * np.exp(-1e-8 * x ** 2) * np.exp(0.005 * z)

    @staticmethod
    def d2V_dz2(x, z):
        return -5e-6 * np.exp(-1e-8 * x ** 2) * np.exp(0.005 * z)

    @staticmethod
    def d2V_dxdz(x, z):
        return 2e-11 * x * np.exp(-1e-8 * x ** 2) * np.exp(0.005 * z)


class Ray(object):
    """
    contains ray attributes.
    """

    def __init__(self, barotropic, theta):

        self.theta = np.deg2rad(theta)
        self.DopplerShift = False
        if barotropic:
            self.BF = BarotropicBackFlow()
            self.Barotropic = True
        else:
            self.BF = BaroclinicBackFlow()
            self.Barotropic = False
            if theta != 0:
                self.DopplerShift = True
                self.ky = 0

    def Cg_x(self, x, z, kx, kz):
        if self.Barotropic:
            return self.BF.N ** 2 * kx / self.BF.f / kz ** 2
        else:
            return self.BF.N(z) ** 2 * kx / self.BF.f / kz ** 2 - self.BF.dV_dz(x, z) / kz

    def Cg_z(self, x, z, kx, kz):
        """
        why is sign(Cg_z)=-sign(kz)??
        """
        if self.Barotropic:
            return -self.BF.N ** 2 * kx ** 2 / self.BF.f / kz ** 3
        else:
            return -self.BF.N(z) ** 2 * kx ** 2 / self.BF.f / kz ** 3 + kx / kz ** 2 * self.BF.dV_dz(x, z)

    def ode_kx(self, x, z, kx, kz):
        if self.Barotropic:
            return -0.5 * self.BF.d2V_dx2(x)
        elif not self.DopplerShift:
            return -0.5 * self.BF.d2V_dx2(x, z) \
                   + kx / kz * self.BF.d2V_dxdz(x, z)
        else:
            return -0.5 * self.BF.d2V_dx2(x, z) \
                   + kx / kz * self.BF.d2V_dxdz(x, z) \
                   - self.ky * self.BF.dV_dx(x, z)

    def ode_kz(self, x, z, kx, kz):
        if self.Barotropic:
            return 0
        elif not self.DopplerShift:
            return -0.5 * self.BF.d2V_dxdz(x, z) \
                   + kx / kz * self.BF.d2V_dz2(x, z) \
                   - kx ** 2 / self.BF.f / kz ** 2 * self.BF.N(z) * self.BF.dN_dz(z)
        else:
            return -0.5 * self.BF.d2V_dxdz(x, z) \
                   + kx / kz * self.BF.d2V_dz2(x, z) \
                   - kx ** 2 / self.BF.f / kz ** 2 * self.BF.N(z) * self.BF.dN_dz(z) \
                   - self.ky * self.BF.dV_dz(x, z)

    def omega(self, x, z, kx, kz):
        # make sure it is constant
        if self.Barotropic:
            self.BF.f_eff(x) + self.BF.N ** 2 * kx ** 2 / 2 / self.BF.f / kz ** 2
        elif not self.DopplerShift:
            self.BF.f_eff(x, z) + self.BF.N(z) ** 2 * kx ** 2 / 2 / self.BF.f / kz ** 2
        else:
            self.BF.f_eff(x, z) + self.BF.N(z) ** 2 * kx ** 2 / 2 / self.BF.f / kz ** 2 \
            + self.ky * self.BF.V(x, z)


class RayTracing(Ray):
    """
    ray tracing simulation using Runge-Kutta 4order algorithm.
    the code writes the output data to a text file.
    """

    def __init__(self, t0=0, t_final=24 * 30 * 12, N_eval=201, barotropic=True, theta=0, file_name='Data.out'):
        Ray.__init__(self, barotropic=barotropic, theta=theta)
        self.file_name = file_name
        self.t0 = t0 * 3600
        self.t_final = t_final * 3600
        self.t_eval = np.linspace(self.t0, self.t_final, N_eval)

    def differential_eq(self, t, var0):
        """
        this equation returns the differential equations for the simulation
        """
        return self.Cg_x(*var0), self.Cg_z(*var0), self.ode_kx(*var0), self.ode_kz(*var0)

    def output_line(self, t, x, z, kx, kz, omega):
        output_line = '{0:6.6f}\t{1:6.6f}\t{2:6.6f}\t{3:6.6f}\t{4:6.6f}\n'.format(t, x, z, kx, kz, omega)
        return output_line

    def open_data_file(self):
        file_output = open(self.file_name, 'w')
        file_output.write('time\tx\tz\tk_x\tk_z\tomega\n')
        return file_output

    def run_simulation(self, x0, z0, Lx0, Lz0, x_lim=[-30e3, 34e3], z_lim=-0.5e3):
        """
        uses integrate.solve_ivp, with method='RK45' to solve the differential equations.
        """
        kx0 = 2 * np.pi / Lx0 * np.cos(self.theta)
        self.ky = 2 * np.pi / Lx0 * np.sin(self.theta)
        kz0 = 2 * np.pi / Lz0

        # defining space limits for the RK45 algorithm
        def z_wall(t, y): return y[1] - z_lim
        z_wall.terminal = True
        def x_wall1(t, y): return y[0] - x_lim[0]
        x_wall1.terminal = True
        def x_wall2(t, y): return y[0] - x_lim[1]
        x_wall2.terminal = True

        # RK45 algorithm
        sol = integrate.solve_ivp(fun=self.differential_eq, t_span=(self.t0, self.t_final), y0=(x0, z0, kx0, kz0),
                                  method='RK45', t_eval=self.t_eval, events=[x_wall1, x_wall2, z_wall])

        # write to output file
        file_output = self.open_data_file()
        for i in range(len(sol.t)):
            x, z, kx, kz = sol.y[:, i]
            file_output.write(self.output_line(sol.t[i], x, z, kx, kz, self.omega(x, z, kx, kz)))
        file_output.close()


def get_data_from_file(file_name):
    with open(file_name, 'r') as f:
        next(iter(f))
        x = []
        z = []
        for line in f:
            line = line.split()
            x.append(float(line[1]) / 1000)
            z.append(float(line[2]) / 1000)

    return x, z


if __name__ == "__main__":

    # simulation parameters
    barotropic = True  # background flow
    theta = 0  # degrees
    coordinates = [(x * 1e3, 0) for x in np.linspace(-30, 20, 11)]  # [(x0, z0),]

    # simulation
    ray_system = RayTracing(barotropic=barotropic, theta=theta)
    legend = []
    for x0, z0 in coordinates: # loop over multiple rays
        ray_system.run_simulation(x0, z0, 40e3, 1e2, x_lim=[-50e3, 60e3], z_lim=-0.6e3)
        plt.plot(*get_data_from_file('Data.out'))
        legend.append('x0: {}km'.format(x0 / 1000))

    # plot background flow
    x = np.linspace(-20, 20, 50)
    z = np.linspace(-0.5, 0, 50)
    X, Z = np.meshgrid(x, z)
    if barotropic:
        V = ray_system.BF.V(X * 1000)
    else:
        V = ray_system.BF.V(X * 1000, Z * 1000)
    plt.contour(X, Z, V, 3, colors='black')

    plt.title('Ray Tracing: ' + ['Barotropic' if barotropic else 'Baroclinic'][0] + ' flow, angle: ' + str(theta))
    plt.ylabel('z [km]')
    plt.xlabel('x [km]')
    plt.legend(legend)
    plt.ylim((-0.5, 0))
    plt.xlim((-40, 40))
    plt.grid()
    plt.show()

    """
    problems:
    alot denser in the z axis.
    baroclinic, theta=180, x=0, ray is not trapped.
    baroclinic, theta=0, only 3 (instead of 4) rays are trapped..
    baroclinic, theta=60, x=-25 is trapped, x=15 is not.
    baroclinic, theta=210, x=-15.
    """
