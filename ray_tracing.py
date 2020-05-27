import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt


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


class Ray(object):
    """
    contains ray attributes.
    """

    def __init__(self, option):
        if option == 1:
            self.BF = BarotropicBackFlow()
            self.Barotropic = True

    def Cg_x(self, x, kx, kz):
        if self.Barotropic:
            return self.BF.N ** 2 * kx / self.BF.f / kz ** 2
        else:
            return self.BF.N ** 2 * kx / self.BF.f / kz ** 2

    def Cg_z(self, x, kx, kz):
        """
        why is sign(Cg_z)=-sign(kz)??
        """
        if self.Barotropic:
            return -self.BF.N ** 2 * kx ** 2 / self.BF.f / kz ** 3
        else:
            return -self.BF.N ** 2 * kx ** 2 / self.BF.f / kz ** 3 + kx / kz ** 2 * self.BF.dV_dz(x)

    def ode_kx(self, x, kx, kz):
        if self.Barotropic:
            return -0.5 * self.BF.d2V_dx2(x)
        else:
            return -0.5 * self.BF.d2V_dx2(x) \
                   + kx / kz * self.BF.d2V_dxdz(x) \
                   - kx ** 2 / self.BF.f / kz ** 2 * self.BF.N(x) * self.BF.dN_dx(x) \
                   - kx * self.BF.dV_dx(x)  # only if doppler shift k*V is relevant

    def ode_kz(self, x, kx, kz):
        if self.Barotropic:
            return 0
        else:
            return -0.5 * self.BF.d2V_dxdz(x) \
                   + kx / kz * self.BF.d2V_dz2(x) \
                   - kx ** 2 / self.BF.f / kz ** 2 * self.BF.N(x) * self.BF.dN_dz(x) \
                   - kz * self.BF.dV_dz(x)  # only if doppler shift k*V is relevant

    def omega(self, x, kx, kz):
        # make sure it is constant
        self.BF.f_eff(x) + self.BF.N ** 2 * kx ** 2 / 2 / self.BF.f / kz ** 2


class RayTracing(Ray):
    # simulates a coupled pendulums system, including pendulums and springs
    def __init__(self, t0=0, t_final=24 * 30 * 12, N_eval=201, option=1, file_name='Data.out'):
        Ray.__init__(self, option=option)
        self.file_name = file_name

        self.t0 = t0 * 3600
        self.t_final = t_final * 3600
        self.t_eval = np.linspace(self.t0, self.t_final, N_eval)

    def differential_eq(self, t, var0):
        """
        this equation returns differential equations
        """
        if self.Barotropic:
            return self.Cg_x(var0[0], var0[2], var0[3]), self.Cg_z(var0[0], var0[2], var0[3]), \
                   self.ode_kx(var0[0], var0[2], var0[3]), self.ode_kz(var0[0], var0[2], var0[3])
        else:
            return self.Cg_x(*var0), self.Cg_z(*var0), self.ode_kx(*var0), self.ode_kz(*var0)

    def output_line(self, t, x, z, kx, kz, omega):
        output_line = '{0:6.6f}\t{1:6.6f}\t{2:6.6f}\t{3:6.6f}\t{4:6.6f}\n'.format(t, x, z, kx, kz, omega)
        return output_line

    def open_data_file(self):
        file_output = open(self.file_name, 'w')
        file_output.write('time\tx\tz\tk_x\tk_z\tomega\n')
        return file_output

    def main(self, x0, z0, lx0, lz0, x_lim=[-30e3, 34e3], z_lim=-0.5e3):
        """
        uses integrate.solve_ivp, with method='RK45' to solve the differential equations.
        """
        kx0, kz0 = 2 * np.pi / lx0, 2 * np.pi / lz0

        def z_wall(t, y): return y[1] - z_lim
        z_wall.terminal = True
        def x_wall1(t, y): return y[0] - x_lim[0]
        x_wall1.terminal = True
        def x_wall2(t, y): return y[0] - x_lim[1]
        x_wall2.terminal = True

        sol = integrate.solve_ivp(fun=self.differential_eq, t_span=(self.t0, self.t_final), y0=(x0, z0, kx0, kz0),
                                  method='RK45', t_eval=self.t_eval, events=[x_wall1, x_wall2, z_wall])

        # write to output file
        file_output = self.open_data_file()
        for i in range(len(sol.t)):
            x, z, kx, kz = sol.y[:, i]
            file_output.write(self.output_line(sol.t[i], x, z, kx, kz, self.omega(x, kx, kz)))
        file_output.close()


def plot_simulation(file_name):
    """
    plots data from file
    """
    with open(file_name, 'r') as f:
        next(iter(f))
        x = []
        z = []
        for line in f:
            line = line.split()
            x.append(float(line[1]) / 1000)
            z.append(float(line[2]) / 1000)

    plt.plot(x, z)


if __name__ == "__main__":
    ray_system = RayTracing()

    for x0, z0 in [(-20e3, 0), (-10e3, 0), (0, 0), (5e3, 0), (10e3, 0)]:
        ray_system.main(x0, z0, 40e3, 1e2)
        plot_simulation('Data.out')

    # plot background flow
    x = np.linspace(-12, 12, 50)
    z = np.linspace(-0.5, 0, 50)
    X, Z = np.meshgrid(x, z)
    V = ray_system.BF.V(X / 1000)
    plt.contour(X, Z, V, 4, colors='black')

    plt.title('ray tracing')
    plt.ylabel('z [km]')
    plt.xlabel('x [km]')
    plt.grid()
    plt.show()
