from ray_tracing import *


def plotting(ax, barotropic, theta, coordinates=[(x * 1e3, 0) for x in np.arange(-24,20,8)]):
    """
    this functions runs and plots ray tracing simulation for a subplot.
    subplot is given via argument 'ax'.
    """
    # simulation
    ray_system = RayTracing(barotropic=barotropic, theta=theta)
    legend = []
    for x0, z0 in coordinates:  # loop over multiple rays
        ray_system.run_simulation(x0, z0, 40e3, 1e2, x_lim=[-50e3, 60e3], z_lim=-0.6e3)
        ax.plot(*get_data_from_file('Data.out'), 'b')
        legend.append('x0: {}km'.format(x0 / 1000))

    # plot background flow
    x = np.linspace(-20, 20, 50)
    z = np.linspace(-0.5, 0, 50)
    X, Z = np.meshgrid(x, z)
    if barotropic:
        V = ray_system.BF.V(X * 1000)
    else:
        V = ray_system.BF.V(X * 1000, Z * 1000)
    ax.contour(X, Z, V, 3, colors='k')

    ax.set_ylim((-0.3, 0))
    ax.set_xlim((-40, 40))
    title = ['Barotropic' if barotropic else 'Baroclinic'][0] + ' flow, angle: ' + str(theta)
    ax.set_title(title, fontsize=subtitle_size)
    ax.set_xlabel('x [km]')
    ax.set_ylabel('z [km]')
    ax.grid()

# subplots parameters
title_size = 14
subtitle_size = 12
plt_axes = [3, 2]
fig, (axes) = plt.subplots(*plt_axes)
fig.suptitle('RAY TRACING', fontsize=title_size)
fig.tight_layout(rect=[0, 0, 1, 0.93])

# plots
i, j = 0, 0
plotting(axes[i][j], barotropic=True, theta=0)

i, j = 0, 1
plotting(axes[i][j], barotropic=True, theta=180)

i, j = 1, 0
plotting(axes[i][j], barotropic=False, theta=0)

i, j = 1, 1
plotting(axes[i][j], barotropic=False, theta=180)

i, j = 2, 0
plotting(axes[i][j], barotropic=False, theta=45)

i, j = 2, 1
plotting(axes[i][j], barotropic=False, theta=-45)

plt.show()
