from matplotlib.pyplot import *
from numpy import *
from poynting import *
from matplotlib import rc

# Enable LaTeX rendering
rc('text', usetex=True)
rc('font', family='serif')


def plot_angular_power_distribution(ax, omega=0.5, radius=100.0, normalize=False, **pl_args):
    """
    Plot the angular power distribution dP/dΩ (which is <Sr> * r^2) as a function of theta.
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axis instance to plot on
    omega : float, optional
        The angular frequency of the particle (default: 0.5)
    radius : float, optional
        The radius at which to compute the Poynting vector (default: 100.0)
    normalize : bool, optional
        Whether to normalize the power distribution to its maximum value (default: False)
    **pl_args : dict
        Additional plotting arguments passed to ax.plot()
    """
    particle = ParticleTrajectory(omega=omega, ell=1.0)
    theta = linspace(0, pi, 1000)  # angle from z-axis (0 to π)
    x = radius * sin(theta)
    z = radius * cos(theta)
    Sr = array([poynting_vector_time_avg(particle, x=x, z=z) for x, z in zip(x, z)])

    # Calculate dP/dΩ = <Sr> * r^2
    dP_dOmega = Sr * radius**2
    if normalize:
        dP_dOmega /= max(dP_dOmega)

    # Plot the angular power distribution
    line, = ax.plot(theta, dP_dOmega, **pl_args)

    ax.set_xlabel(r'Angle $\theta$')
    ax.set_ylabel(r'$\mathrm{d}P/\mathrm{d}\Omega$' + r' (normalized)' if normalize else '')
    ax.set_xticks([0, pi/4, pi/2, 3*pi/4, pi])
    ax.set_xticklabels([r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'])
    ax.grid(True)


def plot_time_avg_poynting_vector_meridion():
    fig = figure(figsize=[8, 6])
    ax = fig.add_subplot(111)
    plot_angular_power_distribution(ax, normalize=True, omega=0.10, c='k', lw=0.5, label=r"$\omega \ell / c = 0.1$")
    plot_angular_power_distribution(ax, normalize=True, omega=0.40, c='k', lw=1.0, label=r"$\omega \ell / c = 0.4$")
    plot_angular_power_distribution(ax, normalize=True, omega=0.80, c='k', lw=1.5, label=r"$\omega \ell / c = 0.8$")
    plot_angular_power_distribution(ax, normalize=True, omega=0.99, c='k', lw=2.0, label=r"$\omega \ell / c = 0.99$")
    ax.legend(loc='lower right')
    fig.suptitle(r'Angular Power Distribution as a Function of $\omega \ell / c$')
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    plot_time_avg_poynting_vector_meridion()
    # show()
    savefig('dP_domega.pdf')
