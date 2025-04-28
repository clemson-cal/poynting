from matplotlib.pyplot import *
from matplotlib import rc
from numpy import *
from poynting import *


def plot_poynting_vector_time_series(ax, omega=0.5, radius=10.0, latitude=0.0, normalize=False, **pl_args):
    """
    Plot the time series of the Poynting vector as a function of time.
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axis instance to plot on
    omega : float, optional
        The angular frequency of the particle (default: 0.5)
    radius : float, optional [lambda]
        The radius at which to compute the Poynting vector (default: 100.0)
    latitude : float, optional [deg]
        The latitude at which to compute the Poynting vector (default: 0.0)
    normalize : bool, optional
        Whether to normalize the Poynting vector (default: False)
    **pl_args : dict
        Additional plotting arguments passed to ax.plot()
    """
    radius *= 2 * pi / omega
    particle = ParticleTrajectory(omega=omega)
    t = linspace(0, 2 * pi / omega, 1000)  # time
    x = radius * cos(latitude * pi / 180)
    z = radius * sin(latitude * pi / 180)
    S = array([poynting_vector(particle, t=t, x=x, z=z) for t in t])
    Sr = array([(S.x * x + S.z * z) / radius for S in S])

    if normalize:
        Sr /= max(Sr)

    ax.plot(t * omega, Sr, **pl_args)
    ax.set_xlabel(r'Phase $\omega t$')
    ax.set_ylabel(r'$\mathrm{S}_r$' + r' (normalized)' if normalize else '')
    ax.set_xticks([0, pi/2, pi, 3*pi/2, 2*pi])
    ax.set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
    # ax.grid(True)


def fig_poynting_vector_time_series(fig):
    # Create subplots with zero vertical spacing between them
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    # First subplot (15° latitude)
    plot_poynting_vector_time_series(ax1, omega=0.10, latitude=15, normalize=True, c='k', lw=0.5, label=r"$\omega \ell / c = 0.1$")
    plot_poynting_vector_time_series(ax1, omega=0.40, latitude=15, normalize=True, c='k', lw=1.0, label=r"$\omega \ell / c = 0.4$")
    plot_poynting_vector_time_series(ax1, omega=0.80, latitude=15, normalize=True, c='k', lw=1.5, label=r"$\omega \ell / c = 0.8$")
    ax1.text(0.05, 0.85, r"Latitude 15°", transform=ax1.transAxes)
    ax1.set_xticklabels([])
    ax1.set_xlabel('')
    ax1.set_xticklabels([])

    # Second subplot (45° latitude)
    plot_poynting_vector_time_series(ax2, omega=0.10, latitude=45, normalize=True, c='k', lw=0.5, label=r"$\omega \ell / c = 0.1$")
    plot_poynting_vector_time_series(ax2, omega=0.40, latitude=45, normalize=True, c='k', lw=1.0, label=r"$\omega \ell / c = 0.4$")
    plot_poynting_vector_time_series(ax2, omega=0.80, latitude=45, normalize=True, c='k', lw=1.5, label=r"$\omega \ell / c = 0.8$")
    ax2.text(0.05, 0.85, r"Latitude 45°", transform=ax2.transAxes)
    ax2.set_xticklabels([])
    ax2.set_xlabel('')
    ax2.set_xticklabels([])

    # Third subplot (75° latitude)
    plot_poynting_vector_time_series(ax3, omega=0.10, latitude=75, normalize=True, c='k', lw=0.5, label=r"$\omega \ell / c = 0.1$")
    plot_poynting_vector_time_series(ax3, omega=0.40, latitude=75, normalize=True, c='k', lw=1.0, label=r"$\omega \ell / c = 0.4$")
    plot_poynting_vector_time_series(ax3, omega=0.80, latitude=75, normalize=True, c='k', lw=1.5, label=r"$\omega \ell / c = 0.8$")
    ax3.text(0.05, 0.85, r"Latitude 75°", transform=ax3.transAxes)
    ax3.legend(loc='upper right')

    # Adjust layout with no spacing between subplots
    fig.subplots_adjust(hspace=0)
    fig.tight_layout(pad=0.1)
    return fig


def plot_poynting_vector_power_spectrum(ax, omega=0.5, radius=10.0, latitude=0.0, normalize=False, **pl_args):
    """
    Plot the power spectrum of the Poynting vector.
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axis instance to plot on
    omega : float, optional
        The angular frequency of the particle (default: 0.5)
    radius : float, optional [lambda]
        The radius at which to compute the Poynting vector (default: 10.0)
    latitude : float, optional [deg]
        The latitude at which to compute the Poynting vector (default: 0.0)
    normalize : bool, optional
        Whether to normalize the power spectrum (default: False)
    **pl_args : dict
        Additional plotting arguments passed to ax.plot()
    """
    radius *= 2 * pi / omega
    particle = ParticleTrajectory(omega=omega)
    t = linspace(0, 2 * pi / omega, 10000, endpoint=False)  # time
    x = radius * cos(latitude * pi / 180)
    z = radius * sin(latitude * pi / 180)
    S = array([poynting_vector(particle, t=t, x=x, z=z) for t in t])
    Sr = array([(S.x * x + S.z * z) / radius for S in S])

    # Compute FFT
    Sr_fft = fft.rfft(Sr)
    Sr_squared = (Sr_fft * Sr_fft.conj()).real
    freqs = fft.rfftfreq(len(t), t[1] - t[0])

    if normalize:
        Sr_squared /= max(Sr_squared)

    ax.plot(freqs / omega * 2 * pi, Sr_squared, **pl_args)
    ax.set_xlabel(r'Frequency ($\omega$)')
    ax.set_ylabel(r'Power' + r' (normalized)' if normalize else '')
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.grid(True)


def fig_poynting_vector_power_spectrum_multi_latitude(fig):
    # Create subplots with zero vertical spacing between them
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    # First subplot (15° latitude)
    plot_poynting_vector_power_spectrum(ax1, omega=0.010, latitude=15, normalize=True, c='k', lw=1.0, label=r"$\omega \ell / c = 0.01$")
    plot_poynting_vector_power_spectrum(ax1, omega=0.100, latitude=15, normalize=True, c='k', lw=1.0, label=r"$\omega \ell / c = 0.1$")
    plot_poynting_vector_power_spectrum(ax1, omega=0.900, latitude=15, normalize=True, c='k', lw=1.5, label=r"$\omega \ell / c = 0.9$")
    plot_poynting_vector_power_spectrum(ax1, omega=0.990, latitude=15, normalize=True, c='k', lw=2.0, label=r"$\omega \ell / c = 0.99$")
    ax1.text(0.75, 0.85, r"Latitude 15°", transform=ax1.transAxes)
    ax1.set_xticklabels([])
    ax1.set_xlabel('')

    # Second subplot (45° latitude)
    plot_poynting_vector_power_spectrum(ax2, omega=0.010, latitude=45, normalize=True, c='k', lw=1.0, label=r"$\omega \ell / c = 0.01$")
    plot_poynting_vector_power_spectrum(ax2, omega=0.100, latitude=45, normalize=True, c='k', lw=1.0, label=r"$\omega \ell / c = 0.1$")
    plot_poynting_vector_power_spectrum(ax2, omega=0.900, latitude=45, normalize=True, c='k', lw=1.5, label=r"$\omega \ell / c = 0.9$")
    plot_poynting_vector_power_spectrum(ax2, omega=0.990, latitude=45, normalize=True, c='k', lw=2.0, label=r"$\omega \ell / c = 0.99$")
    ax2.text(0.75, 0.85, r"Latitude 45°", transform=ax2.transAxes)
    ax2.set_xticklabels([])
    ax2.set_xlabel('')

    # Third subplot (75° latitude)
    plot_poynting_vector_power_spectrum(ax3, omega=0.010, latitude=75, normalize=True, c='k', lw=1.0, label=r"$\omega \ell / c = 0.01$")
    plot_poynting_vector_power_spectrum(ax3, omega=0.100, latitude=75, normalize=True, c='k', lw=1.0, label=r"$\omega \ell / c = 0.1$")
    plot_poynting_vector_power_spectrum(ax3, omega=0.900, latitude=75, normalize=True, c='k', lw=1.5, label=r"$\omega \ell / c = 0.9$")
    plot_poynting_vector_power_spectrum(ax3, omega=0.990, latitude=75, normalize=True, c='k', lw=2.0, label=r"$\omega \ell / c = 0.99$")

    # Set y-axis limits for all plots
    ax1.set_ylim(2e-10, 0.99)
    ax2.set_ylim(2e-10, 0.99)
    ax3.set_ylim(2e-10, 0.99)

    ax3.text(0.75, 0.85, r"Latitude 75°", transform=ax3.transAxes)
    ax1.legend(loc='lower right')

    fig.subplots_adjust(hspace=0)
    fig.tight_layout(pad=0.1)
    return fig


def fig_poynting_vector_power_spectrum(fig):
    ax1 = fig.add_subplot(1, 1, 1)
    plot_poynting_vector_power_spectrum(ax1, omega=0.10, latitude=45, normalize=True, c='k', lw=0.5, label=r"$\omega \ell / c = 0.1$")
    plot_poynting_vector_power_spectrum(ax1, omega=0.40, latitude=45, normalize=True, c='k', lw=1.0, label=r"$\omega \ell / c = 0.4$")
    plot_poynting_vector_power_spectrum(ax1, omega=0.80, latitude=45, normalize=True, c='k', lw=1.5, label=r"$\omega \ell / c = 0.8$")
    plot_poynting_vector_power_spectrum(ax1, omega=0.99, latitude=45, normalize=True, c='k', lw=1.5, label=r"$\omega \ell / c = 0.99$")
    fig.tight_layout(pad=0.1)
    return fig


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
    particle = ParticleTrajectory(omega=omega)
    theta = linspace(0, pi, 1000)  # angle from z-axis (0 to π)
    x = radius * sin(theta)
    z = radius * cos(theta)
    Sr = array([poynting_vector_time_avg(particle, x=x, z=z) for x, z in zip(x, z)])

    dP_dOmega = Sr * radius**2
    if normalize:
        dP_dOmega /= max(dP_dOmega)

    ax.plot(theta, dP_dOmega, **pl_args)
    ax.set_xlabel(r'Angle $\theta$')
    ax.set_ylabel(r'$\mathrm{d}P/\mathrm{d}\Omega$' + r' (normalized)' if normalize else '')
    ax.set_xticks([0, pi/4, pi/2, 3*pi/4, pi])
    ax.set_xticklabels([r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'])
    ax.grid(True)


def fig_angular_power_distribution(fig):
    ax = fig.add_subplot(111)
    plot_angular_power_distribution(ax, normalize=True, omega=0.10, c='k', lw=0.5, label=r"$\omega \ell / c = 0.1$")
    plot_angular_power_distribution(ax, normalize=True, omega=0.40, c='k', lw=1.0, label=r"$\omega \ell / c = 0.4$")
    plot_angular_power_distribution(ax, normalize=True, omega=0.80, c='k', lw=1.5, label=r"$\omega \ell / c = 0.8$")
    plot_angular_power_distribution(ax, normalize=True, omega=0.99, c='k', lw=2.0, label=r"$\omega \ell / c = 0.99$")
    ax.legend(loc='lower right')
    fig.suptitle(r'Angular Power Distribution as a Function of $\omega \ell / c$')
    fig.tight_layout()
    return fig


def plot_poynting_streamplot(ax, t=0.0, omega=0.1, samples=100, density=10.0):
    """
    Plot streamplot of the Poynting vector field in the x-z plane.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axis instance to plot on
    t : float, optional
        The time at which to compute the Poynting vector (default: 0.0)
    omega : float, optional
        The angular frequency of the particle (default: 0.1)
    samples : int, optional
        Number of sample points in each dimension (default: 100)
    density : float, optional
        Density of streamlines in the plot (default: 10.0)
    """
    particle = ParticleTrajectory(omega=omega)
    d = particle.position(t)
    c = 1.0
    L = 2 * pi * c / particle.omega
    x = linspace(0.2 * L, 10 * L, samples)
    z = linspace(-5 * L, 5 * L, samples)
    X, Z = meshgrid(x, z)

    S = vectorize(lambda x, z: poynting_vector(particle, t, x, 0, z))(X, Z)
    Sx = vectorize(lambda S: S.x)(S)
    Sz = vectorize(lambda S: S.z)(S)
    ax.streamplot(X / L, Z / L, Sx, Sz, density=density, linewidth=0.5, arrowsize=0.5, color='k')
    ax.plot(d.x, d.z, 'bo', markersize=5)
    ax.set_xlabel(r'$x / \lambda$')
    ax.set_ylabel(r'$z / \lambda$')
    ax.set_title(rf'Poynting Vector Flow in the $x$-$z$ Plane at $t={t:0.2f}$')
    ax.set_aspect('equal')


def fig_poynting_streamplot(fig):
    ax = fig.add_subplot(111)
    plot_poynting_streamplot(ax, t=0.0)
    fig.tight_layout()
    return fig


def get_available_figures():
    """Get all functions starting with 'fig_' in the current module"""
    return {name: func for name, func in globals().items() if callable(func) and name.startswith('fig_')}


def main():
    import argparse

    available_figs = get_available_figures()
    fig_names = [name[4:] for name in available_figs.keys()]  # Remove 'fig_' prefix

    parser = argparse.ArgumentParser(description='Demonstrates EM radiation calculations')
    parser.add_argument('figures', nargs='*', choices=fig_names + ['all'])
    parser.add_argument('--save', action='store_true', help='Save figures as PDF')
    args = parser.parse_args()

    # If no figures specified or 'all' is specified, generate all figures
    if 'all' in args.figures:
        selected_figures = list(available_figs.values())
    else:
        selected_figures = [available_figs[f'fig_{name}'] for name in args.figures]

    if args.save:
        rc('text', usetex=True)
        rc('font', family='serif')

    figures = []
    for fig_func in selected_figures:
        fig = fig_func(figure(figsize=(8, 6)))
        figures.append((fig, fig_func.__name__[4:]))  # Store figure with name (minus 'fig_' prefix)

    if args.save:
        for fig, name in figures:
            filename = f"{name}.pdf"
            fig.savefig(filename)
            print(f"Saved figure to {filename}")
    else:
        show()


if __name__ == "__main__":
    main()
