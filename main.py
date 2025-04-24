from matplotlib.pyplot import figure, show
from numpy import *
from poynting import *

c = 1.0


def plot_retarded_time_condition(omega, include_solution=False, comment=str()):
    particle = ParticleTrajectory(omega=omega)
    t = 0.7
    r = FourVector(t, 0.0, 0.0, 2.0)
    tprimes = linspace(-3.5 * r.spatial_norm / c, t, 100)
    rprimes = [particle.position(tp) for tp in tprimes]
    deltas = [(r.t - rp.t) - spatial_norm(r - rp) for rp in rprimes]
    tr = find_retarded_time(particle, r)
    fig = figure(figsize=[8, 6])
    ax1 = fig.add_subplot(111)
    ax1.plot(tprimes, deltas)
    if include_solution:
        ax1.axvline(tr, ls="--", c="k")
    ax1.axvline(t, ls=":", c="purple")
    ax1.text(t, 0.5, "Now", rotation=90)
    ax1.axhline(0.0, ls="-", c="grey")
    ax1.set_xlabel(r"$t'$")
    ax1.set_ylabel(r"$(t - t') - |\vec x - \vec x'| / c$")
    ax1.set_title(
        rf"$\omega = {omega} c / \ell$ $\vec r = (0, 0, {r[2]}\ell)$ $t = {t}s$"
    )
    ax1.set_ylim(-0.66 * r.spatial_norm / c, 0.66 * r.spatial_norm / c)
    fig.tight_layout(pad=0.1)
    return fig


def plot_retarded_time_in_xy_plane():
    particle = ParticleTrajectory(omega=1.0)
    t = 0.0
    x = linspace(-12.0, 12.0, 300)
    y = linspace(-12.0, 12.0, 300)
    x, y = meshgrid(x, y)
    r = vectorize(lambda x, y: FourVector(t, x, y, 0).spatial_norm)(x, y)
    tr = vectorize(lambda x, y: find_retarded_time(particle, FourVector(t, x, y, 0)))(x, y)
    fig = figure(figsize=[8, 8])
    ax1 = fig.add_subplot(111)
    ax1.set_aspect("equal")
    cim = ax1.pcolormesh(x, y, (1 + c * tr / r) * r**2)
    ax1.set_xlabel(r"$x/\ell$")
    ax1.set_ylabel(r"$y/\ell$")
    fig.tight_layout(pad=0.1)
    return fig


def plot_retarded_time_in_xz_plane():
    particle = ParticleTrajectory(omega=1.0)
    t = 0.0
    x = linspace(-24.0, 24.0, 300)
    z = linspace(-24.0, 24.0, 300)
    x, z = meshgrid(x, z)
    r = vectorize(lambda x, z: FourVector(t, x, 0, z).spatial_norm)(x, z)
    tr = vectorize(lambda x, z: find_retarded_time(particle, FourVector(t, x, 0, z)))(x, z)
    fig = figure(figsize=[8, 8])
    ax1 = fig.add_subplot(111)
    ax1.set_aspect("equal")
    ax1.set_title(r"$\tilde t + r / c$")
    cim = ax1.pcolormesh(x, z, tr + r / c)
    ax1.set_xlabel(r"$x/\ell$")
    ax1.set_ylabel(r"$z/\ell$")
    fig.tight_layout(pad=0.1)
    fig.colorbar(cim)
    return fig


def main():
    plot_retarded_time_condition(
        0.9,
        include_solution=True,
        comment="Unique solution",
    )
    plot_retarded_time_condition(
        1.2,
        include_solution=False,
        comment="Multiple solutions",
    )
    plot_retarded_time_in_xy_plane()
    plot_retarded_time_in_xz_plane()
    show()

main()
