from math import sqrt, pi
from .ext import _core

__all__ = ['poynting_vector', 'time_averaged_energy_flux']


def poynting_vector(particle, t=0.0, x=0.0, y=0.0, z=0.0, fd_order=4, fd_eps=1e-3):
    """
    Calculate the Poynting vector at a specified spacetime position for a given particle.

    Parameters
    ----------
    particle : object
        Particle object containing electromagnetic field information.
    t : float, optional
        Time coordinate. Default is 0.0.
    x : float, optional
        X spatial coordinate. Default is 0.0.
    y : float, optional
        Y spatial coordinate. Default is 0.0.
    z : float, optional
        Z spatial coordinate. Default is 0.0.
    fd_order : int, optional
        Order of finite difference approximation. Default is 4.
    fd_eps : float, optional
        Step size for finite difference. Default is 1e-3.

    Returns
    -------
    object
        Poynting vector at the specified position and time.
    """
    F = _core.faraday_tensor(particle, _core.FourVector(t, x, y, z), fd_order, fd_eps)
    return _core.poynting_vector(F)


def time_averaged_energy_flux(particle, x=0.0, y=0.0, z=0.0, nhat='rhat'):
    """
    Calculate the time-averaged energy flux through a surface element.

    Parameters
    ----------
    particle : object
        Particle object containing electromagnetic field information.
    x : float, optional
        X spatial coordinate. Default is 0.0.
    y : float, optional
        Y spatial coordinate. Default is 0.0.
    z : float, optional
        Z spatial coordinate. Default is 0.0.
    nhat : str or tuple, optional
        Direction of the normal vector to the surface. If 'rhat', uses the
        radial unit vector. Otherwise, expected to be a 3-tuple representing
        the normal vector components. Default is 'rhat'.

    Returns
    -------
    float
        Time-averaged energy flux through the surface element.
    """
    from scipy.integrate import quad

    if nhat == 'rhat':
        r = sqrt(x**2 + y**2 + z**2)
        nhat = (x / r, y / r, z / r)

    def integrand(t):
        S = poynting_vector(particle, t, x, y, z)
        return S.x * nhat[0] + S.y * nhat[1] + S.z * nhat[2]

    period = 2 * pi / particle.omega
    result, _ = quad(integrand, 0.0, period)
    return result / period
