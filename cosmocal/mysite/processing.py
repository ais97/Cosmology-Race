import numpy as np
from astropy import units as u
from astropy import constants as const
from scipy import integrate


def age_of_universe(z, H0, omega_m0, omega_rad0, omega_lam0, omega_k0):

    def humanyears(B):
        B = float(B)
        KB = float(1000)
        MB = float(KB ** 2)
        GB = float(KB ** 3)
        TB = float(KB ** 4)

        if B < KB:
          return '{0:.3f} yrs'.format(B)
        elif KB <= B < MB:
          return '{0:.3f} kyrs'.format(B/KB)
        elif MB <= B < GB:
          return '{0:.3f} Myrs'.format(B/MB)
        elif GB <= B < TB:
          return '{0:.3f} Gyrs'.format(B/GB)
        elif TB <= B:
          return '{0:.3f} Tyrs'.format(B/TB)

    pc_to_m = const.pc.value
    H0 = H0*(1e3/(1e6*pc_to_m))
    s_to_yr = (1*u.second).to(u.yr).value
    s_to_Gyr = s_to_yr/1e9
    E = lambda z: omega_m0*(1+z)**3 + omega_k0*(1+z)**2 + omega_lam0 + omega_rad0*(1+z)**4

    def modified_z_integrand(z):
        Z = -s_to_Gyr/(H0*(1+z)*np.sqrt(E(z)))
        return Z

    t,dt = integrate.quad(modified_z_integrand,np.inf,z)
    t0,dt = integrate.quad(modified_z_integrand,np.inf,0)
    t_b = t0-t
    t = humanyears(t*1e9)
    t0 = humanyears(t0*1e9)
    t_b = humanyears(t_b*1e9)
    return t,t0,t_b#np.round(t,3),np.round(t0,3)

def distance(z, H0, omega_m0, omega_rad0, omega_lam0, omega_k0):

    def humandist(B):
        B = float(B)
        KB = float(1000)
        MB = float(KB ** 2)
        GB = float(KB ** 3)
        TB = float(KB ** 4)
        mB = float(1e-3)

        if B < KB:
          return '{0:.3f} pc'.format(B)
        elif KB <= B < MB:
          return '{0:.3f} kpc'.format(B/KB)
        elif MB <= B < GB:
          return '{0:.3f} Mpc'.format(B/MB)
        elif GB <= B < TB:
          return '{0:.3f} Gpc'.format(B/GB)
        elif TB <= B:
          return '{0:.3f} Tpc'.format(B/TB)
        elif mB <= B < KB:
          return '{0:.3f} mpc'.format(B/mB)

    def humanscale(B):
        B = float(B)
        KB = float(1000)
        mB = float(1e-3)

        if B < KB:
          return '{0:.3f} pc/arcsec'.format(B)
        elif mB <= B < KB:
          return '{0:.3f} mpc/arcsec'.format(B/mB)
        elif KB <= B:
          return '{0:.3f} kpc/arcsec'.format(B/KB)

    pc_to_m = const.pc.value
    m_to_pc = 1.0/pc_to_m
    m_to_Mpc = 1e-6*m_to_pc
    H0 = H0*(1e3/(1e6*pc_to_m))
    c = const.c.value
    E = lambda z: omega_m0*(1+z)**3 + omega_k0*(1+z)**2 + omega_lam0 + omega_rad0*(1+z)**4

    def integrand(z):
        Z = (c*m_to_Mpc)/(H0*np.sqrt(E(z)))
        return Z
    def horizon_integrand(z):
        Z = (c*m_to_Mpc)/(H0*(1+z)*np.sqrt(E(z)))
        return Z

    rc,drc = integrate.quad(integrand,0,z)
    dA,dL = rc/(1+z), (1+z)*rc
    dAs = (dA*1e3*np.pi)/(180.0*60.0*60.0)
    Hl,dHl = integrate.quad(horizon_integrand,z,np.inf)
    rc = humandist(rc*1e6)
    dA = humandist(dA*1e6)
    dL = humandist(dL*1e6)
    Hl = humandist(Hl*1e6)
    dAs = humanscale(dAs*1e3)
    return rc,dA,dL,dAs,Hl#np.round(rc,3),np.round(dA,3),np.round(dL,3),np.round(dAs,3),np.round(Hl,3)

def temperature(z):
    return 2.725*(1+z)
