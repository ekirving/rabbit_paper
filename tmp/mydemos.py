"""
    One-population demographic models.
"""

import numpy

from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

def snm(notused, ns, pts):
    """
        Standard neutral model.
        
        ns = (n1,)
        
        n1: Number of samples in resulting Spectrum
        pts: Number of grid points to use in integration.
        """
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

def growth(params, ns, pts):
    """
        Exponential growth beginning some time ago.
        
        params = (nu,T)
        ns = (n1,)
        
        nu: Ratio of contemporary to ancient population size
        T: Time in the past at which growth began (in units of 2*Na
        generations)
        n1: Number of samples in resulting Spectrum
        pts: Number of grid points to use in integration.
        """
    nu,T = params
    
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    
    nu_func = lambda t: numpy.exp(numpy.log(nu) * t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)
    
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

def bottlegrowth(params, ns, pts):
    """
        Instantanous size change followed by exponential growth.
        
        params = (nuB,nuF,T)
        ns = (n1,)
        
        nuB: Ratio of population size after instantanous change to ancient
        population size
        nuF: Ratio of contemporary to ancient population size
        T: Time in the past at which instantaneous change happened and growth began
        (in units of 2*Na generations)
        n1: Number of samples in resulting Spectrum
        pts: Number of grid points to use in integration.
        """
    nuB,nuF,T = params
    
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    
    nu_func = lambda t: nuB*numpy.exp(numpy.log(nuF/nuB) * t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)
    
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs


"""
    Two-population demographic models.
"""

import numpy

from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

def snm2(notused, ns, pts):
    """
        ns = (n1,n2)
        
        Standard neutral model, populations never diverge.
        """
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def IM(params, ns, pts):
    """
        ns = (n1,n2)
        params = (s,nu1,nu2,T,m12,m21)
        
        Isolation-with-migration model with exponential pop growth.
        
        s: Size of pop 1 after split. (Pop 2 has size 1-s.)
        nu1: Final size of pop 1.
        nu2: Final size of pop 2.
        T: Time in the past of split (in units of 2*Na generations)
        m12: Migration from pop 2 to pop 1 (2*Na*m12)
        m21: Migration from pop 1 to pop 2
        n1,n2: Sample sizes of resulting Spectrum
        pts: Number of grid points to use in integration.
        """
    s,nu1,nu2,T,m12,m21 = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def Iso(params, ns, pts):
    """
        This model is a special case of the IM model, refer to it for parameter definitions.  Only difference is that migration is cut off here.
    """
    
    s,nu1,nu2,T = params
    return IM((s, nu1, nu2, T, 0, 0), ns, pts)

