import numpy as np

def vol(R):
    return 4./3.*np.pi*R**3.

def CDFR(r,Re,alpha,beta):
    '''
    Ricciardelli+ 2014 CDF
    '''
    return (r/Re)**alpha*np.exp((r/Re)**beta -1.)

def dprofile_Hamaus(r, dc,rs,Rv,rhom):
    '''
    Compute the density profile following Hamaus+ 20214 fit
    Rv is the void radius
    dc, rs are the fit's free parameters at Rv
    '''
    a  = -2.*rs/Rv + 4.
    if rs/Rv <= 0.91:
        b = 17.5*rs/Rv - 6.5
    else:
        b = -9.8*rs/Rv + 18.4

    rho = rhom * (dc * ((1. - (r/rs)**a)/(1. + (r/Rv)**b)) + 1.) #density

    return rho
