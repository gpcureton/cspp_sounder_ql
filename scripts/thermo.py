#!/usr/bin/env python
# encoding: utf-8
"""
thermo.py

Various thermodynamic relationship for ice and water.

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2009-04-04.
Copyright (c) 2009-2015 University of Wisconsin Regents. All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from scipy import log10
from numpy import exp, log
from numpy import min, minimum, where

# -----------------------------------------------------------------------
# Here we go. A set of functions that I use from time to time to calculate
# the basic stuff that I'm sick of doing over and over! I'm going to
# endeavour to include references and global constants to make it all nice
# and legible.
# -----------------------------------------------------------------------

Rs_da=287.05          # Specific gas const for dry air, J kg^{-1} K^{-1}
Rs_v=461.51           # Specific gas const for water vapour, J kg^{-1} K^{-1}
Cp_da=1004.6          # Specific heat at constant pressure for dry air
Cv_da=719.            # Specific heat at constant volume for dry air
Cp_v=1870.            # Specific heat at constant pressure for water vapour
Cv_v=1410.            # Specific heat at constant volume for water vapour
Cp_lw=4218	      # Specific heat at constant pressure for liquid water
Epsilon=0.622         # Epsilon=Rs_da/Rs_v; The ratio of the gas constants
degCtoK=273.15        # Temperature offset between K and C (deg C)
rho_w=1000.           # Liquid Water density kg m^{-3}
grav=9.81             # Gravity, m s^{-2}
Lv=2.5e6              # Latent Heat of vaporisation
boltzmann=5.67e-8     # Stefan-Boltzmann constant
mv=18.0153            # Mean molar mass of water vapor(g/mol)


def rh_to_mr( rh, p, t) :
	'''
	Returns mixing ratio, in g/kg, given relative humidity in %,
	pressure in hPa and temperature in K.
	'''
	return rh * 0.01 * satmix(p, t)

def rh_to_mr_wat( rh, p, t) :
	'''
	Returns mixing ratio over water, in g/kg, given relative humidity in %,
	pressure in hPa and temperature in K.
	'''
	return rh * 0.01 * satmixwat(p, t)

def rh_to_mr_ice( rh, p, t) :
	'''
	Returns mixing ratio over ice, in g/kg, given relative humidity in %,
	pressure in hPa and temperature in K.
	'''
	return rh * 0.01 * satmixice(p, t)

def mr_to_rh( mr,  p,  t) :
	'''
	Returns relative humidity in %, given the mixing ratio in g/kg,
	pressure in hPa and temperature in K.
	'''
	return mr * 100. / satmix(p, t)

def mr_to_rh_wat( mr,  p,  t) :
	'''
	Returns relative humidity in %, given the mixing ratio over water in g/kg,
	pressure in hPa and temperature in K.
	'''
	return mr * 100. / satmixwat(p, t)

def mr_to_rh_ice( mr,  p,  t) :
	'''
	Returns relative humidity in %, given the mixing ratio over ice in g/kg,
	pressure in hPa and temperature in K.
	'''
	return mr * 100. / satmixice(p, t)

def satmix( p, t) :
	'''
	Returns saturation mixing ratio in g/kg, given pressure in hPa and
	temperature in K.
	'''
	if (t > 253.) :
		return satmixwat(p, t)
	else :
		return satmixice(p, t)

def satmixwat( p,  t) :
	'''
	Returns saturation mixing ratio over water, in g/kg, given pressure in hPa and
	temperature in K.
	'''
	es = svpwat(t)
	return (622. * es)/p

def satmixice( p, t) :
	'''
	Returns saturation mixing ratio over ice, in g/kg, given pressure in hPa and
	temperature in K.
	'''
	es = svpice(t);
	return (622. * es) / p;


def svpwat(t) :
    '''
    Returns saturation vapor pressure over water, in hPa, given temperature in K.
    '''

    a0 =  0.999996876e0
    a1 = -0.9082695004e-2
    a2 =  0.7873616869e-4
    a3 = -0.6111795727e-6
    a4 =  0.4388418740e-8
    a5 = -0.2988388486e-10
    a6 =  0.2187442495e-12
    a7 = -0.1789232111e-14
    a8 =  0.1111201803e-16
    a9 = -0.3099457145e-19
    b = 0.61078e+1
    t -= 273.16
    return (b / ((a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*a9)))))))))**8.))

def svpice( t) :
    '''
    Returns saturation vapor pressure over ice, in hPa, given temperature in K.
    The Goff-Gratch equation (Smithsonian Met. Tables,  5th ed., pp. 350, 1984)
    '''
    a = 273.16 / t
    exponent = -9.09718 * (a - 1.) - 3.56654 * log10(a) + 0.876793 * (1. - 1./a) + log10(6.1071)

    return 10.0**exponent


def dewpoint_magnus(T,RH) :
    '''
    Returns the dewpoint temperature in dec C, using the Magnus formula.

    Inputs:

    T : temperature in Kelvin
    RH: relative humidity in %

    '''

    #a = 6.112 # millibars
    b = 17.67
    c = 243.5 # deg C

    #T = T - degCtoK
    gamma = log(RH/100.) + (b*T)/(c+T)

    T_d = (c * gamma) / (b - gamma)

    return T_d


def dewpoint_AB(T_K,RH) :
    '''
    Returns the dewpoint temperature in dec C, using the  Arden Buck formula.

    Inputs:

    T_K : temperature in Kelvin
    RH: relative humidity in %

    '''

    # Convert T to deg Celcius
    T = T_K - degCtoK

    # Contants from 1980 paper by David Bolton in the Monthly Weather Review
    #a = 6.112 # millibars
    #b = 17.67
    #c = 243.5 # deg C

    # Contants from Buck(1996)
    #a = 6.1121 # millibars
    #b = 18.678
    #c = 257.14 # deg C

    # Contants from Sonntag (1990)
    #a = 6.1121 # millibars
    #b = 17.62
    #c = 243.12 # deg C

    if ((T > -40.) and (T <= 0.)) :
        #print("{0:5.2f} is > -40. and {0:5.2f} is  <= 0."). \
                #format(T)
        a = 6.1121 # millibars
        b = 17.368
        c = 238.88 # deg C
    elif ((T > 0.) and (T <= 50.)) :
        #print("{0:5.2f} is > 0. and {0:5.2f} is <= 50."). \
                #format(T)
        a = 6.1121 # millibars
        b = 17.966
        c = 247.15 # deg C
    else :
        #print("Temperature {} degC ({} K) out of range, using defaults..."). \
                #format(T,T_K)
        a = 6.1121 # millibars
        b = 18.678
        c = 257.14 # deg C

    # Extra parameter
    d = 234.5  # deg C

    exponent = (b - (T/d)) * (T/(c + T))

    gamma = log((RH/100.) * exp(exponent))

    T_d = (c * gamma) / (b - gamma)

    return T_d


    #########################################

def Theta(tempk,pres,pref=100000.):
    """Potential Temperature

    INPUTS:
    tempk (K)
    pres (Pa)
    pref: Reference pressure (default 100000 Pa)

    OUTPUTS: Theta (K)

    Source: Wikipedia
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """

    try:
        minpres=min(pres)
    except TypeError:
        minpres=pres

    if minpres<2000:
        print("WARNING: P<2000 Pa; did you input a value in hPa?")

    return tempk*(pref/pres)**(Rs_da/Cp_da)


def TempK(theta,pres,pref=100000.):
    """Inverts Theta function."""

    try:
        minpres=min(pres)
    except TypeError:
        minpres=pres

    if minpres<2000:
        print("WARNING: P<2000 Pa; did you input a value in hPa?")

    return theta*(pres/pref)**(Rs_da/Cp_da)


def ThetaE():
    """Equivalent potential temperature"""
    raise NotImplementedError


def ThetaV(tempk,pres,e):
    """Virtual Potential Temperature

    INPUTS
    tempk (K)
    pres (Pa)
    e: Water vapour pressure (Pa) (Optional)
    """

    mixr=MixRatio(e,pres)
    theta=Theta(tempk,pres)

    return theta*(1+mixr/Epsilon)/(1+mixr)


def GammaW(tempk,pres,e=None):
    """Function to calculate the moist adiabatic lapse rate (deg C/Pa) based
    on the temperature, pressure, and rh of the environment.

    INPUTS:
    tempk (K)
    pres (Pa)
    RH (%)

    RETURNS:
    GammaW: The moist adiabatic lapse rate (Dec C/Pa)
    """

    tempc=tempk-degCtoK
    es=SatVap(tempc)
    ws=MixRatio(es,pres)

    if e is None:
        # assume saturated
        e=es

    w=MixRatio(e,pres)

    tempv=VirtualTempFromMixR(tempk,w)
    latent=Latentc(tempc)

    A=1.0+latent*ws/(Rs_da*tempk)
    B=1.0+Epsilon*latent*latent*ws/(Cp_da*Rs_da*tempk*tempk)
    Rho=pres/(Rs_da*tempv)
    Gamma=(A/B)/(Cp_da*Rho)
    return Gamma


def DensityHumid(tempk,pres,e):
    """Density of moist air.
    This is a bit more explicit and less confusing than the method below.

    INPUTS:
    tempk: Temperature (K)
    pres: static pressure (Pa)
    mixr: mixing ratio (kg/kg)

    OUTPUTS:
    rho_air (kg/m^3)

    SOURCE: http://en.wikipedia.org/wiki/Density_of_air
    """

    pres_da=pres-e
    rho_da=pres_da/(Rs_da*tempk)
    rho_wv=e/(Rs_v*tempk)

    return rho_da+rho_wv


def Density(tempk,pres,mixr):
    """Density of moist air

    INPUTS:
    tempk: Temperature (K)
    pres: static pressure (Pa)
    mixr: mixing ratio (kg/kg)

    OUTPUTS:
    rho_air (kg/m^3)
    """

    virtualT=VirtualTempFromMixR(tempk,mixr)
    return pres/(Rs_da*virtualT)


def VirtualTemp(tempk,pres,e):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    e: vapour pressure (Pa)
    p: static pressure (Pa)

    OUTPUTS:
    tempv: Virtual temperature (K)

    SOURCE: hmmmm (Wikipedia)."""

    tempvk=tempk/(1-(e/pres)*(1-Epsilon))
    return tempvk


def VirtualTempFromMixR(tempk,mixr):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    mixr: Mixing Ratio (kg/kg)

    OUTPUTS:
    tempv: Virtual temperature (K)

    SOURCE: hmmmm (Wikipedia). This is an approximation
    based on a m
    """

    return tempk*(1.0+0.6*mixr)


def Latentc(tempc):
    """Latent heat of condensation (vapourisation)

    INPUTS:
    tempc (C)

    OUTPUTS:
    L_w (J/kg)

    SOURCE:
    http://en.wikipedia.org/wiki/Latent_heat#Latent_heat_for_condensation_of_water
    """

    return 1000*(2500.8 - 2.36*tempc + 0.0016*tempc**2 - 0.00006*tempc**3)


def SatVap(tempc,phase="liquid"):
    """Calculate saturation vapour pressure over liquid water and/or ice.

    INPUTS:
    tempc: (C)
    phase: ['liquid'],'ice'. If 'liquid', do simple dew point. If 'ice',
    return saturation vapour pressure as follows:

    Tc>=0: es = es_liquid
    Tc <0: es = es_ice


    RETURNS: e_sat  (Pa)

    SOURCE: http://cires.colorado.edu/~voemel/vp.html (#2:
    CIMO guide (WMO 2008), modified to return values in Pa)

    This formulation is chosen because of its appealing simplicity,
    but it performs very well with respect to the reference forms
    at temperatures above -40 C. At some point I'll implement Goff-Gratch
    (from the same resource).
    """

    over_liquid=6.112*exp(17.67*tempc/(tempc+243.12))*100.
    over_ice=6.112*exp(22.46*tempc/(tempc+272.62))*100.
    # return where(tempc<0,over_ice,over_liquid)

    if phase=="liquid":
        # return 6.112*exp(17.67*tempc/(tempc+243.12))*100.
        return over_liquid
    elif phase=="ice":
        # return 6.112*exp(22.46*tempc/(tempc+272.62))*100.
        return where(tempc<0,over_ice,over_liquid)
    else:
        raise NotImplementedError


def MixRatio(e,p):
    """Mixing ratio of water vapour
    INPUTS
    e (Pa) Water vapor pressure
    p (Pa) Ambient pressure

    RETURNS
    qv (kg kg^-1) Water vapor mixing ratio`
    """

    return Epsilon*e/(p-e)


def MixR2VaporPress(qv,p):
    """Return Vapor Pressure given Mixing Ratio and Pressure
    INPUTS
    qv (kg kg^-1) Water vapor mixing ratio`
    p (Pa) Ambient pressure

    RETURNS
    e (Pa) Water vapor pressure
    """

    return qv*p/(Epsilon+qv)


def VaporPressure(dwpt):
    """Water vapor pressure
    INPUTS
    dwpt (C) Dew Point Temperature (for SATURATION vapor
	     pressure use tempc)

    RETURNS
    e (Pa) Water Vapor Pressure

    SOURCE:
    Bolton, Monthly Weather Review, 1980, p 1047, eq. (10)
    """

    return 611.2*exp(17.67*dwpt/(243.5+dwpt))


def DewPoint(e):
    """ Use Bolton's (1980, MWR, p1047) formulae to find tdew.
    INPUTS:
    e (Pa) Water Vapor Pressure
    OUTPUTS:
    Td (deg C)
      """

    ln_ratio=log(e/611.2)
    Td=((17.67-ln_ratio)*degCtoK+243.5*ln_ratio)/(17.67-ln_ratio)
    return Td-degCtoK


#
# Code from Nadia Smith and E. Weisz, which originated from
# Hal Woolf.
#

def svpwat(temp):

    """
    -----------------------------------------------------
     Calculate saturation vapor pressure over water given
     a certain temperature
    -------------------------------------------------------
     input:  temperature in Kelvin or degrees Celsius
     output: saturation vapor pressure at given temp over water in  mb or hPa
    -----------------------------------------------------
     NOTE: Due to the prevalence of the super-cooling phenomenon
     in clouds and upper atmosphere relative humidity is calculated
     for conditions over water irrespective of the temperature.
     This is according to WMO standards.
     Consequently, the vp threshold was lowered from 0.636d-1 to 0.1403d-4

     REFERENCE: Hardy, Bob, "ITS-90 Formulations of ..." Papers and
     Abstracts from teh Third International Symposium on Humidity &
     Moisture, London, England, April 1998, vol. 1, 214-222.
    -------------------------------------------------------

    python version of Hal Woolf's fortran code 'svpwat.f'

    """

    tref=273.15

    a0 =  0.999996876e0
    a1 = -0.9082695004e-2
    a2 =  0.7873616869e-4
    a3 = -0.6111795727e-6
    a4 =  0.4388418740e-8
    a5 = -0.2988388486e-10
    a6 =  0.2187442495e-12
    a7 = -0.1789232111e-14
    a8 =  0.1111201803e-16
    a9 = -0.3099457145e-19
    b  =  0.61078e+1

    t = temp

    if (t > 100.):
        t = t - tref

    s = a0 + t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*a9))))))))
    s = b / (s**8)

    return s


def tsvwat(press):
    """
    -----------------------------------------------------
     Calculates temperature [K] at given saturation vapor pressure [hPa]
     zero is returned if the given vapor pressure is out of range
    -------------------------------------------------------
     input: saturation vapor pressure [hPa]
     output: temperature at which saturation vapor pressure occur over water (not ice)
    -----------------------------------------------------
     NOTE: Due to the prevalence of the super-cooling phenomenon
     in clouds the upper atmosphere relative humidity is calculated
     for conditions over water irrespective of the temperature.
     This is according to WMO standards.
     Consequently, the vp threshold was lowered from 0.636d-1 to 0.1403d-4

     REFERENCE: Hardy, Bob, "ITS-90 Formulations of ..." Papers and
     Abstracts from teh Third International Symposium on Humidity &
     Moisture, London, England, April 1998, vol. 1, 214-222.
    -------------------------------------------------------

    matlab version of Hal Woolf's fortran code 'tsvwat.f'

    ------------------------------------------------------------
    """

    a0 = -0.225896152438e+2
    a1 =  0.261012286592e+2
    a2 =  0.30206720594e+1
    a3 =  0.370219024579e+0
    a4 =  0.72838702401e-1

    vp=press

    # original:  if (vp > 0.636d-1 & vp < 0.1233972d+3)
    if (vp > 0.1403e-4 and vp < 0.1233972e+3):
        vp = log10(vp)
        t = a0 + vp*(a1+vp*(a2+vp*(a3+vp*a4)))
        ts = t + 273.16
    else :
        ts=0.

    return ts


def sat_mr_water(p,t):
    """
    * Get saturation mixing ratio over water from press, temp

    input: pressure in mbar,hPa
    input: temperature in Kelvin, deg C

    returns: saturation water vapour mixing ratio in g/kg
    """

    es = svpwat(t)

    ws = 621.946*es/p  # saturation MR [g/kg]

    return ws


def dewhum(p,t,w):
    """
    * Computes relative humidity [%] from wv
    * Checks for oversaturation and return corrected wv values
    * Calculates dewpoint temperature

    NOTE: Due to the prevalence of the super-cooling phenomenon
    in clouds and upper atmosphere relative humidity is calculated
    for conditions over water irrespective of the temperature.
    This is according to WMO standards.

    REFERENCE: Hardy, Bob, "ITS-90 Formulations of ..." Papers and
    Abstracts from teh Third International Symposium on Humidity &
    Moisture, London, England, April 1998, vol. 1, 214-222.

    * Get dewpoint and relative humidity from press, temp, wvmr
    .... version of 22.04.98

    python version of Hal Woolf's fortran code 'dewhum.f'

    -----------------------------------------------------------

    input: pressure in mbar,hPa
    input: temperature in Kelvin, deg C
    input: water vapor mixing ratio in g/kg

    returns: dewpoint temperature in Kelvin
    returns: relative humidity
    returns: water vapour mixing ratio in g/kg
    """

    tt = t

    es = svpwat(tt)

    ws = 621.946*es/p  # saturation MR [g/kg]

    w = minimum(w,ws) # if input w is over-saturated set to saturation value

    hh = w/ws

    h = hh * 100. # RH

    ee = hh * es # partial vapor pressure [mbar], = (w(i)/622)*p(i)

    d = tsvwat(ee)

    return d, h, w

