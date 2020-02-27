"""
Title: Functions for aerodynamic distributed load discretization
"""

import numpy as np
from src.input.input import Input


# ===================  Global inputs to generate arrays  =========================
la = 2.661
stepsize = 0.1  # [m] set the distance between points in trapezoidal rule
load = Input('A').aero_input()

# ======================  8 Functions  ======================================
""" The get_discrete_xxxx functions make discrete functions for the respective xxxxx feature.
  The xxxxx_resultants use these discrete functions/arrays then to derrive the approximated value at an exact input 
  location. 
  NOTE: The input of the get_discrete_xxx functions should always be the total length of the aileron, i.e. length 
  aileron (la) """


def trapezoidal_rule(row, step):
    """ Just trapezoidal rule between set of points"""
    resultant = 0
    for i in range(len(row)-1):
        r_i = (row[i-1] + row[i])*step*0.5
        resultant += r_i
    return resultant

# --------------- Discrete Load -----------------------


def get_discrete_load(x, cont_load, step):
    """ Given a continous load function q(x), this will make an array of the load at different
    locations with a set interval. For the trapezoidal rule"""
    discrete_load = np.empty_like(np.arange(0, x+step, step))
    for i in np.arange(0, x+step, step):
        discrete_load[int(round(i/step))] = cont_load.get_q(-i)
    return discrete_load


# --------------- Discrete resultant -----------------------


def get_discrete_resultant(la, discrete_load, step):
    """ Make discrete resultant function """
    discrete_resultant = np.zeros_like(np.arange(0, la+step, step))
    for i in np.arange(step, la+step, step):
        discrete_resultant[int(round(i/step))] = trapezoidal_rule(discrete_load[0:int(i/step)+1], step)
    return discrete_resultant


def magnitude_resultant(x, discrete_resultant, step):
    """ Finds resultant force of distribution from 0 till x_end according to given span distr.
        First it takes points from that distr and then uses trapezoidal rule. """
    if int((x+step)/step) >= len(discrete_resultant):
        return discrete_resultant[int(x/step)]
    return 0.5*(discrete_resultant[int((x+step)/step)]+discrete_resultant[int(x/step)])


# --------------- Discrete locations -----------------------


def get_discrete_location_resultant(la, discrete_resultant, discrete_load, step):
    """ Finds location of application resultant force. With formula:
            xbar  = integral(x*q(x))/integral(q(x)) """
    discrete_location = np.zeros_like(np.arange(0, la+step, step))
    discrete_resultant_x = discrete_load*np.arange(0, la+step, step)
    for i in np.arange(step, la, step):
        discrete_location[int(round(i/step))] = trapezoidal_rule(discrete_resultant_x[1:int(round(i/step))+2], step) / magnitude_resultant(i, discrete_resultant, step)
    return discrete_location


def location_resultant(x, discrete_location, step):
    """ Finds resultant force of distribution from 0 till x_end according to given span distr.
        First it takes points from that distr and then uses trapezoidal rule. """
    if int((x+step)/step) >= len(discrete_location):
        return discrete_location[int(x/step)]
    return 0.5*(discrete_location[int((x+step)/step)]+discrete_location[int(x/step)])


# --------------- Discrete moments -----------------------


def get_discrete_moment(discrete_resultant, discrete_location):
    """ Finds moment with respect to end point """
    return discrete_resultant*discrete_location


def moment_resultant(x, discrete_moment, step):
    """ Finds resultant force of distribution from 0 till x_end according to given span distr.
        First it takes points from that distr and then uses trapezoidal rule. """
    if int((x+step)/step) >= len(discrete_moment):
        return discrete_moment[int(x/step)]
    return 0.5*(discrete_moment[int((x+step)/step)]+discrete_moment[int(x/step)])


# --------------- Discrete angles -----------------------


def get_discrete_angle(la, discrete_moment, step):
    """ Make discrete resultant function """
    discrete_angle = np.zeros_like(np.arange(0, la+step, step))
    for i in np.arange(step, la+step, step):
        discrete_angle[int(round(i/step))] = trapezoidal_rule(discrete_moment[0:int(i/step)+1], step)
    return discrete_angle


def angle_resultant(x, discrete_angle, step):
    """ Finds resultant force of distribution from 0 till x_end according to given span distr.
        First it takes points from that distr and then uses trapezoidal rule. """
    if int((x+step)/step) >= len(discrete_angle):
        return discrete_angle[int(x/step)]
    return 0.5*(discrete_angle[int((x+step)/step)]+discrete_angle[int(x/step)])


# --------------- Discrete deflections -----------------------


def get_discrete_deflection(la, discrete_angle, step):
    """ Make discrete deflection function """
    discrete_deflection = np.zeros_like(np.arange(0, la+step, step))
    for i in np.arange(step, la+step, step):
        discrete_deflection[int(round(i/step))] = trapezoidal_rule(discrete_angle[0:int(i/step)+1], step)
    return discrete_deflection


def deflection_resultant(x, discrete_deflection, step):
    """ Finds resultant force of distribution from 0 till x_end according to given span distr.
        First it takes points from that distr and then uses trapezoidal rule. """
    if int((x+step)/step) > len(discrete_deflection):
        return discrete_deflection[int(x/step)]
    return 0.5*(discrete_deflection[int((x+step)/step)]+discrete_deflection[int(x/step)])


# ========================= Arrays ===========================
"""  The discrete functions for the respective features.  """
discrete_loads = get_discrete_load(la, load, stepsize)
discrete_resultants = get_discrete_resultant(la, discrete_loads, stepsize)
discrete_locations = get_discrete_location_resultant(la, discrete_resultants, discrete_loads, stepsize)
discrete_moments = get_discrete_moment(discrete_resultants, discrete_locations)
discrete_angles = get_discrete_angle(la, discrete_moments, stepsize)
discrete_deflections = get_discrete_deflection(la, discrete_angles, stepsize)


# ===============================================================================
""" Checks for constant load -55.7 N/m (The load case of the B737) """
DEBUG = False

if DEBUG:
    # inputs
    x_end = 2.661  # [m] end point, set calc will be done for distr from 0 till this point
    stepsize = 0.001  # [m] set the distance between points in trapezoidal rule
    load = Input('B').aero_input()

    # test
    res1 = magnitude_resultant(1, discrete_resultants, stepsize)
    print('Resultant should be -55.7 = ', res1)
    loc1 = location_resultant(1, discrete_locations, stepsize)
    print('location should be 0.5 = ', loc1)
    mom1 = moment_resultant(1, discrete_moments, stepsize)
    print('moment should be ', -55.7/2, ' = ', mom1)
    ang1 = angle_resultant(1, discrete_angles, stepsize)
    print('Angle should be ', -55.7/2/3, ' = ', ang1)
    def1 = deflection_resultant(1, discrete_deflections, stepsize)
    print('Deflection should be ', -55.7/2/3/4, ' = ', def1)