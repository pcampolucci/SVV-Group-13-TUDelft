"""
Title: Functions for aerodynamic distributed load discretization
"""

import numpy as np
from src.input.input import Input


def get_discrete_load(x, cont_load, step):
    """ Given a continous load function q(x), this will make an array of the load at different
    locations with a set interval. For the trapezoidal rule"""
    discrete_load = np.empty_like(np.arange(0, x+step, step))
    for i in np.arange(0, x+step, step):
        discrete_load[int(round(i/step))] = cont_load.get_q(-i)
    return discrete_load


def trapezoidal_rule(row, step):
    """ Just trapezoidal rule between set of points"""
    resultant = 0
    for i in range(len(row)-1):
        r_i = (row[i-1] + row[i])*step*0.5
        resultant += r_i
    return resultant


def magnitude_resultant(x, cont_load, step):
    """ Finds resultant force of distribution from 0 till x_end according to given span distr.
        First it takes points from that distr and then uses trapezoidal rule. """
    return trapezoidal_rule(get_discrete_load(x, cont_load, step), step)


def location_resultant(x, cont_load, step):
    """ Finds location of application resultant force. With formula:
            xbar  = integral(x*q(x))/integral(q(x)) """
    discrete_load = get_discrete_load(x, cont_load, step)
    resultant = magnitude_resultant(x, cont_load, step)
    discrete_load_x = discrete_load*np.arange(0, x+step, step)

    return trapezoidal_rule(discrete_load_x, step)/resultant


def moment_resultant(x, cont_load, step):
    """ Finds moment with respect to end point """
    resultant = magnitude_resultant(x, cont_load, step)
    location = location_resultant(x, cont_load, step)
    return resultant*(location-x)


def angle_distributed(x, cont_load, step):
    discrete_moment = np.zeros(np.linspace(0, x, int(round((x+step)/step, int(-np.log10(step))))).shape)
    for i in np.linspace(0, x, int(round((x+step)/step, int(-np.log10(step))))):
        if i > 0:
            discrete_moment[int(round(i/step))-1] = moment_resultant(i, cont_load, step)
    return trapezoidal_rule(discrete_moment, step)


def deflection_distributed(x, cont_load, step):
    discrete_angle = np.zeros(np.linspace(0, x, int(round((x+step)/step, int(-np.log10(step))))).shape)
    for i in np.linspace(0, x, int(round((x+step)/step, int(-np.log10(step))))):
        discrete_angle[int(round(i/step))-1] = angle_distributed(i, cont_load, step)
    return trapezoidal_rule(discrete_angle, step)


# ===============================================================================
DEBUG = False

if DEBUG:
    # inputs
    x_end = 1  # [m] end point, set calc will be done for distr from 0 till this point
    stepsize_aero = 0.001  # [m] set the distance between points in trapezoidal rule
    load = Input('A').aero_input()

    print(deflection_distributed(0, load, stepsize_aero))
