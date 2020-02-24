import numpy as np
from src.input.input import Input

# inputs
x_end = 1  # [m] end point, set calc will be done for distr from 0 till this point
stepsize_aero = 0.1  # [m] set the distance between points in trapezoidal rule
load = Input('B').aero_input()



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
    # print('mom', 'x:',x, 'res', resultant)
    location = location_resultant(x, cont_load, step)
    return resultant*(location-x)

# Test with boeing data ( which is a cte load)
#Qdisc = get_discrete_load(x_end, load, stepsize_aero)
#Q = magnitude_resultant(x_end, load, stepsize_aero)
#Qx = location_resultant(x_end, load, stepsize_aero)
#QM = moment_resultant(x_end, load, stepsize_aero)


# --------------------------------------------------------------------------------------------

def angle_distributed(x, cont_load, step):
    discrete_moment = np.zeros((int(round(x+step)/step), round(1/step)))
    for i in np.arange(0+step, x+step, step):
        discrete_moment[int(round(i/step))] = moment_resultant(i, cont_load, step)
    return trapezoidal_rule(discrete_moment, step)


def deflection_distributed(x, cont_load, step):
    discrete_angle = np.empty((int(round(x+step)/step), round(1/step)))
    for i in np.arange(0+step, x+step, step):
        discrete_angle[int(round(i/step))] = angle_distributed(i, cont_load, step)
    return trapezoidal_rule(discrete_angle, step)


print(deflection_distributed(0.04, load, 0.01))

"""
def moment_resultant1(x, cont_load, step):
     Finds moment with respect to end point 
    discrete_force = np.zeros_like(np.arange(0, x, step))
    for i in np.arange(0, x, step):
        discrete_force[int(i/step)] = magnitude_resultant(i, cont_load, step)
        print(discrete_force[int(round(i/step))], '?')
    return trapezoidal_rule(discrete_force, step)
"""
