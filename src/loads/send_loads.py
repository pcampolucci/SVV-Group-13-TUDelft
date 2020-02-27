"""
Script will initialize all the loads so that the combination scripts do not need any multiple initialization
"""

from tqdm import tqdm
import numpy as np
from src.input.input import Input
from src.loads.moment import Moment
from src.loads.shear import Shear
from src.loads.torque import Torque
from src.loads.discrete_load import PointLoads
from src.loads.distributed_load import *  # not a class, just sending and executing functions


class Loads:

    def __init__(self, aircraft, steps, step_size):

        print("Initializing load class ...")
        self.steps = steps
        self.step_size = step_size

        self.input_dict = Input(aircraft).get_discrete()
        self.aero_load = Input(aircraft).aero_input()  # functionality get_q(x) for distributed aero load

        self.la = self.input_dict['la'][aircraft]
        print(" ... ")
        # information in discrete load
        self.discrete_input = PointLoads(aircraft).get_discrete_input()
        self.geometry_input = PointLoads(aircraft).get_geometry()
        self.point_loads = PointLoads(aircraft).get_discrete_loads()
        print(" ... ")
        # distributed loads from aerodynamic loads
        self.discrete_loads = get_discrete_load(self.la, self.aero_load, self.step_size)
        self.discrete_resultants = get_discrete_resultant(self.la, self.discrete_loads, self.step_size)
        self.discrete_locations = get_discrete_location_resultant(self.la, self.discrete_resultants, self.discrete_loads, self.step_size)
        self.discrete_moments = get_discrete_moment(self.discrete_resultants, self.discrete_locations)
        self.discrete_angles = get_discrete_angle(self.la, self.discrete_moments, self.step_size)
        self.discrete_deflections = get_discrete_deflection(self.la, self.discrete_angles, self.step_size)
        print(" ... ")
        self.x_location = np.linspace(1e-10, self.input_dict["la"][aircraft], self.steps)  # get loads like steps
        print(" ... Done\n")

    def send_moment_y(self):
        moment_lst = [float(Moment(self.discrete_input, self.point_loads,
                                   self.discrete_moments, self.step_size).M_y(i))
                      for i in tqdm(self.x_location, desc="Moment Y")]

        return moment_lst

    def send_moment_z(self):
        moment_lst = [float(Moment(self.discrete_input, self.point_loads,
                                   self.discrete_moments, self.step_size).M_z(i))
                      for i in tqdm(self.x_location, desc="Moment Z")]

        return moment_lst

    def send_shear_y(self):
        shear_lst = [float(Shear(self.discrete_input, self.point_loads,
                                 self.discrete_resultants, self.step_size, self.aero_load).V_y(i))
                      for i in tqdm(self.x_location, desc="Shear Y")]

        return shear_lst

    def send_shear_z(self):
        shear_lst = [float(Shear(self.discrete_input, self.point_loads,
                                 self.discrete_resultants, self.step_size, self.aero_load).V_z(i))
                     for i in tqdm(self.x_location, desc="Shear Z")]

        return shear_lst

    def send_torque(self):
        torque_lst = [float(Torque(self.discrete_input, self.point_loads,
                                   self.geometry_input).T(i))
                     for i in tqdm(self.x_location, desc="Torque")]

        return torque_lst


# =========================================================================================================
# DEBUGGING IN GENERAL

DEBUG = False

if DEBUG:

    load_init = Loads('A', 1000, 0.1)  # by doing this we initialize every single load

    # get the time to send all the completed loads to destination
    load_init.send_moment_y()
    load_init.send_moment_z()
    load_init.send_shear_y()
    load_init.send_shear_z()
    load_init.send_torque()
