"""
Script will initialize all the loads so that the combination scripts do not need any multiple initialization
"""

from tqdm import tqdm
from src.loads.moment import Moment
from src.loads.shear import Shear
from src.loads.torque import Torque
from src.loads.discrete_load import PointLoads
from src.loads.distributed_load import *  # not a class, just sending and executing functions
import matplotlib.pyplot as plt


class Loads:

    def __init__(self, aircraft, steps, step_size):

        print("Initializing load class ...")
        self.a = aircraft
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

    def plot_forces(self):

        print("=" * 100)
        print(f"Plotting Forces Report for Aircraft: {self.a}")
        print("=" * 100, "\n")

        # plot input dictionary of all the forces obtained
        force_lst = self.point_loads
        string_lst = ["F_z1", "F_z2", "F_z3", "F_a", "F_y1", "F_y2", "F_y3"]

        for i in range(7):
            print(f"| {string_lst[i]}: {force_lst[i]}")
        print()

        # plot distribution of forces along span
        plt.figure(figsize=[17, 10])
        x_axis = self.x_location

        plt.suptitle(f"Load Distribution for Aircraft: {self.a}, (steps={self.steps})")

        plt.subplot(3, 2, 1)
        plt.title(f"Moment Z")
        plt.plot(x_axis, self.send_moment_z(), color='b')
        plt.grid()

        plt.subplot(3, 2, 2)
        plt.title(f"Moment Y")
        plt.plot(x_axis, self.send_moment_y(), color='b')
        plt.grid()

        plt.subplot(3, 2, 3)
        plt.title(f"Shear Z")
        plt.plot(x_axis, self.send_shear_z(), color='b')
        plt.grid()

        plt.subplot(3, 2, 4)
        plt.title(f"Shear Y")
        plt.plot(x_axis, self.send_shear_y(), color='b')
        plt.grid()

        plt.subplot(3, 2, 5)
        plt.title(f"Torque")
        plt.plot(x_axis, self.send_torque(), color='b')
        plt.grid()

        plt.show()


# =========================================================================================================
# DEBUGGING IN GENERAL

DEBUG = False

if DEBUG:

    load_init = Loads('A', 1000, 0.001)  # by doing this we initialize every single load

    # get the time to send all the completed loads to destination
    load_init.send_moment_y()
    load_init.send_moment_z()
    load_init.send_shear_y()
    load_init.send_shear_z()
    load_init.send_torque()

    # plot results
    load_init.plot_forces()
