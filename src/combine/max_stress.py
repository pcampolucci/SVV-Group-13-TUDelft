"""
Title: Max Shear Stress tool
"""

from src.input.cross_section.cross_section import CrossSection
from src.input.general.discrete_input import input_dict
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.loads.torque import Torque
from src.loads.moment import Moment
from tqdm import tqdm


class ShearStress:

    def __init__(self, aircraft, steps):
        self.input = input_dict
        self.aircraft = aircraft
        self.q_1 = CrossSection(input_dict, self.aircraft).twist_of_aileron([1], self.input['G'][self.aircraft])[0]
        self.q_2 = CrossSection(input_dict, self.aircraft).twist_of_aileron([1], self.input['G'][self.aircraft])[1]
        self.q_lst = CrossSection(input_dict, self.aircraft).get_shear_center()[1]
        self.x_location = np.linspace(0.0001, self.input["la"][self.aircraft], steps)
        self.n_points = steps
        self.T = [660.44862356,  359.64170576,  278.01919764, - 252.7818348,  3979.4813368, - 186.49317495,  100.51646629, - 276.95420611, - 314.78817203, - 839.02635342]#[float(Torque(self.aircraft, steps).T(i)) for i in tqdm(self.x_location, desc="Torque")]
        self.My = [1279.48835396,  6219.7912426,  17825.58106292, 29427.9314608, 48893.12680388, 54825.33865843, 39017.58694589, 24526.26153543, 10348.0636794,   1101.02127179]# [float(Moment(self.aircraft, steps).M_y(i)) for i in tqdm(self.x_location, desc="Moment Y")]
        self.Mz = [390.08069176, - 3988.5854398, - 11680.01329549, - 18662.70792584, - 21746.73532002, - 13060.85756981, - 10813.19731313, - 7099.87885604, - 2946.2971255,    743.22999279]# [float(Moment(self.aircraft, steps).M_z(i)) for i in tqdm(self.x_location, desc="Moment Z")]


    def shear_stress_due_to_shear(self):  # compute shear stress in skins and spar by dividing by thickness, note spar shear flow not included in shear_flow_lst
        # Shear stress distribution per section due to Sy and Sz computed at the middle of each section

        qb1, qb2, qb3, qb4, qb5, qb6 = self.q_lst
        t_sk = self.input['tsk'][self.aircraft]
        t_sp = self.input['tsp'][self.aircraft]

        tau_1 = [qb1[i] / t_sk for i in range(len(qb1))]
        tau_2 = [qb2[i] / t_sp for i in range(len(qb2))]
        tau_3 = [qb3[i] / t_sk for i in range(len(qb3))]
        tau_4 = [qb4[i] / t_sk for i in range(len(qb4))]
        tau_5 = [qb5[i] / t_sp for i in range(len(qb5))]
        tau_6 = [qb6[i] / t_sk for i in range(len(qb6))]
        return tau_1, tau_2, tau_3, tau_4, tau_5, tau_6

    def twist_of_aileron(self):

        """Solves a set of 3 equations for unit torque applied, output q1, q2 and twsit_rate_times_G
        For T_lst use the value of the torque at each location, for examply for adding a loop or using input out of a list of torques
        For verification, one can change the perimiter of the circle and triangle by chaning the geometry and calculate the q1,q2 and twist rate and see wheter it makes sense or not
        After computing twist rate, take distance from hinge line to shear center to compute the deflection of the hinge line"""

        # input
        T_lst = self.T
        t_sk = self.input['tsk'][self.aircraft]
        t_sp = self.input['tsp'][self.aircraft]
        l_a = self.input['la'][self.aircraft]
        G = self.input['G'][self.aircraft]
        A1, A2 = CrossSection(self.input, self.aircraft).cell_area()
        h = self.input['h'][self.aircraft]/2
        per_triangle = CrossSection(self.input, self.aircraft).stiffener_spacing()[3]
        per_semicircle = CrossSection(self.input, self.aircraft).stiffener_spacing()[2]

        # function
        k = 1 / (2 * A1) * (per_semicircle / t_sk + 2 * h / t_sp)
        l = 1 / (2 * A1) * -2 * h / t_sp
        m = 1 / (2 * A2) * -2 * h / t_sp
        n = 1 / (2 * A2) * (per_triangle / t_sk + 2 * h / t_sp)
        B = np.array([[2 * A1, 2 * A2, 0], [k, l, -1], [m, n, -1]])
        twist_rate_lst = []
        q1_lst = []
        q2_lst = []
        x_theta_0 = l_a / 2  # due to assumption around x, x_sc in middle [m]
        theta_0 = 0  # this is a boundary condition [rad]

        for i in tqdm(range(len(T_lst)), desc="twist of aileron"):
            w = np.array([T_lst[i], 0, 0])
            solution = np.linalg.solve(B, w)
            q1 = solution[0]
            q2 = solution[1]
            q1_lst.append(q1)
            q2_lst.append(q2)
            twist_rate_times_G = solution[2]
            twist_rate = twist_rate_times_G / G
            twist_rate_lst.append(twist_rate)

        J = T_lst[-1] / twist_rate_times_G  # calculate the J for a combination of torque and twist rate
        dx = l_a / (len(T_lst) - 1)  # step in x direction between the points where the torque is computed and thus where twist_rate is known
        n_steps = math.floor(x_theta_0 / dx)  # number of full steps untill location of boundary condition reached, returns an integer
        twist_before_bc = sum([twist_rate_lst[j] for j in range(n_steps)]) * dx + theta_0  # twist of first section
        twist_lst = [twist_before_bc]
        twist_after_bc = 0

        for i in range(1, len(T_lst)):
            if i < n_steps:
                twist_before_bc = twist_before_bc - twist_rate_lst[
                    i - 1] * dx  # compute the twist of each section between two points (positive for positive twist rate)
                twist_lst.append(twist_before_bc)

            if i == n_steps:  # this is the section where the boundary condition is applied
                twist_lst.append(
                    theta_0)  # now the section of the boundary condition is reached, this entire section attains this value (neglecting the twist along the even smaller subsection if point of boundary condition falls in between two points)
            if i > n_steps:
                twist_after_bc = twist_after_bc + twist_rate_lst[
                    i] * dx  # or -, plot if torque distribution is known. At the boundary condition, the sign of the twist should change
                twist_lst.append(twist_after_bc)

        return q1_lst, q2_lst, J, twist_rate_lst, twist_lst  # J, twist rate and twist at every x location taken

    def shear_stress_due_to_torsion(self):

        # input values
        t_sk = self.input['tsk'][self.aircraft]
        t_sp = self.input['tsp'][self.aircraft]
        q1, q2 = self.twist_of_aileron()[:2]

        tau_skin_cell_1_lst = [q1[i] / t_sk for i in range(len(q1))]
        tau_skin_cell_2_lst = [q2[i] / t_sk for i in range(len(q1))]
        tau_spar_lst = [(q2[i] - q1[i]) / t_sp for i in range(len(q1))]

        return tau_skin_cell_1_lst, tau_skin_cell_2_lst, tau_spar_lst

    def total_shear_stress(self):
        """Calculates the total shear stress distribution tau_yz per section for every x location along the span and puts this in a list, so 6 lists for each spanwise location
           it returns a list of lists with every list a shear stress distribution at a specific x location."""

        tau_1, tau_2, tau_3, tau_4, tau_5, tau_6 = self.shear_stress_due_to_shear()
        tau_skin_cell_1_lst, tau_skin_cell_2_lst, tau_spar_lst = self.shear_stress_due_to_torsion()

        total_shear_stress_distribution_at_every_x_loc = []

        for j in range(len(tau_skin_cell_1_lst)):

            tau_total_1_at_x_loc = [tau_1[i] + tau_skin_cell_1_lst[j] for i in range(len(tau_1))]
            tau_total_2_at_x_loc = [tau_2[i] + tau_spar_lst[j] for i in range(len(tau_2))]
            tau_total_3_at_x_loc = [tau_3[i] + tau_skin_cell_2_lst[j] for i in range(len(tau_3))]
            tau_total_4_at_x_loc = [tau_4[i] + tau_skin_cell_2_lst[j] for i in range(len(tau_4))]
            tau_total_5_at_x_loc = [tau_5[i] + tau_spar_lst[j] for i in range(len(tau_5))]
            tau_total_6_at_x_loc = [tau_6[i] + tau_skin_cell_1_lst[j] for i in range(len(tau_6))]

            total_shear_stress_distribution_at_every_x_loc.append(
                tau_total_1_at_x_loc + tau_total_2_at_x_loc + tau_total_3_at_x_loc + tau_total_4_at_x_loc + tau_total_5_at_x_loc + tau_total_6_at_x_loc)

        return total_shear_stress_distribution_at_every_x_loc  # example output[0] = shear stress distribution at first x location along the span, in the order of the sections 1,2,3,4,5,6

    def z_location(self):

        # input
        n1 = self.input['n_points'][self.aircraft]
        h = self.input['h'][self.aircraft]/2
        C_a = self.input['Ca'][self.aircraft]

        s_co2 = np.linspace(0, h, 2 * n1)

        ds = 0
        segment = [0, 1, 2]
        zco1 = []
        yco1 = []
        h_seg = 0.5 * np.pi / ((len(segment) * n1) - 1)
        for segment in segment:
            for i in range(n1):
                b = h_seg * (i + ds)
                zco1.append(-(h - h * np.cos(b)))
                yco1.append(h * np.sin(b))
            ds += n1

        zco6 = zco1[::-1]
        yco6 = [-i for i in yco1[::-1]]

        zco2 = 2 * n1 * [-h]
        yco2 = s_co2
        zco5 = zco2
        yco5 = np.linspace(0, -h, 2 * n1)
        zco3 = np.linspace(-h, -C_a, 7 * n1)
        yco3 = np.linspace(h, 0, 7 * n1)
        zco4 = zco3[::-1]
        yco4 = [-i for i in yco3[::-1]]

        return yco1, zco1, yco2, zco2, yco3, zco3, yco4, zco4, yco5, zco5, yco6, zco6

    def direct_stress_distribution(self):  # for a unit moment in x and y direction

        yco1, zco1, yco2, zco2, yco3, zco3, yco4, zco4, yco5, zco5, yco6, zco6 = self.z_location()
        Iyy = CrossSection(self.input, self.aircraft).get_moments_inertia()[2]
        zc = CrossSection(self.input, self.aircraft).get_centroid()[1]
        My = np.array(self.My)
        Mz = np.array(self.Mz)

        direct_stress_per_x = []

        for j in range(self.n_points):

            """Computes the Direct stress distrubtion along the cross-section at each point where the shear flow is calculated based on the Mx and My of a specific location along the span"""

            sigma_xx_1 = [My * (zco1[i] - zc) / Iyy + Mz * yco1[i] for i in tqdm(range(len(zco1)), desc="direct stress distribution")]
            sigma_xx_2 = [My * (zco2[i] - zc) / Iyy + Mz * yco2[i] for i in range(len(zco2))]
            sigma_xx_3 = [My * (zco3[i] - zc) / Iyy + Mz * yco3[i] for i in range(len(zco3))]
            sigma_xx_4 = [My * (zco4[i] - zc) / Iyy + Mz * yco4[i] for i in range(len(zco4))]
            sigma_xx_5 = [My * (zco5[i] - zc) / Iyy + Mz * yco5[i] for i in range(len(zco5))]
            sigma_xx_6 = [My * (zco6[i] - zc) / Iyy + Mz * yco6[i] for i in range(len(zco6))]

            direct_stress_per_x.append(sigma_xx_1 + sigma_xx_2 + sigma_xx_3 + sigma_xx_4 + sigma_xx_5 + sigma_xx_6)


        return direct_stress_per_x

    def von_mises_stress_distribution(self):  # use total shear stresses!
        """This function calculates the Von Mises stress distrubtion for every x-location taken
        along the span, based on the direct stress calculation and the total shear stress calculation"""

        # input
        direct_stress_distribution = self.direct_stress_distribution()
        shear_stress_distribution = self.total_shear_stress()
        T_lst = self.T  # TODO: CHANGE THIS
        sigma_vm_distribution_at_every_x_loc = []

        for j in range(self.n_points):
            sigma_vm = [np.sqrt(direct_stress_distribution[j-1][i] ** 2 + 3 * shear_stress_distribution[j-1][i] ** 2) for i in range(len(direct_stress_distribution[0]))]
            sigma_vm_distribution_at_every_x_loc.append(sigma_vm)

        return sigma_vm_distribution_at_every_x_loc  # example sigma_vm_distribution[0] is the distribution at the first x-location taken along the span

    def plot_shear_2d(self, x):

        yco1, zco1, yco2, zco2, yco3, zco3, yco4, zco4, yco5, zco5, yco6, zco6 = self.z_location()
        von_mises = np.array(self.von_mises_stress_distribution())

        for i, j in zip([yco2, yco3, yco4, yco5, yco6], [zco2, zco3, zco4, zco5, zco6]):

            yco1.extend(i)
            zco1.extend(j)

        plt.scatter(zco1, yco1, c=von_mises[x, :, x], cmap='jet')
        plt.xlabel('z[m]', fontsize=16)
        plt.ylabel('y[m]', fontsize=16)
        # plt.legend()
        c = plt.colorbar()
        c.set_label('Von Mises Stress [N/m^2]')
        plt.title('Von Mises Stress Distribution')
        plt.show()

        return 0

    def plot_shear_3d(self):

        yco1, zco1, yco2, zco2, yco3, zco3, yco4, zco4, yco5, zco5, yco6, zco6 = self.z_location()
        von_mis = np.array(self.von_mises_stress_distribution())

        l_a = self.input['la'][self.aircraft]
        x_axis = np.linspace(0, l_a, len(von_mis))

        for i, j in zip([yco2, yco3, yco4, yco5, yco6], [zco2, zco3, zco4, zco5, zco6]):

            yco1.extend(i)
            zco1.extend(j)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for i in range(len(von_mis)):
            p = ax.scatter(zco1, yco1, x_axis[i], c=von_mis[i, :, i], cmap='jet')

        c = fig.colorbar(p)
        c.set_label('Von Mises Stress [N/m^2]')
        ax.set_xlabel('Y [m]')
        ax.set_ylabel('Z [m]')
        ax.set_zlabel('Aileron Span [m]')
        plt.title('Von Mises Stress Distribution')

        plt.show()
        return 0

# debugging ========================================
DEBUG = True

if DEBUG:
    shear = ShearStress('A', 10)
    shear.plot_shear_3d()

