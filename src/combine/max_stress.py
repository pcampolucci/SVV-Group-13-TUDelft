"""
Title: Max Shear Stress tool
"""

import numpy as np
import matplotlib.pyplot as plt
from src.input.input import Input


class MaxStress:

    def __init__(self, cross_section, torque, moment_y, moment_z, input, aircraft, steps):

        self.cross_section = cross_section
        self.T = torque
        self.My = np.array(moment_y)
        self.Mz = np.array(moment_z)
        self.G = input['G'][aircraft]
        self.t_sk = input['tsk'][aircraft]
        self.t_sp = input['tsp'][aircraft]
        self.la = input['la'][aircraft]
        self.n1 = input['n_points'][aircraft]
        self.h = input['h'][aircraft]/2
        self.ca = input['Ca'][aircraft]
        self.q_lst = self.cross_section.get_shear_center()[1]
        self.i_yy = self.cross_section.get_moments_inertia()[2]
        self.zc = self.cross_section.get_centroid()[1]
        self.steps = steps//100

    # ============================================================================================================
    # GET AILERON TWIST

    def twist_of_aileron(self):

        q1_lst, q2_lst, J, twist_rate_lst, twist_lst = self.cross_section.twist_of_aileron(self.T, self.G)

        return q1_lst, q2_lst, J, twist_rate_lst, twist_lst  # J, twist rate and twist at every x location taken

    # ============================================================================================================
    # GET STRESS DUE TO TORSION AND SHEAR

    def shear_stress_due_to_shear(self):
        """ compute shear stress in skins and spar by dividing by thickness, note spar shear flow not included in shear_flow_lst
        Shear stress distribution per section due to Sy and Sz computed at the middle of each section """

        # input values
        qb1, qb2, qb3, qb4, qb5, qb6 = self.q_lst

        tau_1 = [qb1[i] / self.t_sk for i in range(len(qb1))]
        tau_2 = [qb2[i] / self.t_sp for i in range(len(qb2))]
        tau_3 = [qb3[i] / self.t_sk for i in range(len(qb3))]
        tau_4 = [qb4[i] / self.t_sk for i in range(len(qb4))]
        tau_5 = [qb5[i] / self.t_sp for i in range(len(qb5))]
        tau_6 = [qb6[i] / self.t_sk for i in range(len(qb6))]

        return tau_1, tau_2, tau_3, tau_4, tau_5, tau_6

    def shear_stress_due_to_torsion(self):

        # input values
        q1, q2 = self.twist_of_aileron()[:2]

        tau_skin_cell_1_lst = [q1[i] / self.t_sk for i in range(len(q1))]
        tau_skin_cell_2_lst = [q2[i] / self.t_sk for i in range(len(q1))]
        tau_spar_lst = [(q2[i] - q1[i]) / self.t_sp for i in range(len(q1))]

        return tau_skin_cell_1_lst, tau_skin_cell_2_lst, tau_spar_lst

    def total_shear_stress(self):
        """Calculates the total shear stress distribution tau_yz per section for every x location along the span
           and puts this in a list, so 6 lists for each spanwise location
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

        return total_shear_stress_distribution_at_every_x_loc

    # ============================================================================================================
    # GET VON MISES STRESS

    def z_location(self):

        # input

        s_co2 = np.linspace(0, self.h, 2 * self.n1)

        ds = 0
        segment = [0, 1, 2]
        zco1 = []
        yco1 = []
        h_seg = 0.5 * np.pi / ((len(segment) * self.n1) - 1)
        for sub_segment in segment:
            for i in range(self.n1):
                b = h_seg * (i + ds)
                zco1.append(-(self.h - self.h * np.cos(b)))
                yco1.append(self.h * np.sin(b))
            ds += self.n1

        zco6 = zco1[::-1]
        yco6 = [-i for i in yco1[::-1]]

        zco2 = 2 * self.n1 * [-self.h]
        yco2 = s_co2
        zco5 = zco2
        yco5 = np.linspace(0, -self.h, 2 * self.n1)
        zco3 = np.linspace(-self.h, -self.ca, 7 * self.n1)
        yco3 = np.linspace(self.h, 0, 7 * self.n1)
        zco4 = zco3[::-1]
        yco4 = [-i for i in yco3[::-1]]

        return yco1, zco1, yco2, zco2, yco3, zco3, yco4, zco4, yco5, zco5, yco6, zco6

    def direct_stress_distribution(self):  # for a unit moment in x and y direction

        yco1, zco1, yco2, zco2, yco3, zco3, yco4, zco4, yco5, zco5, yco6, zco6 = self.z_location()

        direct_stress_per_x = []

        for j in range(self.steps):
            """Computes the Direct stress distrubtion along the cross-section at each point
             where the shear flow is calculated based on the Mx and My of a specific location along the span"""

            sigma_xx_1 = [self.My * (float(zco1[i]) - self.zc) / self.i_yy + self.Mz * float(zco1[i]) for i in range(len(zco1))]
            sigma_xx_2 = [self.My * (float(zco2[i]) - self.zc) / self.i_yy + self.Mz * float(zco2[i]) for i in range(len(zco2))]
            sigma_xx_3 = [self.My * (float(zco3[i]) - self.zc) / self.i_yy + self.Mz * float(zco3[i]) for i in range(len(zco3))]
            sigma_xx_4 = [self.My * (float(zco4[i]) - self.zc) / self.i_yy + self.Mz * float(zco4[i]) for i in range(len(zco4))]
            sigma_xx_5 = [self.My * (float(zco5[i]) - self.zc) / self.i_yy + self.Mz * float(zco5[i]) for i in range(len(zco5))]
            sigma_xx_6 = [self.My * (float(zco6[i]) - self.zc) / self.i_yy + self.Mz * float(zco6[i]) for i in range(len(zco6))]

            direct_stress_per_x.append(sigma_xx_1 + sigma_xx_2 + sigma_xx_3 + sigma_xx_4 + sigma_xx_5 + sigma_xx_6)

        return direct_stress_per_x

    def von_mises_stress_distribution(self):  # use total shear stresses!
        """This function calculates the Von Mises stress distrubtion for every x-location taken
        along the span, based on the direct stress calculation and the total shear stress calculation"""

        # input
        direct_stress_distribution = self.direct_stress_distribution()
        shear_stress_distribution = self.total_shear_stress()
        sigma_vm_distribution_at_every_x_loc = []

        for j in range(self.steps):
            sigma_vm = [np.sqrt(direct_stress_distribution[j-1][i] ** 2 + 3 * shear_stress_distribution[j-1][i]
                                ** 2) for i in range(len(direct_stress_distribution[0]))]
            sigma_vm_distribution_at_every_x_loc.append(sigma_vm)

        return sigma_vm_distribution_at_every_x_loc

    # ============================================================================================================
    # PLOTTING

    def plot_shear_2d(self, x):

        yco1, zco1, yco2, zco2, yco3, zco3, yco4, zco4, yco5, zco5, yco6, zco6 = self.z_location()
        von_mises = np.array(self.von_mises_stress_distribution())

        for i, j in zip([yco2, yco3, yco4, yco5, yco6], [zco2, zco3, zco4, zco5, zco6]):

            yco1.extend(i)
            zco1.extend(j)

        plt.scatter(zco1[x], yco1[x], c=von_mises[x, :, x], cmap='jet')
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
        x_axis = np.linspace(0, self.la, len(von_mis))

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


