"""
Title: Cross Section Information

Author: Pietro Campolucci
"""

# import packages
import numpy as np
import matplotlib.pyplot as plt

# debugging setting
DEBUG = False
from src.input.general.discrete_input import input_dict


class CrossSection:
    """ get all the geometrical information for the aircraft chosen and output properties """

    def __init__(self, input_dict, aircraft):
        self.input = input_dict
        self.aircraft = aircraft

    def stiffener_spacing(self):

        # input required
        h = self.input['h'][self.aircraft]
        C_a = self.input['Ca'][self.aircraft]
        spacing_circle = self.input['nst_circle'][self.aircraft] + 1
        spacing_triangle = self.input['nst_triangle'][self.aircraft] + 2

        # function
        perimeter_semicircle = np.pi * h
        perimeter_triangle = 2 * np.sqrt((C_a - h) * 2 + h * 2)
        stiffener_spacing_semicircle = perimeter_semicircle / spacing_circle
        stiffener_spacing_triangle = perimeter_triangle / spacing_triangle

        return stiffener_spacing_triangle, stiffener_spacing_semicircle, perimeter_semicircle, perimeter_triangle

    def boom_locations(self):

        # input required
        h = self.input['h'][self.aircraft]
        C_a = self.input['Ca'][self.aircraft]
        stiffener_spacing_semicircle = self.stiffener_spacing()[1]
        spacing_circle = self.input['nst_circle'][self.aircraft] + 1
        spacing_triangle = self.input['nst_triangle'][self.aircraft] + 2

        # function
        y_lst = []
        z_lst = []

        for i in range(spacing_circle//2):
            theta = np.radians(i * stiffener_spacing_semicircle / (np.pi * h) * 180)
            delta_y = h * np.sin(theta)
            delta_z = h - h * np.cos(theta)
            y_lst.append(delta_y)
            z_lst.append(-delta_z)

        for i in range(1, spacing_triangle//2):
            y = h - i * h / (spacing_triangle//2)
            z = -h - i * (C_a - h) / (spacing_triangle//2)
            y_lst.append(y)
            z_lst.append(z)

        reverse_z = z_lst[::-1]
        reverse_y = y_lst[::-1]

        for i in range(len(reverse_z) - 1):
            z_lst.append(reverse_z[i])
            y_lst.append(-reverse_y[i])

        return z_lst, y_lst

    def stiffener_area(self):

        # input required
        w = self.input['wst'][self.aircraft]
        h = self.input['hst'][self.aircraft]
        t = self.input['tst'][self.aircraft]

        # function
        area = w * t + (h-t) * t

        return area

    def centroid(self):

        # input required
        C_a = self.input['Ca'][self.aircraft]
        t_sk = self.input['tsk'][self.aircraft]
        h = self.input['h'][self.aircraft]
        z_lst = self.boom_locations()[0]
        A_st = self.stiffener_area()
        n_st = self.input['nst_circle'][self.aircraft] + self.input['nst_triangle'][self.aircraft]

        # function
        yc = 0
        stff = sum([i * A_st for i in z_lst])
        triangle = -2 * (h + (C_a - h) / 2) * np.sqrt((C_a - h) * 2 + h * 2) * t_sk
        semic = -h * (1 - 2 / np.pi) * np.pi*(h * 2 - (h - t_sk) * 2) / 2
        spar = -t_sk * h * 2 * h

        zc = (stff + semic + triangle + spar) / (n_st * A_st + np.sqrt((C_a - h) * 2 + h * 2) * t_sk * 2 + np.pi*(
            h * 2 - (h - t_sk) * 2) / 2 + t_sk * 2 * h)

        return yc, zc

    def MMI(self):

        # input information
        C_a = self.input['Ca'][self.aircraft]
        t_sk = self.input['tsk'][self.aircraft]
        h = self.input['h'][self.aircraft]
        la = self.input['la'][self.aircraft]
        z_lst = self.boom_locations()[0]
        y_lst = self.boom_locations()[1]
        A_st = self.stiffener_area()
        yc = self.centroid()[0]
        zc = self.centroid()[1]

        # functions
        inclined_section = np.sqrt((C_a - h) * 2 + h * 2)

        stff_z = sum([(i - yc) ** 2 * A_st for i in y_lst])
        triangle_z = ((2 * inclined_section) * 3 * t_sk / 12)*(
            h / inclined_section) * 2 + 2 * inclined_section * t_sk * (h / 2) * 2
        spar_z = (2 * h / 12)*t_sk * 3
        semic_z = np.pi * h ** 3 * t_sk / 2
        Izz = stff_z + triangle_z + spar_z + semic_z

        stff_y = sum([(i - zc) ** 2 * A_st for i in z_lst])
        triangle_y = ((2 * inclined_section) * 3 * t_sk / 12) * (
                    (C_a - h) / inclined_section) * 2 + 2 * inclined_section * t_sk * (-(h + (C_a - h) / 2) - zc) ** 2
        spar_y = (t_sk / 12) * (2 * h) * 3 + 2 * h * t_sk*(-h - zc) ** 2
        semic_y = np.pi * h * 3 * t_sk / 2 + (-h*(1 - 2 / np.pi) - zc) * 2 * np.pi*(h * 2 - (h - t_sk) * 2) / 2
        Iyy = stff_y + triangle_y + spar_y + semic_y

        Ixx = 1 / 12 * la * C_a ** 3

        return Izz, Ixx, Iyy

    def plot_cross_section(self):
        x_axis = self.boom_locations()[0]
        y_axis = self.boom_locations()[1]
        centroid = self.centroid()

        # plot boom points
        plt.figure(1)
        plt.title("Aileron Cross Section Plot")
        plt.scatter(x_axis, y_axis, color='g', label='stiffeners')
        plt.scatter(centroid[1], centroid[0], color='r', label='centroid')
        plt.plot([0, 0], [-0.25, 0.25], color='k', alpha=0.4)
        plt.plot([-0.55, 0.05], [0, 0], color='k', alpha=0.4)
        plt.legend()
        plt.show()

        return 0

    def get_shear_center(self):
        return 0

    def get_all(self):
        # return a cross section report as dict
        all_dict = {
            'area_st': CrossSection.stiffener_area(self),
            'SC': CrossSection.get_shear_center(self),
            'Centroid': CrossSection.centroid(self),
            'Izz': CrossSection.MMI(self)[0],
            'Iyy': CrossSection.MMI(self)[1],
            'Ixx': CrossSection.MMI(self)[2]
        }

        return all_dict

    def get_report(self):
        print("\nCross section input information for aircraft aileron\n")
        print(self.get_all())
        print("\nPlotting ...\n")
        self.plot_cross_section()
        return 0

# debugging
if DEBUG:
    cross_section = CrossSection(input_dict, 'A')
    cross_section.get_report()

