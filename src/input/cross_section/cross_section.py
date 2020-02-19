"""
Title: Cross Section Information

Author: Pietro Campolucci
"""

# import packages
import numpy as np

# debugging setting
DEBUG = True


class CrossSection:
    """ get all the geometrical information for the aircraft chosen and output properties """

    def __init__(self, input_dict, aircraft):
        self.input = input_dict
        self.aircraft = aircraft

    def stiffener_spacing(self):

        # input required
        h = self.input['h'][self.aircraft]
        C_a = self.input['C_a'][self.aircraft]

        # function
        perimeter_semicircle = np.pi * h
        perimeter_triangle = 2 * np.sqrt((C_a - h) * 2 + h * 2)
        stiffener_spacing_semicircle = perimeter_semicircle / 6
        stiffener_spacing_triangle = perimeter_triangle / 14

        return stiffener_spacing_triangle, stiffener_spacing_semicircle, perimeter_semicircle, perimeter_triangle

    def boom_locations(self):

        # input required
        h = self.input['h'][self.aircraft]
        C_a = self.input['C_a'][self.aircraft]
        stiffener_spacing_semicircle = self.stiffener_spacing()[1]

        # function
        y_lst = []
        z_lst = []

        for i in range(3):
            theta = np.radians(i * stiffener_spacing_semicircle / (np.pi * h) * 180)
            delta_y = h * np.sin(theta)
            delta_z = h - h * np.cos(theta)
            y_lst.append(delta_y)
            z_lst.append(-delta_z)
        for i in range(1, 7):
            y = h - i * h / 7
            z = -h - i * (C_a - h) / 7
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
        w_st = self.input['w_st'][self.aircraft]
        h_st = self.input['h_st'][self.aircraft]
        t_st = self.input['t_st'][self.aircraft]

        # function
        area_st = w_st * t_st + (h_st - t_st) * t_st

        return area_st

    def centroid(self):

        # input required
        C_a = self.input['C_a'][self.aircraft]
        t_sk = self.input['t_sk'][self.aircraft]
        h = self.input['h'][self.aircraft]
        z_lst = self.boom_locations()[0]
        A_st = self.stiffener_area()

        # function
        yc = 0
        stff = sum([i * A_st for i in z_lst])
        triangle = -2 * (h + (C_a - h) / 2) * np.sqrt((C_a - h) * 2 + h * 2) * t_sk
        semic = -h * (1 - 2 / np.pi) * np.pi(h * 2 - (h - t_sk) * 2) / 2
        spar = -t_sk * h * 2 * h

        zc = (stff + semic + triangle + spar) / (17 * A_st + np.sqrt((C_a - h) * 2 + h * 2) * t_sk * 2 + np.pi(
            h * 2 - (h - t_sk) * 2) / 2 + t_sk * 2 * h)

        return yc, zc

    def MMI(self): #z_lst, y_lst, A_st, C_a, t_sk, h, zc, yc):

        # input information
        C_a = self.input['C_a'][self.aircraft]
        t_sk = self.input['t_sk'][self.aircraft]
        h = self.input['h'][self.aircraft]
        z_lst = self.boom_locations()[0]
        y_lst = self.boom_locations()[1]
        A_st = self.stiffener_area()
        yc = self.centroid()[0]
        zc = self.centroid()[1]

        # functions
        inclined_section = np.sqrt((C_a - h) * 2 + h * 2)

        stff_z = sum([(i - yc) ** 2 * A_st for i in y_lst])
        triangle_z = ((2 * inclined_section) * 3 * t_sk / 12)(
            (h) / inclined_section) * 2 + 2 * inclined_section * t_sk * (h / 2) * 2
        spar_z = (2 * h / 12)(t_sk) * 3
        semic_z = np.pi * h ** 3 * t_sk / 2
        Izz = stff_z + triangle_z + spar_z + semic_z

        stff_y = sum([(i - zc) ** 2 * A_st for i in z_lst])
        triangle_y = ((2 * inclined_section) * 3 * t_sk / 12) * (
                    (C_a - h) / inclined_section) * 2 + 2 * inclined_section * t_sk * (-(h + (C_a - h) / 2) - zc) ** 2
        spar_y = (t_sk / 12) * (2 * h) * 3 + 2 * h * t_sk(-h - zc) ** 2
        semic_y = np.pi * h * 3 * t_sk / 2 + (-h(1 - 2 / np.pi) - zc) * 2 * np.pi(h * 2 - (h - t_sk) * 2) / 2
        Iyy = stff_y + triangle_y + spar_y + semic_y

        Ixx = 1 / 12 * l_a * C_a ** 3

        return Izz, Ixx, Iyy

    def get_Ixx(self):
        return 0

    def get_Iyy(self):
        return 0

    def get_shear_center(self):
        return 0

    def get_all(self):
        all_dict = {
            'Ixx': CrossSection.get_Ixx(self),
            'Izz': CrossSection.get_Iyy(self),
            'SC': CrossSection.get_shear_center(self),\
        }

        return all_dict