"""
Title: Max Shear Stress tool
"""

from src.input.cross_section.cross_section import CrossSection
from src.input.general.discrete_input import input_dict


class ShearStress:

    def __init__(self):
        self.input = input_dict
        self.aircraft = 'A'
        self.q_1 = CrossSection(input_dict, self.aircraft).twist_of_aileron([1], self.input['G'][self.aircraft])[0]
        self.q_2 = CrossSection(input_dict, self.aircraft).twist_of_aileron([1], self.input['G'][self.aircraft])[1]
        self.q_lst = CrossSection(input_dict, self.aircraft).get_shear_center()[1]

    def shear_stress_due_to_shear(self):
        """compute shear stress in skins and spar by dividing by thickness,
        note spar shear flow not included in shear_flow_lst. Shear_flow_lsts are
        the shear flow calculated at all intermediate points between the booms"""

        # input values
        t_sk = self.input['tsk'][self.aircraft]
        t_sp = self.input['tsp'][self.aircraft]
        qb1, qb2, qb3, qb4, qb5, qb6 = self.q_lst[0], self.q_lst[1], self.q_lst[2], self.q_lst[3], self.q_lst[4], self.q_lst[5]

        # function
        tau_1 = [qb1[i] / t_sk for i in range(len(qb1))]
        tau_2 = [qb2[i] / t_sp for i in range(len(qb2))]
        tau_3 = [qb3[i] / t_sk for i in range(len(qb3))]
        tau_4 = [qb4[i] / t_sk for i in range(len(qb4))]
        tau_5 = [qb5[i] / t_sp for i in range(len(qb5))]
        tau_6 = [qb6[i] / t_sk for i in range(len(qb6))]

        return tau_1, tau_2, tau_3, tau_4, tau_5, tau_6


    def shear_stress_due_to_torsion(self):

        # input values
        t_sk = self.input['tsk'][self.aircraft]
        t_sp = self.input['tsp'][self.aircraft]

        # function
        tau_skin_cell_1 = self.q_1 / t_sk
        tau_skin_cell_2 = self.q_2 / t_sk
        tau_spar = (self.q_2 - self.q_1) / t_sp

        return tau_skin_cell_1, tau_skin_cell_2, tau_spar

    def total_shear_stress(self):
        """Calculate the total shear stress by combining the shear stress due to shear and torsion"""

        # input values
        tau_1, tau_2, tau_3, tau_4, tau_5, tau_6 = self.shear_stress_due_to_shear()
        tau_skin_cell_1, tau_skin_cell_2, tau_spar = self.shear_stress_due_to_torsion()

        # function
        tau_total_1 = [tau_1[i] + tau_skin_cell_1 for i in range(len(tau_1))]
        tau_total_2 = [tau_2[i] + tau_spar for i in range(len(tau_2))]
        tau_total_3 = [tau_3[i] + tau_skin_cell_2 for i in range(len(tau_3))]
        tau_total_4 = [tau_4[i] + tau_skin_cell_2 for i in range(len(tau_4))]
        tau_total_5 = [tau_5[i] + tau_spar for i in range(len(tau_5))]
        tau_total_6 = [tau_6[i] + tau_skin_cell_1 for i in range(len(tau_6))]

        return tau_total_1, tau_total_2, tau_total_3, tau_total_4, tau_total_5, tau_total_6

# debugging ========================================
DEBUG = True

if DEBUG:
    shear = ShearStress().total_shear_stress()
