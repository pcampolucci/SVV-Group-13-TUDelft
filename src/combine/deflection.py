"""
Title: Max Shear Stress tool
"""

from src.input.cross_section.cross_section import CrossSection
from src.input.general.discrete_input import input_dict
from src.input.aero_load.aero_load import AeroLoad


class ShearStress:

    # Deflection due to Mz in y direction
    def __init__(self, dict, aircraft):
        self.section = CrossSection(dict, aircraft)
        self.i_zz = self.section.get_all()['i_zz']
        self.E = self.section.get_all()['E']
        self.P = self.section.get_all()['P']
        self.aero_load = AeroLoad(aircraft)

    def Defl_Mz(self, x):

        Izz = self.i_zz
        E = self.E
        qx = self.aero_load.get_q(x)
        P = self.P

        v = -1 / E / Izz * (-qx) + C1 * x + C2

        if x > x1:
            v += -1 / E / Izz * (F_y1 / 6 * (x - x1) ** (3))
        if x > x2 - xa / 2:
            v += -1 / E / Izz * (F_ya1 / 6 * (x - x2 + xa / 2) ** (3))
        if x > x2:
            v += -1 / E / Izz * (-F_y2 / 6 * (x - x2) ** (3))
        if x > x2 + xa / 2:
            v += -1 / E / Izz * (-P * np.sin(alpha) / 6 * (x - x2 - xa / 2) ** (3))
        if x > x3:
            v += -1 / E / Izz * (F_y3 / 6 * (x - x3) ** (3))

        return v

    def Defl_My(x):

        w = C3 * x + C4

        if x > x1:
            w += -1 / E / Izz * (F_z1 / 6 * (x - x1) ** (3))
        if x > x2 - xa / 2:
            w += -1 / E / Izz * (F_za1 / 6 * (x - x2 + xa / 2) ** (3))
        if x > x2:
            w += -1 / E / Izz * (F_z2 / 6 * (x - x2) ** (3))
        if x > x2 + xa / 2:
            w += -1 / E / Izz * (-P * np.cos(alpha) / 6 * (x - x2 - xa / 2) ** (3))
        if x > x3:
            w += -1 / E / Izz * (F_z3 / 6 * (x - x3) ** (3))

        return w