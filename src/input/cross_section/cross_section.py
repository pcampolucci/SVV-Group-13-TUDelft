"""
Title: Cross Section Information

Author: Pietro Campolucci
"""

# import packages

class CrossSection:
    """ get all the geometrical information for the aircraft chosen and output properties """

    def __init__(self, input_dict, aircraft):
        self.input = input_dict
        self.aircraft = aircraft

    def get_Ixx(self):
        return 0

    def get_Iyy(self):
        return 0

    def get_shear_center(self):
        return 0

    def get_centroid(self):
        return 0

    def get_all(self):
        all_dict = {
            'Ixx': CrossSection.get_Ixx(self),
            'Izz': CrossSection.get_Iyy(self),
            'SC': CrossSection.get_shear_center(self),
            'CG': CrossSection.get_centroid(self)
        }

        return all_dict