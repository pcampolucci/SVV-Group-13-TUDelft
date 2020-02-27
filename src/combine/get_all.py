"""
Initialize input values and loads in one go and get all the information
"""

from src.loads.send_loads import Loads
from src.input.input import Input
from src.combine.deflection import Deflection
from src.combine.max_stress import MaxStress
from src.combine.twist import Twist


class GetAll:

    def __init__(self, aircraft, steps, step_size):

        self.a = aircraft
        self.init = Loads(aircraft, steps, step_size)
        self.cross_section = Input(aircraft).cross_section_input()
        self.input_dict = Input(aircraft).get_discrete()

    def get_deflection(self):

        init_deflection = Deflection(self.init.geometry_input, self.init.discrete_input, self.init.point_loads,
                                     self.init.discrete_angles, self.init.discrete_deflections, self.init.step_size,
                                     self.init.la, self.init.steps, self.a)

        init_deflection.plot_deflection()

        return init_deflection

    def get_max_stress(self, plot=True):

        init_max_stress = MaxStress(self.cross_section, self.init.send_torque(), self.init.send_moment_y(),
                                    self.init.send_moment_z(), self.input_dict, self.a, self.init.steps)

        if plot:
            init_max_stress.plot_shear_3d()

        return init_max_stress

    def get_twist(self):

        twist_aileron = self.get_max_stress(False).twist_of_aileron()

        init_twist = Twist(twist_aileron, self.a, self.init.la, self.init.steps)

        init_twist.plot_twist()

        return init_twist
