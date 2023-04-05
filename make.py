#! /usr/bin/python3

import make_helper
from make_helper import SpecializedEnvVarOption
from make_helper import SpecializedCMakeOption
from make_helper import SpecializedOptionType

specialized_opts = [

    SpecializedCMakeOption(config_var='bind-python', help_str='Build Python bindings for module',
                           option_type=SpecializedOptionType.CMakeBool,
                           modules_to_vars_map={'comp_mesh_core': ['SHEFF_BuildPython']}),

    SpecializedCMakeOption(config_var='bind-fortran', help_str='Build FORTRAN bindings for module',
                           option_type=SpecializedOptionType.CMakeBool,
                           modules_to_vars_map={'comp_mesh_core': ['SHEFF_BuildFortran']})
                           
]

make_helper.make(specialized_opts)
