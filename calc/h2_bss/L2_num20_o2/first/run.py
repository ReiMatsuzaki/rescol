# Common
arg_basis = ["-fem_type", "bss",
             "-bss_order", 2,
             "-bps_num_zs", 21,
             "-bps_zmax", 20.0,
             "-bps_type", "exp",
             "-y2s_rot", "sigma",
             "-y2s_parity", "gerade",
             "-y2s_mirror", "plus",
             "-y2s_lmax", 0]

arg_cscaling = ["-cscaling_type", "sharp_ecs",
                "-cscaling_r0", 40.0,
                "-cscaling_theta", 0.0]

# for he_guess
arg_guess = ["-guess_type", "fit", 
             "-z", 1.0, 
             "-vec_viewer", "binary:guess.dat"]
arg1 = ["../../../../build/complex/he_guess.out"] + arg_basis + arg_guess

# for h2mole.out
arg_eps = ["-eps_nev", 10,
           "-eps_ncv", 100,
           #        "-eps_mpd", 300, 
           "-eps_target", -4.0,
#           "-eps_type", "jd",
           "-eeps_view_values", "ascii:val.dat::append",
           "-init_path", "guess.dat",
           "-eps_error_backward", "::ascii_info_detail",
#           "-eps_conv_eig",
           "-eps_converged_reason"]

#           "-eps_error_backward", "::ascii_info_detail",
#           "-eps_conv_eig",

arg_sys = ["-bondlength", 0.0,
           "-num_bondlength", 2,
           "-d_bondlength", 0.5,
           "-malloc_dump"]

arg2 = ["../../../../build/complex/h2mole.out"] + arg_basis + arg_eps + arg_sys

import os
import subprocess
print arg1
subprocess.call(map(str, arg1))

os.system("echo > val.dat")
print arg2
subprocess.call(map(str, arg2))


