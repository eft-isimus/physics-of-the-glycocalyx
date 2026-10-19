import numpy as np
import sys
import os
import subprocess

# to parse inputs
def parse_array(s, dtype): 
    return np.array(s.strip("[]").split(","), dtype=dtype)

# inputs: total_sqrt_Np \ mixed \ g_rho \ Nm \ Np \ 
# bending_const \ Temperature \ dt \ eq_steps \ 
# eq_traj_dump_freq \ prod_steps \ 
# prod_traj_dump_freq \ output_dir \ n_seeds 

total_sqrt_Np = sys.argv[1]
mixed         = sys.argv[2]
g_rho         = sys.argv[3]
Nm            = sys.argv[4]
Np            = sys.argv[5]
bending_const = sys.argv[6]

Temperature         = sys.argv[7]
dt                  = sys.argv[8]
eq_steps            = sys.argv[9]
eq_traj_dump_freq   = sys.argv[10]
prod_steps          = sys.argv[11]
prod_traj_dump_freq = sys.argv[12]
output_dir          = sys.argv[13]
n_seeds             = int(sys.argv[14])

for seed in range(1, n_seeds + 1):
    seed_dir = os.path.join(output_dir, f"seed_{seed}")

    # create data file
    subprocess.run(
        [sys.executable, "create_brush_data_file.py",
        total_sqrt_Np, mixed, Nm, Np, g_rho, bending_const, seed_dir],
        check=True,
    )

    # create input files
    subprocess.run(
        [sys.executable, "create_input_files.py",
        bending_const, Temperature, dt, eq_steps, eq_traj_dump_freq, prod_steps, prod_traj_dump_freq, seed_dir],
        check=True,
    )