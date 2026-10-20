import numpy as np
import sys
import os

# to parse inputs
def parse_array(s, dtype): 
    return np.array(s.strip("[]").split(","), dtype=dtype)

# inputs: [bending_const_array] Temperature dt eq_steps eq_traj_dump_freq prod_steps prod_traj_dump_freq output_dir
bending_const = parse_array(sys.argv[1], float)
num_polymer_types = len(bending_const)
num_atom_types = int(2*num_polymer_types + 1)
Temperature = parse_array(sys.argv[2], float)[0]
dt = parse_array(sys.argv[3], float)[0]
eq_steps = parse_array(sys.argv[4], float)[0]
eq_traj_dump_freq = parse_array(sys.argv[5], float)[0]
prod_steps = parse_array(sys.argv[6], float)[0]
prod_traj_dump_freq = parse_array(sys.argv[7], float)[0]
output_dir    = sys.argv[8]

# NOTE: !!!
# There is currently only a WCA poly-poly 
# interaction, no separate virus potential
# NOTE: change this when adding viruses

input_string = f"""
#---------------
# Basic attributes of the simulation
#---------------
units lj
boundary p p f
dimension 3

#---------------
# styles
#---------------
atom_style molecular
bond_style fene
angle_style cosine/squared

#---------------
# pair
#---------------
pair_style lj/cut 1.12246   # cutoff = 2^(1/6) sigma
pair_modify shift yes               # energy shift

#---------------
# read-data
#---------------
read_data data_file

#---------------
# bond and angle coefficients
#---------------
bond_coeff 1 30.0 1.5 1.0 1.0   # 
"""

for k in range(num_polymer_types):
    input_string += f"""
angle_coeff {k+1} {bending_const[k]} 180.0
"""

input_string += f"""
#---------------
# pair coefficients
#---------------
"""

# NOTE: !!!
# The pair-coefficients section below 
# currently sets virus-polymer interactions 
# to be same as polymer-polymer int.s
# NOTE: change this when adding viruses

for i in range(num_atom_types):
    for j in range(num_atom_types-i):
        input_string += f"""
pair_coeff {i+1} {num_atom_types-j} 1.0 1 1.12246
"""

input_string += f"""
#---------------
# groups & other settings
#---------------
group base_atoms type"""

for i in range(int((num_atom_types-1)/2)):
    input_string+=f""" {2*i+1}"""

input_string += f"""
group mobile_atoms subtract all base_atoms

neighbor 1.0 bin
neigh_modify delay 100 every 5 check yes                                               # hard requirement for minimize
comm_modify cutoff 3.0
timestep {dt}

#---------------
# fixes and forces
#---------------
fix zwall all wall/reflect zlo EDGE zhi EDGE
velocity all create {Temperature} {np.random.randint(100_000, 999_999)} mom yes rot yes
fix 1 all langevin {Temperature} {Temperature} {100*dt} {np.random.randint(100_000, 999_999)}
fix 2 all nve
fix freeze_force base_atoms setforce 0.0 0.0 0.0                                       # fix just the base atoms to always be attached
fix_modify freeze_force energy no
velocity base_atoms set 0.0 0.0 0.0

# energy minimisation
minimize 1.0e-4 1.0e-6 1000 10000
reset_timestep 0

# equilibration phase
dump traj_dump all custom {eq_traj_dump_freq} eq_dump.lammpstraj id mol type x y z
dump_modify traj_dump sort id
run {eq_steps}      

undump traj_dump

# production phase
dump traj_dump all custom {prod_traj_dump_freq} prod_dump.lammpstraj id mol type x y z
dump_modify traj_dump sort id
run {prod_steps}      

"""

os.makedirs(output_dir, exist_ok=True)
input_file_name = "input_file"
with open(os.path.join(output_dir, input_file_name), "w") as f:
    f.write(input_string)

# creating a log with the settings used
log_string = f"""Temperature: {Temperature}
dt: {dt}
Bending constants: {bending_const}

Equilibriation steps: {eq_steps}
Equilibriation Time: {eq_steps*dt}
Equilibriation dump freq.: {eq_traj_dump_freq}

Production steps: {prod_steps}
Production Time: {prod_steps*dt}
Production dump freq.: {prod_traj_dump_freq}
"""
log_file_name = "input_log.txt"
with open(os.path.join(output_dir, log_file_name), "w") as f:
    f.write(log_string)