import numpy as np
import sys
import os

# to parse inputs
def parse_array(s, dtype): 
    return np.array(s.strip("[]").split(","), dtype=dtype)

# inputs: sqrt_Np mixed(True/False) [Nm array] [Np array] grafting_density [bending_const_array] output_dir
total_sqrt_Np = int(sys.argv[1])
mixed         = sys.argv[2].lower() == "true"
Nm            = parse_array(sys.argv[3], int)
Np            = parse_array(sys.argv[4], int)
g_rho         = float(sys.argv[5])
bending_const = parse_array(sys.argv[6], float)
output_dir    = sys.argv[7]

# calculating constants based on inputs
total_Np = total_sqrt_Np**2
box_len  = np.sqrt(total_Np/g_rho)

# to randomly distribute polymers of diff. lengths
sequence = np.repeat(np.arange(len(Nm)), Np)
np.random.shuffle(sequence)
sequence = sequence.reshape(total_sqrt_Np, total_sqrt_Np)


data_file = f"""
{np.sum(Nm * Np)} atoms 
{np.sum((Nm - 1) * Np)} bonds
{np.sum((Nm - 2) * Np)} angles

{len(Nm)*2 + 1} atom types
1 bond types
{len(Nm)} angle types
"""

bond_length = 1.12246                         # hard-coded here
z0          = 0.5                             # wall position, grafted atoms start here
zhi         = 1.5 * bond_length * np.max(Nm)
n_types     = 2*len(Nm) + 1

data_file += f"""
0 {box_len} xlo xhi
0 {box_len} ylo yhi
0 {zhi} zlo zhi

Masses

""" + "\n".join(f"{t} 1" for t in range(1, n_types + 1)) + "\n\nAtoms # molecular\n\n"

# addding Atoms, Bonds, and Angles sections
spacing = box_len / total_sqrt_Np
atom_lines, bond_lines, angle_lines = [], [], []
atom_id  = 1
mol_id   = 1
bond_id  = 1
angle_id = 1
for row in range(total_sqrt_Np):
    for col in range(total_sqrt_Np):
        t = sequence[row, col]
        x = (col + 0.5) * spacing
        y = (row + 0.5) * spacing
        first = atom_id                      # ID of the grafted atom of this polymer
        for k in range(Nm[t]):
            atom_type = 2*t + 1 if k == 0 else 2*t + 2
            atom_lines.append(f"{atom_id} {mol_id} {atom_type} {x} {y} {z0 + k*bond_length}")
            atom_id += 1
        for k in range(Nm[t] - 1):           # bonds: (k, k+1), last atom has no k+1
            bond_lines.append(f"{bond_id} 1 {first+k} {first+k+1}")
            bond_id += 1
        for k in range(Nm[t] - 2):           # angles: (k, k+1, k+2)
            angle_lines.append(f"{angle_id} {t+1} {first+k} {first+k+1} {first+k+2}")
            angle_id += 1
        mol_id += 1

data_file += "\n".join(atom_lines) + "\n"
data_file += "\nBonds\n\n" + "\n".join(bond_lines) + "\n"
data_file += "\nAngles\n\n" + "\n".join(angle_lines) + "\n"

os.makedirs(output_dir, exist_ok=True)
data_file_name = "data_file"
with open(os.path.join(output_dir, data_file_name), "w") as f:
    f.write(data_file)

# creating a log with the settings used
log_string = f"""
sqrt(Np) / grid length = {total_sqrt_Np}
grafting density: {g_rho}

Mixed: {mixed}
Nm array: {Nm}
Np array: {Np}

Bending constants: {bending_const}
"""
log_file_name = "data_log.txt"
with open(os.path.join(output_dir, log_file_name), "w") as f:
    f.write(log_string)