import numpy as np
import sys
import os
import shutil

# to parse inputs
def parse_array(s, dtype):
    return np.array(s.strip("[]").split(","), dtype=dtype)

# inputs: [bending_const_array] Temperature dt n_steps dump_freq virus_diameter n_virus force input_dir output_dir [seed_numbers]
# input_dir must contain data_file and eq_snapshot_{seed_number}.lammpstraj for every seed number
# virus_diameter = contact diameter, i.e. distance at which two viruses touch (WCA minimum)
# force is the magnitude of the constant force on each virus, applied in -z
bending_const = parse_array(sys.argv[1], float)
num_polymer_types = len(bending_const)
num_atom_types = int(2*num_polymer_types + 1)
virus_type = num_atom_types
Temperature = parse_array(sys.argv[2], float)[0]
dt = parse_array(sys.argv[3], float)[0]
n_steps = parse_array(sys.argv[4], int)[0]
dump_freq = parse_array(sys.argv[5], int)[0]
d_v = parse_array(sys.argv[6], float)[0]
n_virus = parse_array(sys.argv[7], int)[0]
force = parse_array(sys.argv[8], float)[0]
input_dir = sys.argv[9]
output_dir = sys.argv[10]
seed_numbers = parse_array(sys.argv[11], int)

sigma_p = 1.0                # polymer atom WCA sigma
rcut = 1.12246               # WCA cutoff = 2^(1/6) sigma = contact distance
d_p = rcut*sigma_p           # polymer atom (contact) diameter
gap = 1.0                    # polymer atom tops -> virus bottom, and virus top -> wall
virus_mass = d_v/d_p         # linear in diameter, mass 1 for a polymer atom
dump_group = "all"           # change to "virus_atoms" to dump only the viruses

# WCA sigma for a pair of atom types (i <= j), so that the WCA minimum = contact distance
def sigma(i, j):
    if j == virus_type and i == virus_type:
        return d_v/rcut                  # virus-virus: contact distance d_v
    if j == virus_type:
        return (d_p + d_v)/2/rcut        # virus-polymer: contact distance (d_p + d_v)/2
    return sigma_p                       # polymer-polymer

# read timestep, box bounds and max z from a single-frame dump (id mol type x y z)
def read_snapshot(path):
    with open(path) as f:
        lines = f.read().splitlines()
    timestep = int(lines[1])
    n_atoms = int(lines[3])
    box = [[float(v) for v in lines[5+k].split()[:2]] for k in range(3)]
    z_col = lines[8].split()[2:].index("z")
    z_max = max(float(l.split()[z_col]) for l in lines[9:9+n_atoms])
    return timestep, box, z_max

# evenly spaced xy grid (periodic box), same height for all viruses
def place_viruses(n, d, xlo, xhi, ylo, yhi):
    Lx, Ly = xhi - xlo, yhi - ylo
    nx = min(n, int(np.ceil(np.sqrt(n*Lx/Ly))))
    ny = int(np.ceil(n/nx))
    dx, dy = Lx/nx, Ly/ny
    if n > 1 and min(dx, dy) < d:
        raise RuntimeError(f"Box too small for {n} viruses of diameter {d} (grid spacing {min(dx, dy):.2f})")
    grid = [(xlo + (i + 0.5)*dx, ylo + (j + 0.5)*dy) for j in range(ny) for i in range(nx)]
    return grid[:n]

# parts of the input file which are the same for every seed
angle_string = ""
for k in range(num_polymer_types):
    angle_string += f"angle_coeff {k+1} {bending_const[k]} 180.0\n"

pair_string = ""
for i in range(1, num_atom_types+1):
    for j in range(i, num_atom_types+1):
        s = sigma(i, j)
        pair_string += f"pair_coeff {i} {j} 1.0 {s:.5f} {rcut*s:.5f}\n"

base_types = " ".join(str(2*i+1) for i in range(num_polymer_types))

for seed in seed_numbers:
    seed_dir = os.path.join(output_dir, f"seed_{seed}")
    os.makedirs(seed_dir, exist_ok=True)
    snapshot = f"eq_snapshot_{seed}.lammpstraj"
    shutil.copy(os.path.join(input_dir, "data_file"), seed_dir)
    shutil.copy(os.path.join(input_dir, snapshot), seed_dir)

    # brush extent from the copied snapshot, and virus / wall placement
    snap_step, box, z_max = read_snapshot(os.path.join(seed_dir, snapshot))
    (xlo, xhi), (ylo, yhi), (zlo, _) = box
    z_brush_top = z_max + d_p/2                 # top of the highest polymer atom
    z_v = z_brush_top + gap + d_v/2             # virus centre height
    z_wall = z_v + d_v/2 + gap                  # wall above the viruses
    rng = np.random.default_rng(int(seed))      # LAMMPS random seeds
    positions = place_viruses(n_virus, d_v, xlo, xhi, ylo, yhi)

    virus_string = ""
    for x, y in positions:
        virus_string += f"create_atoms {virus_type} single {x:.4f} {y:.4f} {z_v:.4f} units box\n"

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
pair_style lj/cut {rcut}   # cutoff = 2^(1/6) sigma
pair_modify shift yes               # energy shift

#---------------
# read-data, equilibrated brush, viruses
#---------------
read_data data_file
read_dump {snapshot} {snap_step} x y z box no
change_box all z final {zlo} {z_wall:.4f} units box     # box top = wall above the viruses
mass {virus_type} {virus_mass:.5f}                      # virus mass, linear in diameter
{virus_string}
#---------------
# bond and angle coefficients
#---------------
bond_coeff 1 30.0 1.5 1.0 1.0   #

{angle_string}
#---------------
# pair coefficients (WCA, virus contact diameter {d_v})
#---------------
{pair_string}
#---------------
# groups & other settings
#---------------
group base_atoms type {base_types}
group mobile_atoms subtract all base_atoms
group virus_atoms type {virus_type}

neighbor 1.0 bin
neigh_modify delay 0 every 1 check yes
comm_modify cutoff 3.0
timestep {dt}

#---------------
# fixes and forces
#---------------
fix zwall all wall/reflect zlo EDGE zhi EDGE
fix zwall_virus virus_atoms wall/reflect zhi {z_wall - d_v/2:.4f} units box     # virus surface, not centre, stops at the wall
velocity all create {Temperature} {rng.integers(100_000, 999_999)} mom yes rot yes
fix 1 all langevin {Temperature} {Temperature} {100*dt} {rng.integers(100_000, 999_999)}
fix 2 all nve
fix freeze_force base_atoms setforce 0.0 0.0 0.0                                       # fix just the base atoms to always be attached
fix_modify freeze_force energy no
velocity base_atoms set 0.0 0.0 0.0

# energy minimisation
minimize 1.0e-4 1.0e-6 1000 10000
reset_timestep 0

# constant downward force on the viruses
fix push virus_atoms addforce 0.0 0.0 {-abs(force)}

# run
fix prog all print {max(n_steps//10, 1)} "RUN $(100*elapsed/{n_steps}:%.0f)% done, step $(step), est. time left $(cpuremain:%.0f) s" file progress.txt screen no
dump traj_dump {dump_group} custom {dump_freq} dump.lammpstraj id mol type x y z
dump_modify traj_dump sort id
run {n_steps}
unfix prog

"""

    with open(os.path.join(seed_dir, "input_file"), "w") as f:
        f.write(input_string)

    # creating a log with the settings used
    log_string = f"""Temperature: {Temperature}
dt: {dt}
Bending constants: {bending_const}

Steps: {n_steps}
Time: {n_steps*dt}
Dump freq.: {dump_freq}
Dump group: {dump_group}

Virus contact diameter: {d_v}
Virus mass: {virus_mass}
Number of viruses: {n_virus}
Force per virus (-z): {abs(force)}

Snapshot: {snapshot} (timestep {snap_step})
Max polymer atom z (centre): {z_max}
Brush top (polymer atom top): {z_brush_top}
Virus centre z: {z_v}
Wall z: {z_wall}
Virus xy positions: {[(round(x, 3), round(y, 3)) for x, y in positions]}
"""
    with open(os.path.join(seed_dir, "input_log.txt"), "w") as f:
        f.write(log_string)