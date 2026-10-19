#!/usr/bin/env python3
# Usage: python create_job_files.py FOLDER WALLTIME N_DUMP_STEPS
# Note: This script should be placed in the same place as the analysis scripts
import sys
from pathlib import Path

INPUT = "input_file"
SCRIPTS = Path(__file__).resolve().parent

folder, walltime, n_steps = Path(sys.argv[1]).resolve(), sys.argv[2], sys.argv[3]

for sd in sorted(d for d in folder.glob("seed_*") if d.is_dir()):
    (sd / "job.pbs").write_text(f"""#!/bin/bash
#PBS -N {sd.name}
#PBS -l nodes=1:ppn=8
#PBS -l walltime={walltime}
#PBS -l mem=4gb
#PBS -q cpu
#PBS -j oe
#PBS -o {sd}/job.log

cd {sd}
module load compiler/gcc-12.0
module load compiler/openmpi-3.1.6
module load apps/lammps-stable_2Aug2023
module load compiler/python3.8
mpirun -np 4 lmp_mpi -in {INPUT} && \\
python3 {SCRIPTS}/heights_calculate.py {sd}/eq_dump.lammpstraj && \\
python3 {SCRIPTS}/heights_calculate.py {sd}/prod_dump.lammpstraj && \\
python3 {SCRIPTS}/pz_calculate.py {sd}/prod_dump.lammpstraj {n_steps}
""")
    print(f"wrote {sd / 'job.pbs'}")