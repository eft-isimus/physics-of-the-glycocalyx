#!/usr/bin/env python3
# Usage: python create_job_files.py FOLDER WALLTIME N_DUMP_STEPS SIM_TYPE
# SIM_TYPE: brush or virus
# Note: This script should be placed in the same place as the analysis scripts
import sys
from pathlib import Path

INPUT = "input_file"
SCRIPTS = Path(__file__).resolve().parent

folder, walltime, n_steps, sim_type, job_name = Path(sys.argv[1]).resolve(), sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]

for sd in sorted(d for d in folder.glob("seed_*") if d.is_dir()):
    analysis = ""
    if sim_type == "brush":
        analysis = f""" && \\
module load compiler/python3.8
$HOME/venvs/lammps38/bin/python {SCRIPTS}/heights_calculate.py {sd}/eq_dump.lammpstraj && \\
$HOME/venvs/lammps38/bin/python {SCRIPTS}/heights_calculate.py {sd}/prod_dump.lammpstraj && \\
$HOME/venvs/lammps38/bin/python {SCRIPTS}/pz_calculate.py {sd}/prod_dump.lammpstraj {n_steps}"""
    elif sim_type == "virus":
        pass  # virus analysis scripts to be added later

    (sd / "job.pbs").write_text(f"""#!/bin/bash
#PBS -N {job_name}_{sd.name}
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
mpirun -np 8 lmp < {INPUT}{analysis}
""")
    print(f"wrote {sd / 'job.pbs'}")