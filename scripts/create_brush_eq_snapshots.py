import os
import shutil
import sys

in_dir, n_snap, out_dir = sys.argv[1], int(sys.argv[2]), sys.argv[3]
path = os.path.join(in_dir, "prod_dump.lammpstraj")
os.makedirs(out_dir, exist_ok=True)
shutil.copyfile(os.path.join(in_dir, "data_file"), os.path.join(out_dir, "data_file"))

with open(path) as f:
    n_frames = sum(line.startswith("ITEM: TIMESTEP") for line in f)
if n_snap > n_frames:
    sys.exit(f"Only {n_frames} frames in {path}")

chosen = {int((k + 0.5) * n_frames / n_snap): k + 1 for k in range(n_snap)}

with open(path) as f:
    for i in range(n_frames):
        header = [f.readline() for _ in range(9)]
        atoms = [f.readline() for _ in range(int(header[3]))]
        if i in chosen:
            print(f"snapshot {chosen[i]}: timestep {header[1].strip()}")
            header[1] = "0\n"
            out = os.path.join(out_dir, f"eq_snapshot_{chosen[i]}.lammpstraj")
            with open(out, "w") as g:
                g.writelines(header + atoms)