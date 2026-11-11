import sys
import numpy as np

path = sys.argv[1]
steps, heights = [], []

with open(path) as f:
    while f.readline():
        steps.append(int(f.readline()))
        f.readline()
        n = int(f.readline())
        for _ in range(4):
            f.readline()
        cols = f.readline().split()[2:]
        data = np.array([f.readline().split() for _ in range(n)], dtype=float)
        v = data[data[:, cols.index("mol")] == 0]
        v = v[np.argsort(v[:, cols.index("id")])]
        heights.append(v[:, cols.index("z")])
        virus_ids = v[:, cols.index("id")].astype(int)

header = "timestep," + ",".join(f"z_{i}" for i in virus_ids)
out = path.rsplit(".", 1)[0] + "_virus_heights.csv"
np.savetxt(out, np.column_stack([steps, heights]), delimiter=",", comments="",
           header=header, fmt=["%d"] + ["%.6f"] * len(virus_ids))