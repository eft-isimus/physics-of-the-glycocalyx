import sys
from collections import deque
from itertools import islice
import numpy as np

traj = sys.argv[1]
N = int(sys.argv[2])
out = traj + "_pz.npz"

frames = deque(maxlen=N)  # keeps only the last N frames: (step, z sorted by id)

with open(traj) as f:
    while True:
        line = f.readline()
        if not line:
            break
        step = int(f.readline())
        f.readline()
        n = int(f.readline())
        for _ in range(4):
            f.readline()
        cols = f.readline().split()[2:]
        data = np.loadtxt(islice(f, n))

        o = np.argsort(data[:, cols.index("id")])  # fixed atom order in every frame
        typ = data[o, cols.index("type")].astype(int)
        z = data[o, cols.index("z")]
        frames.append((step, z))

# atom type 1,2 -> polymer type 1; 3,4 -> 2; ... (types are fixed, so use the last frame)
ptype = (typ + 1) // 2

res = {"timesteps": np.array([s for s, _ in frames])}
for k in np.unique(ptype):
    res[f"type_{k}"] = np.array([z[ptype == k] for _, z in frames])

np.savez(out, **res)