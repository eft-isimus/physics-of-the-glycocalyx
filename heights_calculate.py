import sys
from itertools import islice
import numpy as np

traj = sys.argv[1]
out = traj + "_heights.npz"

steps = []
heights = {}  # polymer type (1, 2, ...) -> list of per-frame arrays

with open(traj) as f:
    while True:
        line = f.readline()
        if not line:
            break
        steps.append(int(f.readline()))
        f.readline()
        n = int(f.readline())
        for _ in range(4):  # BOX BOUNDS header + 3 bound lines
            f.readline()
        cols = f.readline().split()[2:]  # column names after "ITEM: ATOMS"
        data = np.loadtxt(islice(f, n))

        aid = data[:, cols.index("id")]
        mol = data[:, cols.index("mol")]
        typ = data[:, cols.index("type")].astype(int)
        z = data[:, cols.index("z")]

        # last atom of each polymer = highest id within each mol
        o = np.lexsort((aid, mol))
        is_last = np.r_[mol[o][1:] != mol[o][:-1], True]
        sel = o[is_last]
        t, zl = typ[sel], z[sel]

        # polymer atoms (even types only): atom type 2k -> polymer type k
        for tt in np.unique(t[t % 2 == 0]):
            heights.setdefault(tt // 2, []).append(zl[t == tt])

np.savez(out, timesteps=np.array(steps),
         **{f"type_{k}": np.array(v) for k, v in sorted(heights.items())})