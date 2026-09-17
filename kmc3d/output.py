"""
output.py
=========
Writers that reproduce the EXACT output formats of the C++ code, so every
existing post-processing module (z-warp, gaussian density field, morphology
extraction, ML/SHAP pipelines) keeps working unchanged:

  output/trajectory/kmc-coords-<step>.xyz
      Extended-XYZ.  Header carries the 3x3 Lattice and
      Properties=species:S:1:pos:R:3 ; trailing columns: index, charge, name.

  Data2Excel.txt  /  Data2ExcelTrimmed.txt
      One whitespace row per counted kMC step with the same column order as
      print_coordinates / printStepfscreen2 in the original.
"""

from __future__ import annotations
import os
import numpy as np
from .lattice import CODE_NAME

# column header for the reaction log (documented; the C++ wrote no header row,
# we add one as a leading comment so downstream parsers can opt-in)
DATA_COLUMNS = [
    "stop", "MaxToStop", "currentStep", "cycleNumber", "stepCurrentV",
    "timeCurrentV", "currentTime", "presentV", "type", "Ospcs", "Ea",
    "expValue", "reactRate", "time", "RxnPlating", "RxnStripping",
    "RxnLiSurface", "RxnFSI", "RxnSFO", "RxnSOL", "RxnF5D", "RxnSOL2",
    "RxnPlatingSEI", "LiMetal", "LiIon", "LiMetalS", "LiIonS", "nO", "nF",
]


class Output:
    def __init__(self, root: str = "output", write_header: bool = True):
        self.root = root
        self.traj = os.path.join(root, "trajectory")
        os.makedirs(self.traj, exist_ok=True)
        self._data = open(os.path.join(root, "Data2Excel.txt"), "w")
        self._log = open(os.path.join(root, "kmc_info_log.txt"), "w")
        if write_header:
            self._data.write("# " + " ".join(DATA_COLUMNS) + "\n")

    # -- trajectory ---------------------------------------------------------
    def write_xyz(self, step: int, lat, occ, species_code, charge, code2sym):
        busy = np.where(occ == 1)[0]
        path = os.path.join(self.traj, f"kmc-coords-{step}.xyz")
        b = lat.box
        cart = lat.frac[busy] @ b
        with open(path, "w") as fh:
            fh.write(f"{busy.size}\n")
            fh.write(
                'Lattice="{:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} '
                '{:.6f} {:.6f} {:.6f}" '
                'Properties=species:S:1:pos:R:3 pbc="T T T"\n'.format(
                    b[0, 0], b[0, 1], b[0, 2], b[1, 0], b[1, 1], b[1, 2],
                    b[2, 0], b[2, 1], b[2, 2]))
            names = lat.name_code[busy]
            idx = lat.index1[busy]
            for k, s in enumerate(busy):
                sym = code2sym[species_code[s]]
                fh.write(f"{sym}{cart[k,0]:15.8f}{cart[k,1]:15.8f}"
                         f"{cart[k,2]:15.8f}{idx[k]:15d}{charge[s]:15d}"
                         f"{CODE_NAME[int(names[k])]:>15}\n")

    # -- reaction log -------------------------------------------------------
    def write_step(self, row: dict):
        vals = []
        for c in DATA_COLUMNS:
            v = row.get(c, 0)
            if isinstance(v, float):
                vals.append(f"{v:.6g}")
            else:
                vals.append(str(v))
        self._data.write(" ".join(vals) + "\n")
        self._data.flush()

    def log(self, msg: str):
        self._log.write(msg + "\n")
        self._log.flush()

    def close(self):
        self._data.close()
        self._log.close()
