#!/usr/bin/env python3

# author: Ole Schuett

import argparse
from pathlib import Path
import numpy as np
import ase.io
from pao.io import AtomicKind, write_pao_file

# ======================================================================================
# Note that Oxygen has be first because it appears first in the coordinates.
KINDS = {
    "O": AtomicKind(
        name="O",
        atomic_number=8,
        prim_basis_name="DZVP-MOLOPT-GTH",
        prim_basis_size=13,
        pao_basis_size=4,
        nparams=13 * 4,
    ),
    "H": AtomicKind(
        name="H",
        atomic_number=1,
        prim_basis_name="DZVP-MOLOPT-GTH",
        prim_basis_size=5,
        pao_basis_size=2,
        nparams=5 * 2,
    ),
}


# ======================================================================================
def xyz2yzx(xyz):
    assert xyz.size == 3
    return np.array([xyz[1], xyz[2], xyz[0]])


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Synthesizes PAO basis for water molecules."
    )
    parser.add_argument("coords", type=Path)
    args = parser.parse_args()

    atoms = ase.io.read(args.coords)
    assert len(atoms) % 3 == 0
    n_waters = len(atoms) // 3

    xblocks = []
    for i in range(n_waters):
        # Check assumptions.
        assert atoms.symbols[i * 3 + 0] == "O"
        assert atoms.symbols[i * 3 + 1] == "H"
        assert atoms.symbols[i * 3 + 2] == "H"
        oxygen = atoms.positions[i * 3 + 0]
        hydrogen1 = atoms.positions[i * 3 + 1]
        hydrogen2 = atoms.positions[i * 3 + 2]

        # Construct xblock for Oxygen.
        oxygen_x = np.zeros([4, 13])
        oxygen_x[0, 0] = 1.0  # one s-function
        oxygen_x[1:4, 2:5] = np.eye(3)  # complete set of p-functions
        xblocks.append(oxygen_x)

        # Construct xblock for Hydrogens (both are getting the same block).
        molecular_normal = np.cross(oxygen - hydrogen1, oxygen - hydrogen2)
        hydrogen_x = np.zeros([2, 5])
        hydrogen_x[0, 0] = 1.0  # one s-function
        hydrogen_x[1, 2:5] = xyz2yzx(molecular_normal)  # single p-function for HOMO
        xblocks.append(hydrogen_x)
        xblocks.append(hydrogen_x)

    output_path = args.coords.with_name("synthetic.pao")
    write_pao_file(
        path=output_path,
        kinds=KINDS,
        atom2kind=atoms.symbols,
        coords=atoms.positions,
        xblocks=xblocks,
    )
    print(f"Wrote: {output_path}")


# ======================================================================================
main()

# EOF
