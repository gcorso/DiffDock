"""
Symmetry-corrected RMSD calculations in Python
"""

if __name__ == "__main__":
    import argparse as ap
    import importlib.util
    import sys

    from spyrmsd import io
    from spyrmsd.rmsd import rmsdwrapper

    parser = ap.ArgumentParser(
        prog="python -m spyrmsd",
        description="Symmetry-corrected RMSD calculations in Python.",
    )

    parser.add_argument("reference", type=str, help="Reference file")
    parser.add_argument("molecules", type=str, nargs="+", help="Input file(s)")
    parser.add_argument("-m", "--minimize", action="store_true", help="Minimize (fit)")
    parser.add_argument(
        "-c", "--center", action="store_true", help="Center molecules at origin"
    )
    parser.add_argument("--hydrogens", action="store_true", help="Keep hydrogen atoms")
    parser.add_argument(
        "-n", "--nosymm", action="store_false", help="No graph isomorphism"
    )

    args = parser.parse_args()

    if (
        importlib.util.find_spec("openbabel") is None
        and importlib.util.find_spec("rdkit") is None
    ):
        raise ImportError(
            "OpenBabel or RDKit not found. Please install OpenBabel or RDKit to use sPyRMSD as a standalone tool."
        )

    try:
        ref = io.loadmol(args.reference)
    except OSError:
        print("ERROR: Reference file not found.", file=sys.stderr)
        exit(-1)

    # Load all molecules
    try:
        mols = [mol for molfile in args.molecules for mol in io.loadallmols(molfile)]
    except OSError:
        print("ERROR: Molecule file(s) not found.", file=sys.stderr)
        exit(-1)

    # Loop over molecules within fil
    RMSDlist = rmsdwrapper(
        ref,
        mols,
        symmetry=args.nosymm,  # args.nosymm store False
        center=args.center,
        minimize=args.minimize,
        strip=not args.hydrogens,
    )

    for RMSD in RMSDlist:
        print(f"{RMSD:.5f}")
