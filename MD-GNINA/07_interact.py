import MDAnalysis as mda
import prolif as plf
import numpy as np
import os
from prolif.plotting.network import LigNetwork

workDir = os.getcwd()
u = mda.Universe("system_amber.prmtop", "system.pdb")
lig = u.select_atoms("not protein and not (resname HOH) and not (resname CL) and not (resname NA) and not (resname K)")
prot = u.select_atoms("protein")

# create RDKit-like molecules for visualisation
lmol = plf.Molecule.from_mda(lig)
pmol = plf.Molecule.from_mda(prot)

fp = plf.Fingerprint()
fp.run(u.trajectory[::10], lig, prot)
df = fp.to_dataframe(return_atoms=True)

net = LigNetwork.from_ifp(df, lmol,
                          # replace with `kind="frame", frame=0` for the other depiction
                          kind="aggregate", frame=0,
                          rotation=270)
net.save(os.path.join(workDir, "initial.html"))
