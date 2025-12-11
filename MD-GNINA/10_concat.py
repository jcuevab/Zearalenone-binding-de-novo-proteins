#!/usr/bin/env python3
# =============================================
# CONCATENACIÓN, CENTRADO Y ALINEAMIENTO
# DEL COMPLEJO PROTEÍNA–LIGANDO (robusto)
# =============================================

import os
import MDAnalysis as mda
import MDAnalysis.transformations as trans
from MDAnalysis.analysis import align
import parmed as pmd
import pytraj as pt

# ==== CONFIGURACIÓN ====
workDir = os.getcwd()
amber_prmtop = 'system_amber.prmtop'
amber_inpcrd = 'system_amber.inpcrd'
Equilibrated_PDB = 'prot_lig_equil.pdb'
Jobname = "prot_lig_prod"

first_stride = 1
Number_of_strides = 5
Skip = 1
Remove_waters = True
Output_format = "dcd"

# ==== ARCHIVOS ====
nstride = int(Number_of_strides)
output_prefix = f"{first_stride}-{first_stride + nstride - 1}"
template = os.path.join(workDir, f"{Jobname}_%s.{Output_format}")
flist = [template % str(i) for i in range(first_stride, first_stride + nstride)]

print("==============================================")
print("  CONCATENACIÓN, CENTRADO Y ALINEAMIENTO")
print("==============================================")
for f in flist:
    print("  -", os.path.basename(f))
print("==============================================")

# ==== TOPOLOGÍA ====
pprm = pmd.load_file(amber_prmtop, amber_inpcrd)
mprm_from_parmed = mda.Universe(pprm)
u = mprm_from_parmed.select_atoms('not (resname HOH)')
top_nw = u.convert_to('PARMED')

# ==== CARGA TRAJECTORIAS ====
u1 = mda.Universe(Equilibrated_PDB, flist)
u2 = mda.Universe(Equilibrated_PDB)

protein = u1.select_atoms("protein")
ligand = u1.select_atoms("not protein and not resname HOH")
complex_sel = protein + ligand

# ==== TRANSFORMACIONES CORRECTAS ====
# 1️⃣ Unwrap toda la caja completa
# 2️⃣ Centra el complejo (proteína + ligando)
# 3️⃣ (opcional) reenvuelve
print("> Aplicando unwrap y centrado del complejo (sin romper ligando)...")
transformations = [
    trans.unwrap(u1.atoms),
    trans.center_in_box(complex_sel, wrap=False),
    trans.wrap(u1.atoms)
]
u1.trajectory.add_transformations(*transformations)

# ==== ALINEAMIENTO ====
print("> Alineando según la proteína (ligando acompaña)...")
align.AlignTraj(u1, u2, select="name CA", in_memory=True).run()

# ==== EXPORTACIÓN ====
nw_dcd = os.path.join(workDir, f"{Jobname}_{output_prefix}_nw.{Output_format}")
nw_pdb = os.path.join(workDir, f"{Jobname}_nw.pdb")

if Remove_waters:
    sel = u1.select_atoms("not resname HOH")
    print("> Escribiendo trayectoria sin agua...")
    with mda.Writer(nw_dcd, sel.n_atoms) as W:
        for ts in u1.trajectory[::Skip]:
            W.write(sel)
    sel.write(nw_pdb)
    traj = nw_dcd
else:
    traj = os.path.join(workDir, f"{Jobname}_{output_prefix}_whole.{Output_format}")
    print("> Escribiendo trayectoria completa...")
    with mda.Writer(traj, u1.atoms.n_atoms) as W:
        for ts in u1.trajectory[::Skip]:
            W.write(u1)

# ==== VALIDACIÓN ====
if os.path.exists(traj):
    print(f"\n✅ Trayectorias concatenadas correctamente: {os.path.basename(traj)}")
else:
    print("\n❌ ERROR: No se generó el archivo final.")

# ==== PRUEBA CON PYTRAJ ====
print("\n> Probando carga en PyTraj...")
traj_load = pt.load(traj, Equilibrated_PDB)
print(traj_load)
print("\n🎯 ¡Concatenación, centrado y alineamiento del complejo completados exitosamente!")
