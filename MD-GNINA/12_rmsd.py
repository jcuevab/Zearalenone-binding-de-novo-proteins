import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ===============================
# CONFIGURACIÓN
# ===============================
topology = "system_amber.prmtop"
trajectory = "prot_lig_prod_1-5_whole.dcd"
ref_structure = "system.pdb"   # referencia para alinear SOLO la proteína

protein_sel = "backbone"
ligand_sel = "resname UNK and not name H*"

output_csv = "ligand_rmsd_noFit_frame0.csv"
output_plot = "ligand_rmsd_noFit_frame0.png"

print("=== RMSD del ligando SIN FIT usando frame 0 como referencia ===\n")

# ===============================
# CARGA UNIVERSOS
# ===============================
u = mda.Universe(topology, trajectory)          # trayectoria
ref = mda.Universe(topology, ref_structure)     # referencia externa para proteína

ligand = u.select_atoms(ligand_sel)

print(f"Ligando seleccionado: {ligand_sel}")
print(f"Átomos del ligando: {len(ligand)}\n")

# ===============================
# ALINEAR SOLO LA PROTEÍNA CONTRA complex.pdb
# ===============================
print("Alineando SOLO la proteína al complex.pdb...")

aligner = align.AlignTraj(
    u,
    ref,
    select=protein_sel,
    in_memory=True
)
aligner.run()

print("Alineación completada.\n")

# ===============================
# REFERENCIA DEL LIGANDO = FRAME 0 TRAS ALINEAR
# ===============================
u.trajectory[0]  # ir al frame 0 ya alineado
ref_positions = ligand.positions.copy()  # referencia interna

print("Referencia del ligando: posiciones en frame 0\n")

# ===============================
# CALCULAR RMSD DEL LIGANDO SIN FIT LOCAL
# ===============================
print("Calculando RMSD del ligando (sin fit, ref = frame 0)...")

rmsd_values = []

for ts in u.trajectory:
    diff = ligand.positions - ref_positions
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    rmsd_values.append(rmsd)

frames = np.arange(len(rmsd_values))
print(f"Frames analizados: {len(frames)}\n")

# ===============================
# GUARDAR SALIDA
# ===============================
df = pd.DataFrame({"Frame": frames, "RMSD_noFit_frame0": rmsd_values})
df.to_csv(output_csv, index=False)

print(f"✅ Datos guardados en: {output_csv}")

plt.figure(figsize=(7,4))
plt.plot(frames, rmsd_values, lw=1.2)
plt.xlabel("Frame")
plt.ylabel("RMSD Ligando (Å)")
plt.title("Ligand RMSD sin fit (referencia = frame 0)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(output_plot, dpi=300)

print(f"✅ Gráfico guardado como: {output_plot}")
print("\n=== Análisis completado ===")
