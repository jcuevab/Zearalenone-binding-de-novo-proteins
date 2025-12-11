import os
import sys
import csv
import subprocess

# === Ajustar rutas locales ===
workDir = os.getcwd()
receptor = os.path.join(workDir, "receptor.pdb")

# Ruta real a P2Rank
prank_bin = "/home/jvaldiviezo/bin/p2rank_2.4/prank"

# Carpeta de salida
output_p2rank = os.path.join(workDir, "output_p2rank")

# Crear comando exacto
p2rank_cmd = f"{prank_bin} predict -f {receptor} -o {output_p2rank}"

# Guardar como script ejecutable (igual que en Colab)
script_path = os.path.join(workDir, "p2rank.sh")
with open(script_path, "w") as f:
    f.write(p2rank_cmd + "\n")

# Permisos
subprocess.run(["chmod", "700", script_path])

# Ejecutar P2Rank usando bash
subprocess.run(["bash", script_path])

# === Procesar el CSV ===
csv_file = os.path.join(output_p2rank, "receptor.pdb_predictions.csv")

with open(csv_file, "r") as f:
    csvreader = csv.reader(f)
    rows = list(csvreader)

# Columnas
residue = []
score = []
cx = []
cy = []
cz = []

for row in rows:
    residue.append(row[9:10])
    score.append(row[2:3])
    cx.append(row[6:7])
    cy.append(row[7:8])
    cz.append(row[8:9])

# Imprimir resultados de cada bolsillo
for i in range(1, len(residue)):
    file = str((residue[i])[0]).split()
    score_end = str((score[i])[0]).split()
    center_x_end = str((cx[i])[0]).split()
    center_y_end = str((cy[i])[0]).split()
    center_z_end = str((cz[i])[0]).split()

    print(f"Pocket {i}")
    print("Score =", score_end[0])

    final_res = []
    for r in file:
        final_res.append(int(r[2:]))  # Quita “A:” o “B:” etc.

    print("Selected Residues =", final_res)
    print(f"Center x = {center_x_end[0]}  y = {center_y_end[0]}  z = {center_z_end[0]}\n")
