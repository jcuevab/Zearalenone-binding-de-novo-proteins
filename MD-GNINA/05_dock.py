import os, sys
import locale
import py3Dmol

# ============================
# 0. CONFIG GENERAL
# ============================

workDir = os.getcwd()
receptor = "receptor.pdb"
receptor_pdbqt = "receptor.pdbqt"
ligand = os.path.join(workDir, "ligand_min.pdbqt")   # ✅ tu ligando optimizado
res_box = "selection_merge.pdb"

# ✅ Coordenadas que ya me diste
centerX = -6.6712
centerY = -4.9850
centerZ = -5.2102
sizeX = 10
sizeY = 10
sizeZ = 10

# ============================
# 1. PARÁMETROS UNI DOCK
# ============================

seed = "0"
scoring = "vina"
search_mode = "detail"
exhaustiveness = 500
max_step = 40
num_modes = 10

locale.getpreferredencoding = lambda do_setlocale=True: "UTF-8"

# ============================
# 2. COMANDO UNIDOCK
# ============================
# Ruta real del ejecutable
unidock_exec = "/home/jvaldiviezo/bin/unidock"

unidock_cmd = (
    f"{unidock_exec}"
    f" --receptor {receptor_pdbqt} "
    f" --gpu_batch {ligand} "
    f" --search_mode {search_mode} "
    f" --scoring {scoring} "
    f" --center_x {centerX} --center_y {centerY} --center_z {centerZ} "
    f" --size_x {sizeX} --size_y {sizeY} --size_z {sizeZ} "
    f" --num_modes {num_modes} "
    f" --dir {workDir} "
    f" --seed {seed} "
    f" --max_step {max_step} "
    f" --exhaustiveness {exhaustiveness}"
)

# Guardar script
with open("unidock.sh", "w") as f:
    f.write(unidock_cmd + "\n")

os.system("chmod +x unidock.sh")
print("Ejecutando UniDock…")
os.system("./unidock.sh")

# ? Redirigir salida y error al archivo unidock.txt
print("Ejecutando UniDock")
os.system("./unidock.sh > unidock.txt 2>&1")

# ============================
# 3. CONVERSIÓN SALIDA
# ============================

ligand_out = os.path.join(workDir, "ligand_min_out.pdbqt")
ligand_out_sdf = os.path.join(workDir, "ligand_out.sdf")

os.system(f"obabel -i pdbqt {ligand_out} -o sdf -O {ligand_out_sdf} -xh")

# ============================
# 4. VISUALIZACIÓN (opcional)
# ============================

v = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js')
v.addModel(open(receptor).read())
v.setStyle({'cartoon': {}, 'stick': {'colorscheme':'white','radius':.1}})
v.addModel(open(res_box).read())
v.setStyle({'model':1},{'stick':{'colorscheme':'dimgrayCarbon','radius':.175}})
v.addModelsAsFrames(open(ligand_out_sdf,'rt').read())
v.setStyle({'model':2},{'stick':{'colorscheme':'greenCarbon'}})
v.zoomTo({'model':1})
v.show()
