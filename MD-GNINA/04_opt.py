# =============================
# IMPORTS
# =============================
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from openbabel import pybel
import os
import torch
import torchani
from ase import Atoms, io
from ase.optimize import BFGS
from torchani.units import HARTREE_TO_KCALMOL
from rdkit.Geometry import Point3D

workDir = os.getcwd()
Type = "smiles"   # smiles, pdb o mol
smiles_or_filename = "O=C1OC(CCCC(=O)CCCC=Cc2cc(O)cc(O)c12)C"

# =============================
# PARTE 1 — GENERAR ligand.mol Y ligand.pdb
# =============================

print("\n=== Generando estructura inicial del ligando ===")

if Type == "smiles":
    Smiles = smiles_or_filename
    smiles_fig = Chem.MolFromSmiles(Smiles)
    hmol = Chem.AddHs(smiles_fig)
    AllChem.EmbedMolecule(hmol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(hmol)
    AllChem.MolToMolFile(hmol, os.path.join(workDir, "ligand.mol"))
    AllChem.MolToPDBFile(hmol, os.path.join(workDir, "ligand.pdb"))

elif Type == "pdb":
    pdb_name = os.path.join(workDir, smiles_or_filename)
    mol = [m for m in pybel.readfile(filename=pdb_name, format='pdb')][0]
    out = pybel.Outputfile(filename="mol_tmp.mol", format="mol", overwrite=True)
    out.write(mol)
    out.close()

    mol_rd = Chem.MolFromMolFile("mol_tmp.mol")
    Smiles = Chem.MolToSmiles(mol_rd)

    smiles_fig = Chem.MolFromSmiles(Smiles)
    hmol = Chem.AddHs(smiles_fig)
    AllChem.EmbedMolecule(hmol)
    AllChem.MMFFOptimizeMolecule(hmol)
    AllChem.MolToMolFile(hmol, os.path.join(workDir, "ligand.mol"))
    AllChem.MolToPDBFile(hmol, os.path.join(workDir, "ligand.pdb"))

else:  # MOL input
    mol_name = os.path.join(workDir, smiles_or_filename)
    mol_rd = Chem.MolFromMolFile(mol_name)
    Smiles = Chem.MolToSmiles(mol_rd)

    smiles_fig = Chem.MolFromSmiles(Smiles)
    hmol = Chem.AddHs(smiles_fig)
    AllChem.EmbedMolecule(hmol)
    AllChem.MMFFOptimizeMolecule(hmol)
    AllChem.MolToMolFile(hmol, os.path.join(workDir, "ligand.mol"))
    AllChem.MolToPDBFile(hmol, os.path.join(workDir, "ligand.pdb"))

print("SMILES detectado:", Smiles)

# Imagen 2D
Draw.MolToFile(smiles_fig, os.path.join(workDir, Smiles + ".png"))
img = mpimg.imread(os.path.join(workDir, Smiles + ".png"))
plt.imshow(img)
plt.axis("off")
plt.show()

# =============================
# PARTE 2 — OPTIMIZACIÓN CON TORCHANI
# =============================

print("\n=== OPTIMIZACIÓN GEOMÉTRICA ANI-2x ===")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = torchani.models.ANI2x(periodic_table_index=True).to(device)
calculator = torchani.models.ANI2x().ase()

def mol2arr(mol, device=device):
    """Convierte RDKit Mol a tensores TorchANI"""
    pos = mol.GetConformer().GetPositions().tolist()
    atomnums = [a.GetAtomicNum() for a in mol.GetAtoms()]
    coordinates = torch.tensor([pos], requires_grad=True, device=device)
    species = torch.tensor([atomnums], device=device)
    return coordinates, species

# Leer el ligand.mol creado
mol_rd2 = Chem.MolFromMolFile(os.path.join(workDir, "ligand.mol"), removeHs=False)

# Cargar en ASE
atoms = io.read(os.path.join(workDir, "ligand.mol"))
atoms.center(vacuum=3.0)
atoms.set_calculator(calculator)

# Optimización
opt = BFGS(atoms)
opt.run(fmax=0.0001)

io.write(os.path.join(workDir, "ligand_min.xyz"), atoms)

# Actualizar coordenadas en RDKit
xyz = []
with open(os.path.join(workDir, "ligand_min.xyz")) as f:
    lines = f.readlines()[2:]
    for line in lines:
        parts = line.split()
        xyz.append([float(parts[1]), float(parts[2]), float(parts[3])])

conf = mol_rd2.GetConformer()
for i in range(mol_rd2.GetNumAtoms()):
    x, y, z = xyz[i]
    conf.SetAtomPosition(i, Point3D(x, y, z))

# Guardar formatos finales
AllChem.MolToMolFile(mol_rd2, os.path.join(workDir, "ligand_min.mol"))
AllChem.MolToPDBFile(mol_rd2, os.path.join(workDir, "ligand_min.pdb"))

# Convertir a pdbqt
os.system(
    "obabel -i mol ligand_min.mol -o pdbqt -O ligand_min.pdbqt -xh --partialcharge"
)

# Energía ANI
coords, species = mol2arr(mol_rd2, device)
energy = model((species, coords)).energies

print("Energy (hartree):", energy.item())
print("Energy (kcal/mol):", energy.item() * HARTREE_TO_KCALMOL)

print("\nArchivos generados:")
print(" - ligand.mol")
print(" - ligand.pdb")
print(" - ligand_min.xyz")
print(" - ligand_min.mol")
print(" - ligand_min.pdb")
print(" - ligand_min.pdbqt\n")
