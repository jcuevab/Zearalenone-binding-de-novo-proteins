#!/usr/bin/env python3
import os
import numpy as np
import parmed as pmd
from rdkit import Chem

# OpenMM imports
try:
    import openmm
    from openmm import app, unit, Vec3
except ImportError:
    from simtk import openmm, app, unit
    from simtk.openmm import Vec3

from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit.topology import Molecule
from openff.units.openmm import to_openmm

# -----------------------------
# CONFIGURACIÓN DE ARCHIVOS
# -----------------------------
workDir = os.getcwd()
ligand_out_sdf = os.path.join(workDir, "ligand_out.sdf")
receptor = os.path.join(workDir, "receptor.pdb")

# -----------------------------
# PARÁMETROS GENERALES
# -----------------------------
Force_field = "AMBER19SB"
Water_type = "TIP3P"
Padding_distance = 4.0       # angstroms
Ions = "NaCl"
Concentration = 0.15         # M
pH = 7.4

# -----------------------------
# LEER LIGANDO
# -----------------------------
results = Chem.SDMolSupplier(ligand_out_sdf, removeHs=False)
ligand_mol = results[0]
Chem.MolToMolFile(ligand_mol, os.path.join(workDir, "ligand_prepared.sdf"))

ligand = Molecule.from_file(os.path.join(workDir, "ligand_prepared.sdf"))
ligand_positions = ligand.conformers[0]
ligand_topology = ligand.to_topology()

# -----------------------------
# CARGAR FORCE FIELDS
# -----------------------------
ff_protein = "amber19/protein.ff19SB.xml"
ff_water = "amber19/tip3p.xml"
omm_forcefield = app.ForceField(ff_protein, ff_water)

# Ligando con GAFF 2.11
ligand_generator = GAFFTemplateGenerator(molecules=[ligand])
omm_forcefield.registerTemplateGenerator(ligand_generator.generator)

print("Archivos XML de force field cargados:")

# -----------------------------
# CARGAR PDB DE RECEPTOR
# -----------------------------
pdb = app.PDBFile(receptor)
modeller = app.Modeller(pdb.topology, pdb.positions)

# -----------------------------
# AÑADIR LIGANDO
# -----------------------------
modeller.add(ligand_topology.to_openmm(), to_openmm(ligand_positions))

# -----------------------------
# HIDROGENAR SEGÚN PH
# -----------------------------
modeller.addHydrogens(omm_forcefield, pH=pH)

# -----------------------------
# SOLVATACIÓN E IONES
# -----------------------------
positive_ion = 'Na+' if Ions == 'NaCl' else 'K+'

modeller.addSolvent(
    omm_forcefield,
    model='tip3p',
    padding=Padding_distance*unit.angstrom,
    ionicStrength=float(Concentration)*unit.molar,
    positiveIon=positive_ion,
    negativeIon='Cl-'
)

# -----------------------------
# OBTENER TOPOLOGÍA Y POSICIONES
# -----------------------------
topology = modeller.getTopology()
positions = modeller.getPositions()

# -----------------------------
# DEFINIR CAJA PERIÓDICA (Vec3)
# -----------------------------
pos_nm = np.array(positions.value_in_unit(unit.nanometer))
min_pos = pos_nm.min(axis=0)
max_pos = pos_nm.max(axis=0)
vec = max_pos - min_pos + Padding_distance/10.0  # padding extra en nm

box_vectors = (Vec3(vec[0], 0, 0),
               Vec3(0, vec[1], 0),
               Vec3(0, 0, vec[2]))

modeller.topology.setPeriodicBoxVectors(box_vectors)
print(f"Caja periódica (nm): {vec}")

# -----------------------------
# EXPORTAR PDB
# -----------------------------
system_pdb = os.path.join(workDir, "system.pdb")
app.PDBFile.writeFile(topology, positions, open(system_pdb, 'w'))
print(f"Archivo PDB generado: {system_pdb}")

# -----------------------------
# CREAR SISTEMA PARA PARMED
# -----------------------------
system = omm_forcefield.createSystem(modeller.topology,
                                     nonbondedMethod=app.PME,
                                     nonbondedCutoff=1.0*unit.nanometer,
                                     constraints=app.HBonds,
                                     flexibleConstraints=True)
pmd_struct = pmd.openmm.load_topology(topology, system, positions)

# -----------------------------
# EXPORTAR AMBER
# -----------------------------
amber_prmtop = os.path.join(workDir, "system_amber.prmtop")
amber_inpcrd = os.path.join(workDir, "system_amber.inpcrd")
pmd_struct.save(amber_prmtop, overwrite=True)
pmd_struct.save(amber_inpcrd, overwrite=True)

print("Archivos AMBER generados: ", amber_prmtop, amber_inpcrd)
print("¡Topología periódica generada con éxito!")
