#!/usr/bin/env python3
"""
MD Production NPT (strides) con OpenMM y GAFF 2.11
- Carga el PDB equilibrado y archivo de estado (.rst XML).
- Mantiene periodicidad de la caja.
- Guarda DCD, PDB, estado (.rst) y checkpoint en cada stride.
- Uso recomendado: lanzar con 'nohup python 08_prod_md.py &'
"""

import os

os.environ["CUDA_VISIBLE_DEVICES"] = "1"

from sys import stdout
import openmm
from openmm import unit, XmlSerializer, MonteCarloBarostat, LangevinIntegrator
from openmm.app import PDBFile, DCDReporter, StateDataReporter, CheckpointReporter, Simulation, PME, HBonds, ForceField
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator

# -----------------------------
# Configuración de archivos
# -----------------------------
workDir = os.getcwd()
Equilibrated_PDB = 'prot_lig_equil.pdb'
Equilibrated_RST = 'prot_lig_equil.rst'  # XML de OpenMM
Ligand_file = 'ligand_prepared.sdf'
Jobname = 'prot_lig_prod'

# -----------------------------
# Parámetros de la simulación
# -----------------------------
Stride_Time_ns = 2       # tiempo de cada stride [ns]
Number_of_strides = 5     # cantidad de strides
Integration_timestep_fs = 2  # fs
Temperature_K = 300       # K
Pressure_bar = 1          # bar
Frames_per_stride = 100   # frames que queremos por stride
Checkpoint_ps = 100       # frecuencia de checkpoint en ps

# -----------------------------
# Convertir unidades OpenMM
# -----------------------------
stride_time_ps = Stride_Time_ns*1000
stride_time = stride_time_ps*unit.picoseconds
nstride = Number_of_strides
dt = Integration_timestep_fs*unit.femtoseconds
temperature = Temperature_K*unit.kelvin
pressure = Pressure_bar*unit.bar

# Pasos por stride
nsteps = int(stride_time.value_in_unit(unit.picoseconds)/dt.value_in_unit(unit.picoseconds))

# Guardar exactamente 100 frames por stride
nsavcrd = nsteps // Frames_per_stride
nchk = int(Checkpoint_ps*unit.picoseconds/dt)  # checkpoint en pasos
total_steps = nsteps * nstride

# -----------------------------
# Cargar sistema equilibrado
# -----------------------------
pdbfile = os.path.join(workDir, Equilibrated_PDB)
rstfile = os.path.join(workDir, Equilibrated_RST)

pdb = PDBFile(pdbfile)
topology = pdb.topology
positions = pdb.positions

# -----------------------------
# Cargar ligando GAFF 2.11
# -----------------------------
ligand = Molecule.from_file(os.path.join(workDir, Ligand_file))
ligand_generator = GAFFTemplateGenerator(molecules=[ligand])

# -----------------------------
# Crear forcefield
# -----------------------------
ff_protein = "amber19/protein.ff19SB.xml"
ff_water   = "amber19/tip3p.xml"
omm_forcefield = ForceField(ff_protein, ff_water)
omm_forcefield.registerTemplateGenerator(ligand_generator.generator)

# -----------------------------
# Crear sistema
# -----------------------------
system = omm_forcefield.createSystem(
    topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005
)
system.addForce(MonteCarloBarostat(pressure, temperature))

# -----------------------------
# Integrador
# -----------------------------
friction = 1.0/unit.picoseconds
integrator = LangevinIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(1e-6)

simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)

if pdb.topology.getPeriodicBoxVectors() is not None:
    simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# -----------------------------
# Loop de strides
# -----------------------------
for stride in range(1, nstride+1):
    print(f"\n>>> Stride #{stride} <<<")

    dcd_file = f"{Jobname}_{stride}.dcd"
    log_file = f"{Jobname}_{stride}.log"
    rst_file = f"{Jobname}_{stride}.rst"
    chk_file = f"{Jobname}_{stride}.chk"
    prev_rst = f"{Jobname}_{stride-1}.rst"

    # Cargar estado previo
    if stride == 1:
        print(f"> Cargando estado de equilibración: {rstfile}")
        with open(rstfile, 'r') as f:
            simulation.context.setState(XmlSerializer.deserialize(f.read()))
    else:
        print(f"> Cargando estado de stride anterior: {prev_rst}")
        with open(prev_rst, 'r') as f:
            simulation.context.setState(XmlSerializer.deserialize(f.read()))

    # Reporters
    simulation.reporters.append(DCDReporter(dcd_file, nsavcrd))
    simulation.reporters.append(CheckpointReporter(chk_file, nchk))
    simulation.reporters.append(StateDataReporter(
        stdout, nsavcrd,  # imprimir cada frame guardado
        step=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=total_steps,
        separator='\t\t',
        kineticEnergy=True,
        potentialEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True
    ))
    simulation.reporters.append(StateDataReporter(
        log_file, nsavcrd,
        step=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=total_steps,
        separator='\t\t',
        kineticEnergy=True,
        potentialEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True
    ))

    print(f"> Simulando {nsteps} pasos (Stride #{stride})...")
    simulation.step(nsteps)

    # Guardar estado final
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(rst_file, 'w') as f:
        f.write(XmlSerializer.serialize(state))

    # Guardar coordenadas finales
    positions = state.getPositions()
    pdb_out = f"{Jobname}_{stride}.pdb"
    PDBFile.writeFile(simulation.topology, positions, open(pdb_out, 'w'))

    simulation.reporters.clear()  # Limpiar reporters para el siguiente stride

print("\n>>> Simulación completada con éxito <<<")
