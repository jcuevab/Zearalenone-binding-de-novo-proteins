#!/usr/bin/env python3
import os

# *** FORZAR USO DE GPU 1 ***
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

from sys import stdout
import parmed as pmd
from openmm.app import PDBFile, AmberPrmtopFile, AmberInpcrdFile
from openmm.app import Simulation, DCDReporter, StateDataReporter, PME, HBonds
from openmm import LangevinIntegrator, unit, CustomExternalForce, XmlSerializer, MonteCarloBarostat


# -----------------------------
# CONFIGURACIÓN DE ARCHIVOS
# -----------------------------
workDir = os.getcwd()
prmtop_file = os.path.join(workDir, "system_amber.prmtop")
inpcrd_file = os.path.join(workDir, "system_amber.inpcrd")
pdb_file_out = os.path.join(workDir, "prot_lig_equil.pdb")
dcd_file = os.path.join(workDir, "prot_lig_equil.dcd")
log_file = os.path.join(workDir, "prot_lig_equil.csv")
rst_file = os.path.join(workDir, "prot_lig_equil.rst")

# -----------------------------
# PARÁMETROS DE SIMULACIÓN
# -----------------------------
temperature = 300.0 * unit.kelvin
pressure = 1.0 * unit.bar
friction = 1.0 / unit.picosecond
timestep = 2.0 * unit.femtoseconds
nsteps_minimization = 20000
equilibration_ns = 0.5
save_stride_ps = 10.0
restraint_fc = 50.0  # kJ/mol, poner >0 para activar restraints

# -----------------------------
# CARGAR TOPOLOGÍA AMBER
# -----------------------------
prmtop = AmberPrmtopFile(prmtop_file)
inpcrd = AmberInpcrdFile(inpcrd_file)

# Crear sistema
system = prmtop.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=HBonds,
    rigidWater=True
)

# Agregar barostato NPT
system.addForce(MonteCarloBarostat(pressure, temperature))

# Integrador
integrator = LangevinIntegrator(temperature, friction, timestep)

# -----------------------------
# CREAR SIMULACIÓN
# -----------------------------
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# -----------------------------
# APLICAR RESTRICCIONES (OPCIONAL)
# -----------------------------
if restraint_fc > 0:
    restraint = CustomExternalForce('k*periodicdistance(x,y,z,x0,y0,z0)^2')
    restraint.addPerParticleParameter('k')
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')
    for i, atom in enumerate(prmtop.topology.atoms()):
        if atom.element.symbol != 'H':
            pos = inpcrd.positions[i]
            restraint.addParticle(i, [restraint_fc, pos.x, pos.y, pos.z])
    system.addForce(restraint)

# -----------------------------
# MINIMIZACIÓN DE ENERGÍA
# -----------------------------
print(">> Minimización de energía ...")
simulation.minimizeEnergy(maxIterations=nsteps_minimization)
state = simulation.context.getState(getEnergy=True)
print("Energía potencial mínima:", state.getPotentialEnergy())

# Configurar velocidades
simulation.context.setVelocitiesToTemperature(temperature)

# -----------------------------
# REPORTERS
# -----------------------------
stride_steps = int((save_stride_ps * unit.picoseconds) / timestep)
simulation.reporters.append(DCDReporter(dcd_file, stride_steps))
simulation.reporters.append(StateDataReporter(
    log_file, stride_steps,
    step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    temperature=True, density=True, separator=","
))
simulation.reporters.append(StateDataReporter(
    stdout, stride_steps,
    step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    temperature=True, density=True, separator=","
))

# -----------------------------
# EQUILIBRACIÓN
# -----------------------------
nsteps_equil = int((equilibration_ns * 1000 * unit.picoseconds) / timestep)
print(f">> Equilibrando {equilibration_ns} ns = {nsteps_equil} steps ...")
simulation.step(nsteps_equil)

# -----------------------------
# GUARDAR ARCHIVOS FINALES
# -----------------------------
# Último frame PDB
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open(pdb_file_out, 'w'))

# Último estado XML como .rst
state = simulation.context.getState(getPositions=True, getVelocities=True)
with open(rst_file, 'w') as f:
    f.write(XmlSerializer.serialize(state))

print(">> Equilibración completada")
print(f">> Archivos generados:\n- DCD: {dcd_file}\n- PDB: {pdb_file_out}\n- CSV log: {log_file}\n- RST: {rst_file}")
