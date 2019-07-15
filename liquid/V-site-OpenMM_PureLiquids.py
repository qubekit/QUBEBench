from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
#from dcdreporter import DCDReporter
import numpy as np


def OPLS_LJ(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')
    lorentz.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    lorentz.setUseSwitchingFunction(True)
    lorentz.setSwitchingDistance(1.45 * nanometer)
    lorentz.setUseLongRangeCorrection(True)
    #lorentz.setUseDispersionCorrection(True)
    system.addForce(lorentz)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            # print p1,p2,sig,eps
            sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    return system

pdb = PDBFile('NewBox_PYR.pdb')
modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('PYR.xml')
modeller.addExtraParticles(forcefield)
PDBFile.writeFile(modeller.topology, modeller.positions, open('PYR_modeller.pdb', 'w'))
system = forcefield.createSystem(
    modeller.topology, nonbondedMethod=PME, ewaldErrorTolerance=0.0005, nonbondedCutoff=1.5 * nanometer)
system = OPLS_LJ(system)
# FOR NPT
TEMP = 298.15 * kelvin
system.addForce(MonteCarloBarostat(1 * bar, TEMP))
integrator = LangevinIntegrator(TEMP, 1 / picosecond, 0.001 * picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
#simulation.context.computeVirtualSites()
print('MINIMIZATION STARTED')
simulation.minimizeEnergy()
#print('Energy at Minima is %3.3f kcal/mol' % (energy._value * KcalPerKJ))
print('MINIMIZATION DONE')
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter('liquid.txt', 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.step(3000000)
np_equ_pos = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, np_equ_pos, open('NPT_EQ_FINAL.pdb', 'w'))
