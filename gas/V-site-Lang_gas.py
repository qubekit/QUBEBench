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
        'epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)*4.0')
    lorentz.setNonbondedMethod(nonbonded_force.getNonbondedMethod())
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
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


def Minimize(simulation,iters=0):
    simulation.minimizeEnergy(maxIterations=iters)
    position = simulation.context.getState(getPositions=True).getPositions()
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    PDBFile.writeFile(simulation.topology, position,
                          open('gasmin.pdb', 'w'))
    print ('Energy at Minima is %3.3f kcal/mol' % (energy._value * KcalPerKJ))
    return simulation

temperature = 298.15 * kelvin
pdb = PDBFile('PYR.pdb')

modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('PYR.xml')
modeller.addExtraParticles(forcefield)
PDBFile.writeFile(modeller.topology, modeller.positions, open('gas_modeller.pdb', 'w'))
system = forcefield.createSystem(
    modeller.topology, nonbondedMethod=NoCutoff,  constraints=None)
system = OPLS_LJ(system)
integrator = LangevinIntegrator(
    temperature, 5 / picosecond,  0.0005 * picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation = Minimize(simulation,100)
simulation.reporters.append(app.PDBReporter('gas_output.pdb', 1000))
simulation.reporters.append(app.StateDataReporter('gas.txt', 1000, step=True, temperature=True, potentialEnergy=True, density=True,totalSteps=10000, totalEnergy=True))
simulation.step(6000000)
