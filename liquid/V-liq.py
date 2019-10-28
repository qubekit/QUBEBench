from math import sqrt

import numpy as np

from simtk.openmm import app
from simtk.openmm import *
from simtk import unit


sites = 0
excep_pairs = []
normal_pairs = []
graph = {}


def find_shortest_path(graph, start, end, path=[]):
    path += [start]
    if start == end:
        return path
    if start not in graph.keys():
        return None
    shortest = None
    for node in graph[start]:
        if node not in path:
            newpath = find_shortest_path(graph, node, end, path)
            if newpath:
                if not shortest or len(newpath) < len(shortest):
                    shortest = newpath
    return shortest


def opls_lj(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce(
                 'epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)*4.0')
    lorentz.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    lorentz.setUseSwitchingFunction(True)
    lorentz.setSwitchingDistance(1.25 * unit.nanometer)
    lorentz.setUseLongRangeCorrection(True)
    # lorentz.setUseDispersionCorrection(True)
    system.addForce(lorentz)
    l_j_set = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        print(nonbonded_force.getParticleParameters(index))
        print(index)
        l_j_set[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        # if p1 != 11 and p2 != 11:
        lorentz.addExclusion(p1, p2)
        # print("exclusion %s %s"%(p1, p2))
        if eps._value != 0.0:
            # print p1,p2,sig,eps
            sig14 = sqrt(l_j_set[p1][0] * l_j_set[p2][0])
            eps14 = sqrt(l_j_set[p1][1] * l_j_set[p2][1])
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
        for x in range(len(excep_pairs)):
            for y in range(267):
                if p1 == excep_pairs[x, 0]+(y * (site_no + 1)) and p2 == excep_pairs[x, 1]+(y * (site_no + 1)) or p2 == excep_pairs[x, 0]+(y * (site_no + 1)) and p1 == excep_pairs[x, 1]+(y * (site_no + 1)):
                    charge1, sigma1, epsilon1 = nonbonded_force.getParticleParameters(p1)
                    charge2, sigma2, epsilon2 = nonbonded_force.getParticleParameters(p2)
                    q = charge1*charge2*0.5
                    # print(q)
                    sig14 = sqrt(sigma1 * sigma2)
                    eps = sqrt(epsilon1 * epsilon2)
                    nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
        for x in range(len(normal_pairs)):
            for y in range(267):
                if p1 == normal_pairs[x, 0]+(y*(site_no+1)) and p2 == normal_pairs[x, 1]+(y*(site_no+1)) or p2 == normal_pairs[x, 0]+(y*(site_no+1)) and p1 == normal_pairs[x, 1]+(y*(site_no+1)):
                    charge1, sigma1, epsilon1 = nonbonded_force.getParticleParameters(p1)
                    charge2, sigma2, epsilon2 = nonbonded_force.getParticleParameters(p2)
                    q = charge1*charge2
                    # print(q)
                    sig14 = sqrt(sigma1 * sigma2)
                    eps = sqrt(epsilon1 * epsilon2)
                    nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    # out=open('interactions.dat','w+')
    # for i in range(nonbonded_force.getNumExceptions()):
    # out.write(str(nonbonded_force.getExceptionParameters(i))+'\n')
    # out.close()
    return system


TEMP = 298.15 * unit.kelvin
pdb = app.PDBFile('NewBox_MOL.pdb')
modeller = app.Modeller(pdb.topology, pdb.positions)

forcefield = app.ForceField('MOL_extra.xml')
modeller.addExtraParticles(forcefield)
app.PDBFile.writeFile(modeller.topology, modeller.positions, open('MOL_modeller.pdb', 'w'))
 
system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=app.PME, ewaldErrorTolerance=0.0005, nonbondedCutoff=1.3 * unit.nanometer)

with open('gas_modeller.pdb', 'r') as pdb:

    for line in pdb:
        if 'CONECT' in line:
            graph[int(line.split()[1])-1] = [int(line.split()[2])-1]
            for i in range(3, len(line.split())):
                graph[int(line.split()[1])-1].append(int(line.split()[i])-1)

    print(graph)

with open('1-bromobutane_DDECV_BAD.z', 'r') as zmat:
    for line in zmat:
        if ' X0' in line:
            site_no = int(line.split()[0])-3
            graph[int(line.split()[0])-3] = [int(line.split()[4])-3]    # add site and parent
            graph[int(line.split()[4])-3].append(int(line.split()[0])-3)    # add site to parent list
            for i in range(site_no):

                k = find_shortest_path(graph, site_no, i, path=[])

                if len(k) == 4:
                    excep_pairs.append(i)
                    excep_pairs.append(site_no)
                elif len(k) >= 5:
                    normal_pairs.append(i)
                    normal_pairs.append(site_no)

excep_pairs = np.reshape(excep_pairs, (len(excep_pairs) // 2, 2))
normal_pairs = np.reshape(normal_pairs, (len(normal_pairs) // 2, 2))

system = opls_lj(system)
system.addForce(MonteCarloBarostat(1 * unit.bar, TEMP))
integrator = LangevinIntegrator(TEMP, 5 / unit.picosecond, 0.001 * unit.picoseconds)
simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
# simulation.context.computeVirtualSites()
print('MINIMIZATION STARTED')
simulation.minimizeEnergy()
# print('Energy at Minima is %3.3f kcal/mol' % (energy._value * KcalPerKJ))
print('MINIMIZATION DONE')
simulation.reporters.append(app.PDBReporter('output.pdb', 1000))
simulation.reporters.append(app.StateDataReporter('liquid.txt', 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.step(3000000)
np_equ_pos = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, np_equ_pos, open('NPT_EQ_FINAL.pdb', 'w'))
state = simulation.context.getState(getEnergy=True)
# print(state.getPotentialEnergy())
