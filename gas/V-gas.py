from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

sites = 0
# parent_dist = 1
Excep_pairs = []
normal_pairs = []
graph = {}


def find_shortest_path(graph, start, end, path=[]):
    path = path + [start]
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


def OPLS_LJ(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in
              range(system.getNumForces())}
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
        print(nonbonded_force.getParticleParameters(index))
        LJset[index] = (sigma, epsilon)
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
            sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
        for x in range(len(Excep_pairs)):
            if p1 == Excep_pairs[x, 0] and p2 == Excep_pairs[x, 1] or p2 == Excep_pairs[x, 0] and p1 == Excep_pairs[x, 1]:
                charge1, sigma1, epsilon1 = nonbonded_force.getParticleParameters(p1)
                charge2, sigma2, epsilon2 = nonbonded_force.getParticleParameters(p2)
                q = charge1 * charge2 * 0.5
                # print(q)
                sig14 = sqrt(sigma1 * sigma2)
                eps = sqrt(epsilon1 * epsilon2)
                nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
        for x in range(len(normal_pairs)):
            if p1 == normal_pairs[x, 0] and p2 == normal_pairs[x, 1] or p2 == normal_pairs[x, 0] and p1 == normal_pairs[x, 1]:
                charge1, sigma1, epsilon1 = nonbonded_force.getParticleParameters(p1)
                charge2, sigma2, epsilon2 = nonbonded_force.getParticleParameters(p2)
                q = charge1 * charge2
                # print(q)
                sig14 = sqrt(sigma1 * sigma2)
                eps = sqrt(epsilon1 * epsilon2)
                nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    for i in range(nonbonded_force.getNumExceptions()):
        print(nonbonded_force.getExceptionParameters(i))

    return system


def minimise(simulation, iters=0):
    simulation.minimizeEnergy(maxIterations=iters)
    position = simulation.context.getState(getPositions=True).getPositions()
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    PDBFile.writeFile(simulation.topology, position, open('gasmin.pdb', 'w'))
    print('Energy at Minima is %3.3f kcal/mol' % (energy._value * KcalPerKJ))
    return simulation


temperature = 298.15 * kelvin
pdb = PDBFile('MET.pdb')
modeller = Modeller(pdb.topology, pdb.positions)

forcefield = ForceField('MET_extra.xml')
modeller.addExtraParticles(forcefield)
PDBFile.writeFile(modeller.topology, modeller.positions, open('gas_modeller.pdb', 'w'))

system = forcefield.createSystem(
    modeller.topology, nonbondedMethod=NoCutoff, constraints=None)
pdb = open('gas_modeller.pdb', 'r')
lines = pdb.readlines()
for line in lines:
    if 'CONECT' in line:
        graph[int(line.split()[1]) - 1] = [int(line.split()[2]) - 1]
        for i in range(3, len(line.split())):
            graph[int(line.split()[1]) - 1].append(int(line.split()[i]) - 1)
print(graph)
for line in lines:
    if 'HETATM' in line and ' X' in line:
        v_pos = np.array((float(line.split()[6]), float(line.split()[7]), float(line.split()[8])))
        sites = sites + 1
        site_no = int(line.split()[1]) - 1
        parent_dist = 1.1
        for line in lines:
            if 'HETATM' in line and 'C0' not in line and ' H0' not in line and ' X' not in line:
                parent_pos = np.array((float(line.split()[6]), float(line.split()[7]), float(line.split()[8])))
                dist = np.linalg.norm(parent_pos - v_pos)
                if dist <= parent_dist:
                    parent_dist = dist
                    parent = parent_pos
                    parent_no = int(line.split()[1]) - 1
        graph[site_no] = [parent_no]
        graph[parent_no].append(site_no)
        print(graph)
        # graph = {'A': ['B', 'C'],
        #     'B': ['C', 'D'],
        #     'C': ['D'],
        #     'D': ['C'],
        #     'E': ['F'],
        #     'F': ['C']}
        # print(graph)
        # find all of the paths between the extra site and the rest of the atoms in the molecule
        for i in range(site_no):

            k = find_shortest_path(graph, site_no, i, path=[])
            print(k)
            # print(len(k))
            if len(k) == 4:
                Excep_pairs.append(i)
                Excep_pairs.append(site_no)
            elif len(k) >= 5:
                normal_pairs.append(i)
                normal_pairs.append(site_no)
print(Excep_pairs)
Excep_pairs = np.reshape(Excep_pairs, (int(len(Excep_pairs) / 2), 2))
print(Excep_pairs)
print(normal_pairs)
normal_pairs = np.reshape(normal_pairs, (int(len(normal_pairs) / 2), 2))
print(normal_pairs)
system = OPLS_LJ(system)
integrator = LangevinIntegrator(
    temperature, 5 / picosecond, 0.0005 * picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation = minimise(simulation, 100)
simulation.reporters.append(app.PDBReporter('gas_output.pdb', 1000))
simulation.reporters.append(
    app.StateDataReporter('gas.txt', 1000, step=True, temperature=True, potentialEnergy=True, density=True,
                          totalSteps=10000, totalEnergy=True))
simulation.step(6000000)
# state = simulation.context.getState(getEnergy=True)
# print(state.getPotentialEnergy())
