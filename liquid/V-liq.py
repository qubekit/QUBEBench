from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
#from dcdreporter import DCDReporter
import numpy as np


sites=0
#parent_dist=1
Excep_pairs=[]
normal_pairs=[]
graph={}
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
    forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce(
                            'epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)*4.0')
    lorentz.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    lorentz.setUseSwitchingFunction(True)
    lorentz.setSwitchingDistance(1.25 * nanometer)
    lorentz.setUseLongRangeCorrection(True)
    #lorentz.setUseDispersionCorrection(True)
    system.addForce(lorentz)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        print(nonbonded_force.getParticleParameters(index))
        print(index)
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        #ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        #if p1 != 11 and p2 != 11:
        lorentz.addExclusion(p1, p2)
        #print("exclusion %s %s"%(p1, p2))
        if eps._value != 0.0:
        # print p1,p2,sig,eps
          sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
          eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
          nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
        for x in range(len(Excep_pairs)):
            for y in range(267):
              if p1 == Excep_pairs[x,0]+(y*(site_no+1)) and p2 == Excep_pairs[x,1]+(y*(site_no+1)) or p2 == Excep_pairs[x,0]+(y*(site_no+1)) and p1 == Excep_pairs[x,1]+(y*(site_no+1)):
                 charge1, sigma1, epsilon1 = nonbonded_force.getParticleParameters(p1)
                 charge2, sigma2, epsilon2 = nonbonded_force.getParticleParameters(p2)
                 q = charge1*charge2*0.5
                 #print(q)
                 sig14 = sqrt(sigma1 * sigma2)
                 eps = sqrt(epsilon1 * epsilon2)          
                 nonbonded_force.setExceptionParameters(i, p1, p2, q , sig14, eps)
        for x in range(len(normal_pairs)):
            for y in range(267):
              if p1 == normal_pairs[x,0]+(y*(site_no+1)) and p2 == normal_pairs[x,1]+(y*(site_no+1)) or p2 == normal_pairs[x,0]+(y*(site_no+1)) and p1 == normal_pairs[x,1]+(y*(site_no+1)):
                 charge1, sigma1, epsilon1 = nonbonded_force.getParticleParameters(p1)
                 charge2, sigma2, epsilon2 = nonbonded_force.getParticleParameters(p2)
                 q = charge1*charge2
                 #print(q)
                 sig14 = sqrt(sigma1 * sigma2)
                 eps = sqrt(epsilon1 * epsilon2)          
                 nonbonded_force.setExceptionParameters(i, p1, p2, q , sig14, eps)         
    #out=open('interactions.dat','w+')
    #for i in range(nonbonded_force.getNumExceptions()):
    #   out.write(str(nonbonded_force.getExceptionParameters(i))+'\n')
    #out.close()
    return system





TEMP = 298.15 * kelvin
pdb = PDBFile('NewBox_MOL.pdb')
modeller = Modeller(pdb.topology, pdb.positions)

forcefield = ForceField('MOL_extra.xml')
modeller.addExtraParticles(forcefield)
PDBFile.writeFile(modeller.topology, modeller.positions, open('MOL_modeller.pdb', 'w'))
 
system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=PME, ewaldErrorTolerance=0.0005, nonbondedCutoff=1.3 * nanometer)
pdb=open('gas_modeller.pdb','r')
lines=pdb.readlines()
for line in lines:
    if 'CONECT' in line:
        graph[int(line.split()[1])-1] = [int(line.split()[2])-1]
        for i in range(3,len(line.split())):
            graph[int(line.split()[1])-1].append(int(line.split()[i])-1)
print(graph)
pdb.close()
zmat=open('1-bromobutane_DDECV_BAD.z','r')
lines=zmat.readlines()
for line in lines:
    if ' X0' in line:
        site_no=int(line.split()[0])-3
        graph[int(line.split()[0])-3] = [int(line.split()[4])-3] #add site and parent 
        graph[int(line.split()[4])-3].append(int(line.split()[0])-3) #add site to parent list
        for i in range(site_no):
        
           k=find_shortest_path(graph, site_no, i, path=[])
           print(k)
           #print(len(k))
           if len(k) == 4:
              Excep_pairs.append(i)
              Excep_pairs.append(site_no)
           elif len(k) >= 5:
                normal_pairs.append(i)
                normal_pairs.append(site_no)
print(Excep_pairs)
Excep_pairs=np.reshape(Excep_pairs , (int(len(Excep_pairs)/2),2))
print(Excep_pairs)
print(normal_pairs) 
normal_pairs=np.reshape(normal_pairs , (int(len(normal_pairs)/2),2))   
print(normal_pairs)          
system = OPLS_LJ(system)
system.addForce(MonteCarloBarostat(1 * bar, TEMP))
integrator = LangevinIntegrator(TEMP, 5 / picosecond, 0.001 * picoseconds)
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
state = simulation.context.getState(getEnergy=True)
#print(state.getPotentialEnergy())

