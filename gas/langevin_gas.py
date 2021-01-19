import simtk.openmm as mm
from simtk.openmm import app
from simtk import unit


def opls_lj(system):

    forces = {system.getForce(indx).__class__.__name__: system.getForce(indx) for indx in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = mm.CustomNonbondedForce(
        'epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)*4.0')
    lorentz.setNonbondedMethod(nonbonded_force.getNonbondedMethod())
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    system.addForce(lorentz)
    l_j_set = {}

    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        l_j_set[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(index, charge, sigma, epsilon * 0)

    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            sig14 = (l_j_set[p1][0] * l_j_set[p2][0]) ** 0.5
            eps14 = (l_j_set[p1][1] * l_j_set[p2][1]) ** 0.5
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps14)

    return system


def minimise(simulation, iters=0):

    simulation.minimizeEnergy(maxIterations=iters)
    position = simulation.context.getState(getPositions=True).getPositions()
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    app.PDBFile.writeFile(simulation.topology, position, open('gasmin.pdb', 'w'))

    print('Energy at minima is %3.3f kcal/mol' % (energy._value * mm.KcalPerKJ))

    return simulation


def gas_analysis(mol_name, temp=298.15):

    temp *= unit.kelvin
    pdb = app.PDBFile(mol_name + '.pdb')

    modeller = app.Modeller(pdb.topology, pdb.positions)
    forcefield = app.ForceField(mol_name + '.xml')

    try:
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None
        )
    except ValueError:
        modeller.addExtraParticles(forcefield)
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None
        )

    system = opls_lj(system)
    integrator = mm.LangevinIntegrator(temp, 5 / unit.picosecond,  0.0005 * unit.picosecond)
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    print(f'Results for: {mol_name}.')
    simulation = minimise(simulation, 100)
    simulation.reporters.append(app.DCDReporter('gas_output.dcd', 1000))
    simulation.reporters.append(app.StateDataReporter(
        'gas.txt', 1000, step=True, temperature=True, potentialEnergy=True, density=True, totalSteps=10000,
        totalEnergy=True))
    simulation.step(6000000)
