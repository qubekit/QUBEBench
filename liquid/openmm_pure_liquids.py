import simtk.openmm as mm
from simtk.openmm import app
from simtk import unit


def opls_lj(system, switch_dist):

    forces = {system.getForce(indx).__class__.__name__: system.getForce(indx) for indx in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']

    # Custom combining rule
    lorentz = mm.CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')

    lorentz.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    lorentz.setUseSwitchingFunction(True)
    lorentz.setSwitchingDistance(switch_dist * unit.nanometer)
    lorentz.setUseLongRangeCorrection(True)
    system.addForce(lorentz)

    l_j_set = {}
    for indx in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(indx)
        l_j_set[indx] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(indx, charge, 0, 0)

    for indx in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(indx)
        # ALL THE 12, 13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value > 0:
            sig14 = (l_j_set[p1][0] * l_j_set[p2][0]) ** 0.5
            eps14 = (l_j_set[p1][1] * l_j_set[p2][1]) ** 0.5
            nonbonded_force.setExceptionParameters(indx, p1, p2, q, sig14, eps14)

    return system


def liquid_analysis(mol_name, switch_dist=1.25, temp=298.15):

    pdb = app.PDBFile('new.pdb')
    modeller = app.Modeller(pdb.topology, pdb.positions)
    forcefield = app.ForceField(mol_name + '.xml')
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, ewaldErrorTolerance=0.0005,
                                     nonbondedCutoff=(switch_dist + 0.05) * unit.nanometer)
    system = opls_lj(system, switch_dist)
    # FOR NPT
    temp *= unit.kelvin
    system.addForce(mm.MonteCarloBarostat(1 * unit.bar, temp))
    integrator = mm.LangevinIntegrator(temp, 5 / unit.picosecond, 0.001 * unit.picosecond)
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed', 'CudaDeviceIndex': '0'}

    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    print('MINIMISATION STARTED')
    simulation.minimizeEnergy()

    print('MINIMISATION DONE')
    simulation.reporters.append(app.DCDReporter('output.dcd', 1000))
    simulation.reporters.append(app.StateDataReporter(
        'liquid.txt', 1000, step=True, potentialEnergy=True, temperature=True, density=True))
    simulation.step(3000000)
    np_equ_pos = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, np_equ_pos, open('NPT_EQ_FINAL.pdb', 'w'))
