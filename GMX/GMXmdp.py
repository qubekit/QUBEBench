def grace_fep_job(window,np=1):
    ofile=open('fep_window_%d.sh' % window, 'w+')
    jobscript='''#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=GMX
#SBATCH --ntasks=%d --nodes=%d 
#SBATCH --mem-per-cpu=6000 
#SBATCH --time=24:00:00
module load Apps/Gromacs/5.0.5

#GMX FEP JOB SCRIPT CREATED BY LSD
echo 'STARTING MINIMIZATION'
mpirun -np %d gmx_mpi grompp -f em_steep_0.mdp -c solvated_BNZ.gro -p topol.top -o em_steep_0.tpr 
mpirun -np %d gmx_mpi mdrun -v -deffnm em_steep_0
echo 'STARTING NVT EQUILIBRATION'
mpirun -np %d gmx_mpi grompp -f nvt_equil_0.mdp -c em_steep_0.gro -p topol.top -o  nvt_equil_0.tpr
mpirun -np %d gmx_mpi mdrun -v -deffnm nvt_equil_0
echo 'STARTING NPT EQUILIBRATION w Velocoties '
mpirun -np %d gmx_mpi grompp -f npt_equil_0.mdp -c nvt_equil_0.gro -t nvt_equil_0.cpt -p topol.top -o npt_equil_0.tpr
mpirun -np %d gmx_mpi mdrun -v -deffnm npt_equil_0
echo 'STARTING NPT PRODUCTION'
mpirun -np %d gmx_mpi grompp -f npt_prod_0.mdp -c npt_equil_0.gro -t npt_equil_0.cpt -p topol.top -o npt_prod_0.tpr
mpirun -np %d gmx_mpi mdrun -v -deffnm npt_prod_0
echo 'DONE WITH THE JOB'
echo `date` 
'''%(np,np,np,np,np,np,np,np,np)
    ofile.write(jobscript.replace('_0', '_%d' % window))
    ofile.close()
    return None


def write_fep_job(window):
    ofile=open('fep_window_%d.sh'%window,'w+') 
    jobscript='''#GMX FEP JOB SCRIPT CREATED BY LSD
echo 'STARTING MINIMIZATION'
gmx grompp -f em_steep_0.mdp -c solvated_UNK.gro -p topol.top -o em_steep_0.tpr -maxwarn 2
gmx mdrun -nt 4 -v -deffnm em_steep_0
echo 'STARTING NVT EQUILIBRATION'
gmx grompp -f nvt_equil_0.mdp -c em_steep_0.gro -p topol.top -o  nvt_equil_0.tpr -maxwarn 2
gmx mdrun -nt 4  -v -deffnm nvt_equil_0
echo 'STARTING NPT EQUILIBRATION w Velocoties '
gmx grompp -f npt_equil_0.mdp -c nvt_equil_0.gro -t nvt_equil_0.cpt -p topol.top -o npt_equil_0.tpr -maxwarn 2
gmx mdrun -nt 4  -v -deffnm npt_equil_0
echo 'STARTING NPT PRODUCTION'
gmx grompp -f npt_prod_0.mdp -c npt_equil_0.gro -t npt_equil_0.cpt -p topol.top -o npt_prod_0.tpr -maxwarn 2
gmx mdrun -nt 4  -v -deffnm npt_prod_0
echo 'DONE WITH THE JOB'
echo `date` 
'''
    ofile.write(jobscript.replace('_0', '_%d' % window))
    ofile.close()
    return None


def create_top(resname):
    pgn = open('topol.top', 'w+')
    pgn.write('''#include "oplsaa.ff/forcefield.itp"
#include "%s.itp"
#include "oplsaa.ff/tip4pew.itp"

[ system ]
; Name
%s in water

[ molecules ]
; Compound             #mols
%s                1
''' % (resname, resname, resname)
    )
    pgn.close()
    return None


def setup_box(resname):
    from distutils import spawn 
    try: 
        GMX = spawn.find_executable('gmx')
    except:
        print('GROMACS executable not found')
    solvate_commands = [
        '%s editconf -f %s.gro -o box_%s.gro -c -d 1.5 -bt cubic' % (
            GMX, resname, resname),
        '%s solvate -cp box_%s.gro -cs tip4p -o solvated_%s.gro -p topol.top' % (
            GMX, resname, resname),
    ]
    import os
    try:
        os.system(solvate_commands[0])
        os.system(solvate_commands[1])
        res = True
    except:
        res = False
        print('BOX setup failed')
    return None


def minimize_lbfgs(resname, window=None, anh=None):
    if window is not None:
        ofile = open('em_lbfgs_%d.mdp' % (window), 'w+')
    else:
        ofile = open('em_lbfgs.mdp', 'w+')
    ofile.write(''';Run control
integrator               = l-bfgs
nsteps                   = 5000
define                   = -DFLEXIBLE
; EM criteria and other stuff
emtol                    = 100
emstep                   = 0.01
niter                    = 20
nbfgscorr                = 10
; Output control
nstlog                   = 1
nstenergy                = 1
; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 1
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdwtype                  = cutoff
vdw-modifier             = potential-switch 
rvdw-switch              = 0.95
rvdw                     = 1.0
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
; Temperature and pressure coupling are off during EM
tcoupl                   = no
pcoupl                   = no
''')
    if anh == 'Q':
        ofile.write(Q_FEP(resname, window, init='vdw-q', final='vdw'))
    elif anh == 'Vdw':
        ofile.write(Vdw_FEP(resname, window, init='vdw', final='none'))
    ofile.close()
    return None


def minimize_steep(resname, window=None, anh=None):
    if window is not None:
        ofile = open('em_steep_%d.mdp' % (window), 'w+')
    else:
        ofile = open('em_steep.mdp', 'w+')
    ofile.write('''
; Run control
integrator               = steep
nsteps                   = 5000
; EM criteria and other stuff
emtol                    = 100
emstep                   = 0.01
niter                    = 20
nbfgscorr                = 10
; Output control
nstlog                   = 1
nstenergy                = 1
; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 1
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.2
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdwtype                  = cutoff
vdw-modifier             = potential-switch
rvdw-switch              = 0.95
rvdw                     = 1.0
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
; Temperature and pressure coupling are off during EM
tcoupl                   = no
pcoupl                   = no
; No velocities during EM
gen_vel                  = no
; options for bonds
constraints              = none  ; we only have C-H bonds here
''')
    if anh == 'Q':
        ofile.write(Q_FEP(resname, window, init='vdw-q', final='vdw'))
    elif anh == 'Vdw':
        ofile.write(Vdw_FEP(resname, window, init='vdw', final='none'))
    ofile.close()
    return None


def Vdw_FEP(resname, window=0, init='vdw-q', final='none'):
    fep_control = ('''; Free energy control stuff
free_energy              = yes
init_lambda_state        = %d
delta_lambda             = 0
calc_lambda_neighbors    = 1        ; only immediate neighboring windows
; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15  
vdw_lambdas              = 0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.55 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
coul_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
bonded_lambdas           = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
restraint_lambdas        = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
mass_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
temperature_lambdas      = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no       ; linear interpolation of Coulomb (none in this case)
sc-power                 = 1
sc-sigma                 = 0.3
couple-moltype           = %s  ; name of moleculetype to decouple
couple-lambda0           = %s      ; only van der Waals interactions
couple-lambda1           = %s     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 100
''' % (window, resname, init, final))
    return fep_control


def Q_FEP(resname, window=0, init='vdw-q', final='none'):
    fep_control = ('''; Free energy control stuff
free_energy              = yes
init_lambda_state        = %d
delta_lambda             = 0
calc_lambda_neighbors    = 1        ; only immediate neighboring windows
; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4   
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00
coul_lambdas             = 0.00 0.25 0.50 0.75 1.00
bonded_lambdas           = 0.00 0.00 0.00 0.00 0.00
restraint_lambdas        = 0.00 0.00 0.00 0.00 0.00
mass_lambdas             = 0.00 0.00 0.00 0.00 0.00
temperature_lambdas      = 0.00 0.00 0.00 0.00 0.00
; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no       ; linear interpolation of Coulomb (none in this case)
sc-power                 = 1
sc-sigma                 = 0.3
couple-moltype           = %s  ; name of moleculetype to decouple
couple-lambda0           = %s      ; only van der Waals interactions
couple-lambda1           = %s     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 100
''' % (window, resname, init, final))
    return fep_control


def NVT_Equilibrate(resname, window=None, anh=None,solname='SOL'):
    if window is not None:
        ofile = open('nvt_equil_%d.mdp' % (window), 'w+')
    else:
        ofile = open('nvt_equil.mdp', 'w+')
    nvt_equil = '''; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 100000    ; 100 ps
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxout-compressed       = 0
; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdwtype                  = cutoff
vdw-modifier             = potential-switch
rvdw-switch              = 0.95
rvdw                     = 1.0
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tcoupl      = V-rescale                     ; modified Berendsen thermostat
tc-grps     = non-Water Water    ; two coupling groups - more accurate
tau_t       = 0.1   0.1                     ; time constant, in ps
ref_t       = 300   300                     ; reference temperature, one for each group, in K
; Pressure coupling is off for NVT
Pcoupl                   = No
tau_p                    = 0.5
compressibility          = 4.5e-05
ref_p                    = 1.0
; Generate velocities to start
gen_vel                  = yes
gen_temp                 = 300
gen_seed                 = -1
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Do not constrain the starting configuration
continuation             = no
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
'''
    ofile.write(nvt_equil)
    if anh == 'Q':
        ofile.write(Q_FEP(resname, window, init='vdw-q', final='vdw'))
    elif anh == 'Vdw':
        ofile.write(Vdw_FEP(resname, window, init='vdw', final='none'))
    ofile.close()
    return None


def NPT_Equilibrate(resname, window=None, anh=None,solname='SOL'):
    if window is not None:
        ofile = open('npt_equil_%d.mdp' % (window), 'w+')
    else:
        ofile = open('npt_equil.mdp', 'w+')
    npt_equil = '''; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 100000    ; 200 ps
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxout-compressed       = 0
; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdwtype                  = cutoff
vdw-modifier             = potential-switch
rvdw-switch              = 0.95
rvdw                     = 1.0
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
; Temperature coupling
tcoupl      = V-rescale                     ; modified Berendsen thermostat
tc-grps     = non-Water Water    ; two coupling groups - more accurate
tau_t       = 0.1   0.1                     ; time constant, in ps
ref_t       = 300   300                     ; reference temperature, one for each group, in K
; Pressure coupling is on for NPT
Pcoupl                   = Parrinello-Rahman 
tau_p                    = 1.0
compressibility          = 4.5e-05
ref_p                    = 1.0 
; Do not generate velocities
gen_vel                  = no 
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Constrain the starting configuration
; since we are continuing from NVT
continuation             = yes 
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
'''
    ofile.write(npt_equil)
    if anh == 'Q':
        ofile.write(Q_FEP(resname, window, init='vdw-q', final='vdw'))
    elif anh == 'Vdw':
        ofile.write(Vdw_FEP(resname, window, init='vdw', final='none'))
    ofile.close()
    return None


def NPT_Production(resname, window=None, anh=None,solname='SOL'):
    if window is not None:
        ofile = open('npt_prod_%d.mdp' % (window), 'w+')
    else:
        ofile = open('npt_prod.mdp', 'w+')
    npt_equil = '''; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 1500000   ; 1 ns
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxout-compressed       = 0
; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdwtype                  = cutoff
vdw-modifier             = potential-switch
rvdw-switch              = 0.95
rvdw                     = 1.0
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
; Temperature coupling
tcoupl      = V-rescale                     ; modified Berendsen thermostat
tc-grps     = non-Water Water    ; two coupling groups - more accurate
tau_t       = 0.1   0.1                     ; time constant, in ps
ref_t       = 300   300                     ; reference temperature, one for each group, in K
; Pressure coupling is on for NPT
Pcoupl                   = Parrinello-Rahman 
tau_p                    = 1.0
compressibility          = 4.5e-05
ref_p                    = 1.0 
; Do not generate velocities
gen_vel                  = no 
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Constrain the starting configuration
; since we are continuing from NVT
continuation             = yes 
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
'''
    ofile.write(npt_equil)
    if anh == 'Q':
        ofile.write(Q_FEP(resname, window, init='vdw-q', final='vdw'))
    elif anh == 'Vdw':
        ofile.write(Vdw_FEP(resname, window, init='vdw', final='none'))
    ofile.close()
    return None

