"""script to replace the swapping distances in the mdp file based on heavy atoms"""

from os import system
import sys


tag = sys.argv[1]
atoms = 0
light = 0

with open(f'{tag}.gro') as gro:
    lines = gro.readlines()

for line in lines:
    if f'1{tag}' in line:
        atoms += 1
        if 'H' in line.split()[1]:
            light += 1

heavy = atoms - light
print(f'Atoms = {atoms}')
print(f'light = {light}')
print(f'Heavy = {heavy}')

if heavy < 3:
    system("sed -i 's/rcoulomb                 = 1.0/rcoulomb                 = 1.1/g' GMXmdp.py")
    system("sed -i 's/rvdw-switch              = 0.95/rvdw-switch              = 1.05/g' GMXmdp.py")
    system("sed -i 's/rvdw                     = 1.0/rvdw                     = 1.1/g' GMXmdp.py")

elif 3 <= heavy <= 5:
    system("sed -i 's/rcoulomb                 = 1.0/rcoulomb                 = 1.3/g' GMXmdp.py")
    system("sed -i 's/rvdw-switch              = 0.95/rvdw-switch              = 1.25/g' GMXmdp.py")
    system("sed -i 's/rvdw                     = 1.0/rvdw                     = 1.3/g' GMXmdp.py")

elif heavy > 5:
    system("sed -i 's/rcoulomb                 = 1.0/rcoulomb                 = 1.5/g' GMXmdp.py")
    system("sed -i 's/rvdw-switch              = 0.95/rvdw-switch              = 1.45/g' GMXmdp.py")
    system("sed -i 's/rvdw                     = 1.0/rvdw                     = 1.5/g' GMXmdp.py")
