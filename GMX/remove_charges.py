"""script to remove charges from itp file"""

import sys


tag = sys.argv[1]

with open(tag + '.itp', 'r') as file:
    lines = file.readlines()

for n, line in enumerate(lines):
    if " atoms " in line:
        atoms = n
    elif " bonds " in line:
        bonds = n

with open(tag + '.itp', 'r') as file, open('new.itp', 'w+') as out:
    lines = file.readlines()

    for n, line in enumerate(lines):
        if bonds > n > atoms + 1:
            out.write('%6i   %8s      1    %s   %s      1%11.4f%11.4f\n'%(int(line.split()[0]), str(line.split()[1]), str(line.split()[3]), str(line.split()[4]), 0, float(line.split()[7])))
        else:
            out.write(line)
