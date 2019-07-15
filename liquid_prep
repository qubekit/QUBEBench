import networkx as nx
import os


def create_box(pdb_name, nmol=267):

    x = y = z = 4       # Liquid box dimensions
    resn = 'MOL'        # Enter new molecule description code to rename UNK

    with open(pdb_name, 'rt') as f_in, open(f'{resn}.pdb', 'wt') as f_out:
        for line in f_in:
            f_out.write(line.replace('UNK', resn))

    os.system(f'gmx insert-molecules -ci  {resn}.pdb -box {x} {y} {z}  -nmol {nmol} -o box.pdb')


def generate_connections(single_pdb, all_pdb, n_molecules):

    os.system(f'head -n -2 {all_pdb} > temp.pdb ; mv temp.pdb new.pdb')

    with open(single_pdb, 'r') as s_pdb, open('new.pdb', 'a+') as new_pdb:

        lines = s_pdb.readlines()

        topology = nx.Graph()

        for line in lines:
            if 'CONECT' in line:
                topology.add_node(int(line.split()[1]))
                for i in range(2, len(line.split())):
                    if int(line.split()[i]) != 0:
                        topology.add_node(int(line.split()[i]))
                        topology.add_edge(int(line.split()[1]), int(line.split()[i]))

        n_atoms = len(list(topology.nodes))

        for mol in range(n_molecules):
            for node in topology.nodes:
                bonded = sorted(list(nx.neighbors(topology, node)))
                if len(bonded) > 1:
                    new_pdb.write(f'CONECT{node + (mol * n_atoms):5}{"".join(f"{x + (mol * n_atoms):5}" for x in bonded)}\n')

        new_pdb.write('END\n')
