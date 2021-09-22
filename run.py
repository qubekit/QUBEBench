"""
Should be run from either inside an 10_finalise/11_finalise folder (single execution)
or from inside a folder containing many QUBEKit_name folders (bulk execution).
"""

from gas.langevin_gas import gas_analysis
from liquid.openmm_pure_liquids import liquid_analysis

from liquid_prep import create_box, generate_connections

from forcebalance import write_fb_csv
from properties import calculate_properties
from utils.decorators import exception_catcher
from utils.helpers import generate_bulk_csv, make_and_change_into, mol_data_from_csv

import argparse
import csv
import os
import sys


class Molecule:
    """
    Just a data storing class, similar to to Ligand() in QUBEKit.
    Could probably be a DataClass or Struct but w/e.
    """

    def __init__(self):

        self.name = None
        self.home = None
        self.temp = 298.15
        self.nmol = 500
        self.switch_dist = None
        self.analyse = 'both'
        self.csv_path = None
        self.bulk_run = False
        self.forcebalance = False


class ArgsAndConfigs:
    """
    Parses any args and configs and handles bulk when necessary
    """

    def __init__(self):

        self.args = self.parse_commands()

        if self.args.bulk_run:
            self.handle_bulk()

        else:
            self.molecule = Molecule()

            # Set any arguments which have been changed via the terminal commands
            for key, val in vars(self.args).items():
                if val is not None:
                    setattr(self.molecule, key, val)

            self.molecule.home = os.getcwd()
            self.molecule.csv_path = os.path.join(os.getcwd(), 'results.csv')

            self.start_results_file(self.molecule.csv_path)

            Execute(self.molecule)

    @staticmethod
    def parse_commands():

        def string_to_bool(string):
            """Convert a string to a bool for argparse use when casting to bool"""
            return string.casefold() in ['true', 't', 'yes', 'y']

        class CSVAction(argparse.Action):
            """The csv creation class run when the csv option is used."""

            def __call__(self, pars, namespace, values, option_string=None):
                """This function is executed when csv is called."""

                generate_bulk_csv(values)
                sys.exit()

        parser = argparse.ArgumentParser(prog='QUBEKit Analysis', formatter_class=argparse.RawDescriptionHelpFormatter,
                                         description='')

        parser.add_argument('-t', '--temp', type=float)
        parser.add_argument('-n', '--nmol', type=int)
        parser.add_argument('-a', '--analyse', type=str, choices=['both', 'liquid', 'gas'])
        parser.add_argument('-fb', '--forcebalance', type=string_to_bool, choices=[True, False])

        groups = parser.add_mutually_exclusive_group()

        groups.add_argument('-csv', '--csv_filename', action=CSVAction)
        groups.add_argument('-bulk', '--bulk_run')

        return parser.parse_args()

    def handle_bulk(self):

        bulk_data = mol_data_from_csv(self.args.bulk_run)
        bulk_home = os.getcwd()

        if self.args.forcebalance:
            try:
                os.mkdir('targets')
                os.mkdir('forcefield')
            except FileExistsError:
                pass

        if 'results.csv' not in os.listdir(bulk_home):
            self.start_results_file(os.path.join(bulk_home, 'results.csv'))

        # TODO Flip loops to preserve order
        for root, dirs, files in os.walk(bulk_home, topdown=True):
            for di in dirs:
                for name in list(bulk_data):
                    if f'QUBEKit_{name}_' in di:
                        os.chdir(os.path.join(root, di, 'final_parameters'))
                        mol_home = os.getcwd()
                        try:
                            self.molecule = Molecule()
                            self.molecule.bulk_run = True
                            self.molecule.name = name
                            self.molecule.home = mol_home
                            self.molecule.csv_path = os.path.join(bulk_home, 'results.csv')

                            for key, val in bulk_data[name].items():
                                setattr(self.molecule, key, val)

                            if self.args.forcebalance:
                                self.molecule.forcebalance = True
                                try:
                                    os.mkdir(os.path.join(bulk_home, 'targets', f'{self.molecule.name}_liquid'))
                                except FileExistsError:
                                    pass

                            Execute(self.molecule)

                        except Exception as exc:

                            with open(os.path.join(bulk_home, 'errors.txt'), 'a+') as error_file:
                                error_file.write(f'{name}\n\n{exc}\n\n')

                            print(str(exc))
                            continue
                        os.chdir(bulk_home)

        sys.exit()

    @staticmethod
    def start_results_file(path):

        with open(path, 'w+') as results_file:
            file_writer = csv.writer(results_file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            file_writer.writerow(['name', 'liq_temp', 'liq_energy', 'density', 'gas_temp', 'gas_energy', 'heat_of_vap'])


class Execute:

    def __init__(self, molecule):

        self.molecule = molecule
        self.molecule.name, self.molecule.switch_dist, self.total_atoms = self.mol_data()
        self.single_analysis()

    def mol_data(self):
        """
        Find the name of the molecule and the appropriate switching distance.
        This is run from inside the finalise folder so all info is present in the pdb or xml files.
        """

        name = [file for file in os.listdir('.') if file.endswith('.xml')][0].split('.')[0]

        heavy_atoms = 0
        hydrogens = 0

        with open(f'{name}.xml', 'r') as xml:
            for line in xml:
                # Find atom elements, ignore v-sites
                if ('<Type' in line) and ('v-site' not in line):
                    if line.split('element')[1][2] != 'H':
                        heavy_atoms += 1
                    else:
                        hydrogens += 1
        total_atoms = heavy_atoms + hydrogens

        if heavy_atoms < 3:
            switch_dist = 1.05
        elif heavy_atoms < 5:
            switch_dist = 1.25
        else:
            switch_dist = 1.45

        return name, switch_dist, total_atoms

    @exception_catcher
    def single_analysis(self):
        """
        Use the pdb and xml for the liquid properties calcs.
        temp, nmol etc has already been extracted from the relevant places.
        """

        create_box(f'{self.molecule.name}.pdb', self.molecule.nmol)

        generate_connections(f'{self.molecule.name}.pdb', 'box.pdb', self.molecule.nmol)

        if self.molecule.forcebalance:
            os.system(f'cp new.pdb ../../targets/{self.molecule.name}_liquid/liquid.pdb')
            os.system(f'cp {self.molecule.name}.pdb ../../targets/{self.molecule.name}_liquid/gas.pdb')

            write_fb_csv(f'../../targets/{self.molecule.name}_liquid/')
            return

        if self.molecule.analyse == 'both' or self.molecule.analyse == 'liquid':

            make_and_change_into('liquid')

            os.system(f'cp ../{self.molecule.name}.xml ../new.pdb .')
            liquid_analysis(self.molecule.name, switch_dist=self.molecule.switch_dist, temp=self.molecule.temp)

            os.chdir('../')

        if self.molecule.analyse == 'both' or self.molecule.analyse == 'gas':

            make_and_change_into('gas')

            os.system(f'cp ../{self.molecule.name}.xml ../{self.molecule.name}.pdb .')
            gas_analysis(self.molecule.name, temp=self.molecule.temp)

            os.chdir('../')

        calculate_properties(self.molecule.csv_path, self.molecule.name, self.molecule.nmol, self.total_atoms)


def main():
    ArgsAndConfigs()


if __name__ == '__main__':
    # Don't forget entry point, unlike QUBEKit where it's defined in the setup.py
    main()
