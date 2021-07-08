import csv
import json
import os


def generate_bulk_csv(csv_name):
    """Find results jsons for each execution and add the appropriate data for a qubebench run."""
    workflows = []
    for root, dirs, files in os.walk('.', topdown=True):
        for file in files:
            if 'workflow_result.json' in file and 'backups' not in root:
                workflows.append(os.path.join(root, file))

    if not workflows:
        raise FileNotFoundError(
            'No QUBEKit directories with completed workflow files found. '
            'Perhaps you need to move to the parent directory, '
            'or the QUBEKit executions were unsuccessful.'
        )

    names = []
    for file in workflows:
        with open(file) as results_file:
            results_json = json.load(results_file)
            names.append(results_json['input_molecule']['name'])

    with open(csv_name, 'w') as csv_file:
        file_writer = csv.writer(csv_file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        file_writer.writerow(['name', 'temp', 'nmol', 'analyse'])

        for name in names:
            file_writer.writerow([name, '', '', ''])

    print(f'{csv_name} generated.', flush=True)


def mol_data_from_csv(csv_name):

    with open(csv_name) as csv_file:

        mol_confs = csv.DictReader(csv_file)

        rows = []
        for row in mol_confs:

            row = dict(row)
            row['temp'] = float(row['temp']) if row['temp'] else 298.15
            row['nmol'] = int(row['nmol']) if row['nmol'] else 500
            row['analyse'] = row['analyse'] if row['analyse'] else 'both'
            rows.append(row)

    final = {row['name']: row for row in rows}

    for val in final.values():
        del val['name']

    return final


def make_and_change_into(name):
    try:
        os.mkdir(name)
    except FileExistsError:
        pass
    finally:
        os.chdir(name)
