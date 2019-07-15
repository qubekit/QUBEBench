import csv
import os


def generate_bulk_csv(csv_name):

    log_files = []
    for root, dirs, files in os.walk('.', topdown=True):
        for file in files:
            if 'QUBEKit_log.txt' in file and 'backups' not in root:
                log_files.append(os.path.abspath(f'{root}/{file}'))

    if not log_files:
        print('No QUBEKit directories with log files found. Perhaps you need to move to the parent directory.')
        return

    names = []
    for file in log_files:
        with open(file) as log_file:
            for line in log_file:
                if 'Analysing:' in line:
                    names.append(line.split()[1])
                    break

    with open(csv_name, 'w') as csv_file:

        file_writer = csv.writer(csv_file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        file_writer.writerow(['name', 'temp', 'nmol', 'analyse'])

        for name in names:
            file_writer.writerow([name, '', '', ''])

    print(f'{csv_name} generated.', flush=True)
    return


def mol_data_from_csv(csv_name):

    with open(csv_name) as csv_file:

        mol_confs = csv.DictReader(csv_file)

        rows = []
        for row in mol_confs:

            row = dict(row)
            row['temp'] = float(row['temp']) if row['temp'] else 298.15
            row['nmol'] = int(row['nmol']) if row['nmol'] else 267
            row['analyse'] = row['analyse'] if row['analyse'] else 'both'
            rows.append(row)

    final = {row['name']: row for row in rows}

    for val in final.values():
        del val['name']

    return final
