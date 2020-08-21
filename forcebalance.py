import csv


def write_fb_csv(loc=''):
    with open(f'{loc}data.csv', 'w') as csv_file:
        file_writer = csv.writer(csv_file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        file_writer.writerow(['Global', 'rho_denom', '30'])
        file_writer.writerow(['Global', 'hvap_denom', '3'])
        file_writer.writerow(['Global', 'alpha_denom', '1'])
        file_writer.writerow(['Global', 'kappa_denom', '5'])
        file_writer.writerow(['Global', 'cp_denom', '2'])
        file_writer.writerow(['Global', 'eps0_denom', '2'])
        file_writer.writerow(['Global', 'use_cvib_intra', 'FALSE'])
        file_writer.writerow(['Global', 'use_cvib_inter', 'FALSE'])
        file_writer.writerow(['Global', 'use_cni', 'FALSE'])

        file_writer.writerow(['T', 'P', 'MBAR', 'Rho', 'Rho_wt', 'Hvap', 'Hvap_wt'])
        file_writer.writerow(['298.15', '1.0 atm', 'FALSE', '', '1', '', '1'])
