import pandas as pd
import csv


def calculate_properties(csv_path, name, nmol, num_atoms, liq=True, gas=True):
    """
    Script to calculate the liquid properties of an OpenMM simulation.
    Both liquid and gas results files are needed.
    """
    desired_temp = 298.15
    gas_const = 0.0019858775  # Should not be changed

    # The next section calculates the amount of results used in the calculation based on simulation time
    # It will exclude the first nanosecond of data for both simulations based on the time step used

    time_required = 2  # Enter the number of nanoseconds the results should be calculated over
    # Provided there are enough results in the files

    liq_time_step = 0.001  # Enter the value of the time step used in the liquid simulation in picoseconds
    liq_report_step = 1000  # The number of actual steps between csv data points liquid
    liq_exclude_steps = 1000000 / (liq_time_step * liq_report_step)
    liq_step_limit = (time_required * liq_exclude_steps) + liq_exclude_steps + liq_report_step

    gas_time_step = 0.0005  # Enter the value of the time step used in the gas simulation in picoseconds
    gas_report_step = 1000  # The number of actual steps between csv data points gas
    gas_exclude_steps = 1000000 / (gas_time_step * gas_report_step)
    gas_step_limit = (time_required * gas_exclude_steps) + gas_exclude_steps + gas_report_step

    if liq:
        with open('liquid/liquid.txt', 'r') as f_in:
            lines = f_in.readlines()
            lines[0] = lines[0].replace('#', '')  # Get rid of # at the start of file and rename as liquid.csv

        with open('liquid.csv', 'w+') as f_out:
            for line in lines:
                f_out.write(line)

        df_a = pd.read_csv('liquid.csv')
        df_a = df_a[df_a.Step < liq_step_limit]
        df_a = df_a[df_a.Step > liq_exclude_steps]

        # Find properties
        liq_temp = df_a['Temperature (K)'].mean()
        print(f"Average liquid temperature (K) = {liq_temp: .5f}")
        liq_energy = df_a['Potential Energy (kJ/mole)'].mean() / (4.184 * nmol)
        print(f"Average liquid energy (kcal / mol) = {liq_energy: .5f}")
        density = df_a['Density (g/mL)'].mean()
        print(f"Average liquid density (g / cc) = {density: .5f}")

    if gas:
        # Repeat now for the gas file
        with open('gas/gas.txt', 'r') as f_in:
            lines = f_in.readlines()
            lines[0] = lines[0].replace('#', '')

        with open('gas.csv', 'w+') as f_out:
            for line in lines:
                f_out.write(line)

        df_b = pd.read_csv('gas.csv')
        df_b = df_b[df_b.Step < gas_step_limit]
        df_b = df_b[df_b.Step > gas_exclude_steps]

        gas_temp = df_b['Temperature (K)'].mean()
        print(f"Average gas temperature (K) = {gas_temp: .5f}")
        gas_energy = df_b['Potential Energy (kJ/mole)'].mean() / 4.184
        print(f"Average gas energy (kcal / mol) = {gas_energy: .5f}")

    # Calculate the heat of vap
    if liq and gas:
        temperature_diff = liq_temp - gas_temp
        heat_energy = gas_energy - liq_energy
        heat_of_vap = heat_energy + 0.5 * gas_const * temperature_diff * (3 * num_atoms - 6) + gas_const * desired_temp

        print(f"Heat of vap (kcal/mol) = {heat_of_vap: .5f}")

        # Test on other files and add where the heat of vap equation is from

        with open(csv_path, 'a+') as results:
            file_writer = csv.writer(results, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            file_writer.writerow([name, liq_temp, liq_energy, density, gas_temp, gas_energy, heat_of_vap])
