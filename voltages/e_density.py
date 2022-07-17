def get_energy_density(voltage, mol_weight, x, watt_hour: bool):
    F = 96485.3321
    MTSE = ((x * voltage) / mol_weight) * F

    # check unit
    unit = 'kj/kg'
    if watt_hour:
        MTSE = MTSE / 3.6
        unit = 'Wh/kg'

    print('MTSE =', MTSE, unit)
    return MTSE

# from voltages import get_eq_voltage
# eq_v = get_eq_voltage('PBE')
# Li2O_weight = 29.882
# x = 1.5
#
# get_energy_density(eq_v, Li2O_weight, x, watt_hour=True)

