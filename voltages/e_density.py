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

