def get_mass_energy_density(voltage, mol_mass, x, watt_hour=True):
    F = 96485.3321
    MTSE = ((x * voltage) / mol_mass) * F

    # check unit
    unit = 'kj/kg'
    if watt_hour:
        MTSE = MTSE / 3.6
        unit = 'Wh/kg'

    print('MTSE =', MTSE, unit)
    return MTSE

def get_volumetric_energy_density(voltage, mol_vol, x, watt_hour=True):
    F = 96485.3321
    VED = ((x * voltage) / mol_vol) * F

    # check unit
    unit = 'kj/L'
    if watt_hour:
        VED = VED / 3.6
        unit = 'Wh/L'

    print('Volumetric energy density =', VED, unit)
    return VED

def get_specific_capacity(mol_mass, x):
    F = 96485.3321
    SC = ((x * F) / (3.6*mol_mass))

    print('Specific capacity =', SC, 'Ah/kg')
    return SC

