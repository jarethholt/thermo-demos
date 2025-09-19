import numpy as np
from scipy.special import xlogy

CRITICAL_TEMPERATURE = 1.0
CRITICAL_CHEM_POT = -2.0


def free_energy(density, temperature, chem_pot):
    energy = -density * (2 * density + chem_pot)
    entropy = -xlogy(density, density) - xlogy(1 - density, 1 - density)
    return energy - temperature * entropy


def dfree_energy(density, temperature, chem_pot):
    energy = -4 * density - chem_pot
    entropy = -np.log(density) + np.log(1 - density)
    return energy - temperature * entropy


def ddfree_energy(density, temperature, _chem_pot):
    energy = -4
    entropy = -1 / (density * (1 - density))
    return energy - temperature * entropy


def max_gas_density(temperature):
    return 0.5 * (1 - np.sqrt(1 - temperature))


def min_liq_density(temperature):
    return 0.5 * (1 + np.sqrt(1 - temperature))


def max_gas_chem_pot(temperature):
    rt = np.sqrt(1 - temperature)
    return -2 + 2 * rt + temperature * np.log((2 - temperature - 2 * rt) / temperature)


def min_liq_chem_pot(temperature):
    rt = np.sqrt(1 - temperature)
    return -2 - 2 * rt + temperature * np.log((2 - temperature + 2 * rt) / temperature)
