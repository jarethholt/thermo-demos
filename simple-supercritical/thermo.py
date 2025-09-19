from typing import overload

import numpy as np
from scipy.special import xlogy

CRITICAL_TEMPERATURE = 1.0
CRITICAL_CHEM_POT = -2.0

Scalar = float | np.floating
Array = np.ndarray[tuple[int, ...], np.dtype[np.floating]]
Numeric = Scalar | Array
ArrayOut = np.ndarray[tuple[int, ...], np.dtype[np.float64]]


@overload
def free_energy(density: Scalar, temperature: Scalar, chem_pot: Scalar) -> np.float64:
    pass


@overload
def free_energy(density: Array, temperature: Numeric, chem_pot: Numeric) -> ArrayOut:
    pass


def free_energy(density, temperature, chem_pot):
    energy = -density * (2 * density + chem_pot)
    entropy = -xlogy(density, density) - xlogy(1 - density, 1 - density)
    return energy - temperature * entropy


@overload
def dfree_energy(density: Scalar, temperature: Scalar, chem_pot: Scalar) -> np.float64:
    pass


@overload
def dfree_energy(density: Array, temperature: Numeric, chem_pot: Numeric) -> ArrayOut:
    pass


def dfree_energy(density, temperature, chem_pot):
    energy = -4 * density - chem_pot
    entropy = -np.log(density) + np.log(1 - density)
    return energy - temperature * entropy


@overload
def ddfree_energy(
    density: Scalar, temperature: Scalar, _chem_pot: Scalar
) -> np.float64:
    pass


@overload
def ddfree_energy(density: Array, temperature: Numeric, _chem_pot: Numeric) -> ArrayOut:
    pass


def ddfree_energy(density, temperature, _chem_pot):
    energy = -4
    entropy = -1 / (density * (1 - density))
    return energy - temperature * entropy


@overload
def max_gas_density(temperature: Scalar) -> np.float64:
    pass


@overload
def max_gas_density(temperature: Array) -> ArrayOut:
    pass


def max_gas_density(temperature):
    return 0.5 * (1 - np.sqrt(1 - temperature))


@overload
def min_liq_density(temperature: Scalar) -> np.float64:
    pass


@overload
def min_liq_density(temperature: Array) -> ArrayOut:
    pass


def min_liq_density(temperature):
    return 0.5 * (1 + np.sqrt(1 - temperature))


@overload
def max_gas_chem_pot(temperature: Scalar) -> np.float64:
    pass


@overload
def max_gas_chem_pot(temperature: Array) -> ArrayOut:
    pass


def max_gas_chem_pot(temperature):
    rt = np.sqrt(1 - temperature)
    return -2 + 2 * rt + temperature * np.log((2 - temperature - 2 * rt) / temperature)


@overload
def min_liq_chem_pot(temperature: Scalar) -> np.float64:
    pass


@overload
def min_liq_chem_pot(temperature: Array) -> ArrayOut:
    pass


def min_liq_chem_pot(temperature):
    rt = np.sqrt(1 - temperature)
    return -2 - 2 * rt + temperature * np.log((2 - temperature + 2 * rt) / temperature)
