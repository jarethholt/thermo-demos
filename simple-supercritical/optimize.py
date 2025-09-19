from dataclasses import dataclass
from operator import attrgetter

import numpy as np
import scipy.optimize as spo

import thermo

Scalar = float | np.floating


@dataclass
class State:
    density: np.float64
    free_energy: np.float64


def _find_state(
    temperature: Scalar, chem_pot: Scalar, *, bounds: tuple[Scalar, Scalar]
) -> State:
    result = spo.minimize_scalar(
        thermo.free_energy, bounds=bounds, args=(temperature, chem_pot)
    )
    return State(density=result.x, free_energy=result.fun)


def find_supercritical_state(temperature: Scalar, chem_pot: Scalar) -> State:
    if temperature < 1:
        raise ValueError("No supercritical state for temperature < 1")
    return _find_state(temperature, chem_pot, bounds=(0.0, 1.0))


def _find_gas_state(temperature: Scalar, chem_pot: Scalar) -> State:
    d_gas = thermo.max_gas_density(temperature)
    return _find_state(temperature, chem_pot, bounds=(0.0, d_gas))


def _find_liq_state(temperature: Scalar, chem_pot: Scalar) -> State:
    d_liq = thermo.min_liq_density(temperature)
    return _find_state(temperature, chem_pot, bounds=(d_liq, 1.0))


def find_subcritical_states(temperature: Scalar, chem_pot: Scalar) -> list[State]:
    if temperature > 1:
        raise ValueError("No subcritical state for temperature > 1")
    states = []
    if chem_pot <= thermo.max_gas_chem_pot(temperature):
        states.append(_find_gas_state(temperature, chem_pot))
    if chem_pot >= thermo.min_liq_chem_pot(temperature):
        states.append(_find_liq_state(temperature, chem_pot))
    return states


def find_equilibrium_state(temperature: Scalar, chem_pot: Scalar) -> State:
    if temperature >= 1:
        return find_supercritical_state(temperature, chem_pot)
    states = find_subcritical_states(temperature, chem_pot)
    return min(states, key=attrgetter("free_energy"))
