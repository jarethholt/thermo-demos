import scipy.optimize as spo

import thermo


def _find_state(temperature, chem_pot, *, bounds):
    result = spo.minimize_scalar(
        thermo.free_energy, bounds=bounds, args=(temperature, chem_pot)
    )
    return result.x, result.fun


def find_supercritical_state(temperature, chem_pot):
    if temperature < 1:
        raise ValueError("No supercritical state for temperature < 1")
    return _find_state(temperature, chem_pot, bounds=[0, 1])


def _find_gas_state(temperature, chem_pot):
    d_gas = thermo.max_gas_density(temperature)
    return _find_state(temperature, chem_pot, bounds=[0, d_gas])


def _find_liq_state(temperature, chem_pot):
    d_liq = thermo.max_liq_density(temperature)
    return _find_state(temperature, chem_pot, bounds=[d_liq, 1])


def find_subcritical_states(temperature, chem_pot):
    if temperature > 1:
        raise ValueError("No subcritical state for temperature > 1")
    states = []
    if chem_pot <= thermo.max_gas_chem_pot(temperature):
        states.append(_find_gas_state(temperature, chem_pot))
    if chem_pot >= thermo.min_liq_chem_pot(temperature):
        states.append(_find_liq_state(temperature, chem_pot))
    return states


def _get_energy(state):
    return state[1]


def find_equilibrium_state(temperature, chem_pot):
    if temperature >= 1:
        return find_supercritical_state(temperature, chem_pot)
    states = find_subcritical_states(temperature, chem_pot)
    return min(states, key=_get_energy)
