import numpy as np

def chirp_mass_q_to_components(chirp_mass, mass_ratio):

    total = (chirp_mass) * (1 + mass_ratio)**1.2 / mass_ratio**0.6
    m1 = total / (1+mass_ratio)
    m2 = m1 * mass_ratio

    return m1, m2
