# monte_carlo_tissue.py
# Minimal runnable Monte Carlo depth-dose simulator for layered tissues.
# Usage: python monte_carlo_tissue.py

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class Layer:
    name: str
    thickness_cm: float
    mu_a_cm: float
    mu_s_cm: float

@dataclass
class Geometry:
    layers: List[Layer]

    def total_thickness(self) -> float:
        return sum(layer.thickness_cm for layer in self.layers)

    def locate_layer(self, z: float):
        cum = 0.0
        for i, L in enumerate(self.layers):
            next_cum = cum + L.thickness_cm
            if z < next_cum - 1e-12:
                return i, z - cum, L
            cum = next_cum
        return len(self.layers) - 1, self.layers[-1].thickness_cm, self.layers[-1]

def simulate_depth_dose(geom: Geometry, n_photons=100000, photon_energy=1.0, dz=0.01, rng_seed=42, max_interactions=10000):
    rng = np.random.default_rng(rng_seed)
    z_max = geom.total_thickness()
    n_bins = int(np.ceil(z_max / dz))
    dose_bins = np.zeros(n_bins, dtype=float)
    absorbed_count = 0
    transmitted_count = 0
    for _ in range(n_photons):
        z = 0.0
        energy = photon_energy
        interactions = 0
        while True:
            interactions += 1
            if interactions > max_interactions:
                break
            if z >= z_max - 1e-12:
                transmitted_count += 1
                break
            li, z_in_layer, L = geom.locate_layer(z)
            mu_t = L.mu_a_cm + L.mu_s_cm
            if mu_t <= 0.0:
                remaining = L.thickness_cm - z_in_layer
                z += remaining
                continue
            s = -np.log(rng.random()) / mu_t
            dist_to_boundary = L.thickness_cm - z_in_layer
            if s >= dist_to_boundary:
                z += dist_to_boundary
                continue
            else:
                z_int = z + s
                if rng.random() < (L.mu_a_cm / mu_t):
                    bin_idx = min(int(z_int // dz), n_bins - 1)
                    dose_bins[bin_idx] += energy
                    absorbed_count += 1
                    break
                else:
                    z = z_int
                    continue
    depth_dose = dose_bins / dz
    z_edges = np.linspace(0, z_max, n_bins + 1)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    return z_centers, depth_dose, absorbed_count, transmitted_count, n_photons

if __name__ == '__main__':
    layers = [
        Layer('Epidermis', 0.1, 0.15, 5.0),
        Layer('Dermis',    0.2, 0.12, 7.0),
        Layer('Fat',       0.5, 0.05, 4.0),
        Layer('Muscle',    1.2, 0.10, 9.0),
        Layer('Tumor',     0.5, 0.20, 12.0),
    ]
    geom = Geometry(layers=layers)
    z, dose, absorbed, transmitted, N = simulate_depth_dose(geom, n_photons=80000, photon_energy=1.0, dz=0.01, rng_seed=7)
    import pandas as pd
    df = pd.DataFrame({'depth_cm': z, 'dose_au_per_cm': dose})
    df.to_csv('depth_dose_profile.csv', index=False)
    plt.figure()
    plt.plot(z, dose)
    plt.xlabel('Depth (cm)')
    plt.ylabel('Deposited energy per unit depth (a.u./cm)')
    plt.title('Monte Carlo Depth-Dose in Layered Tissue')
    plt.tight_layout()
    plt.savefig('depth_dose_profile.png', dpi=150)
    print(f'Photons: {N} | Absorbed: {absorbed} | Transmitted: {transmitted}')
