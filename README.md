# Monte-Carlo-Simulation-of-Radiation-Energy-Deposition
Model how radiation (photons) deposit energy in layered tissues (skin, fat, muscle, tumor) using a simple Monte Carlo approach. 

Background
Radiation transport can be modeled as photons undergoing scattering + absorption.
Each tissue layer has a linear attenuation coefficient (μ) which governs interaction probability.
Monte Carlo approach → simulate random photon paths, collect statistics of energy deposition.

Equations
Free path length:

𝐿 = − (1/𝜇) ln(𝑅)
where R∼U(0,1).

Transmission probability after thickness 

I(d)=I0. e^−μd

Absorbed dose (per layer): proportional to number of photons absorbed × energy per photon.

Step Plan
Define layers: thickness + μ (pick simple values).
Launch N photons at surface (e.g., 10⁵).
For each photon:
Sample path length.
If path length < layer thickness → absorbed, deposit energy.
Else → continue to next layer with reduced energy.
Repeat until photon escapes or absorbed.
Aggregate → dose vs. depth profile.

Coding Outline (Python)
Inputs: Layer structure (list of μ, thickness), photon count.
Loop: For each photon, simulate travel until termination.
Outputs: Depth-dose curve (plot: energy deposited vs. depth).

Contents
monte_carlo_tissue.py → main code.
notebook.ipynb → visual explanation with plots.
Figures: Dose vs depth
