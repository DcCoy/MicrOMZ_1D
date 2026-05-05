# MicrOMZ_1D 

MicrOMZ_1D is a trait-based functional type model that simulates the growth, mortality, and grazing of chemoautotrophic and heterotrophic bacterial functional type populations in a one-dimensional oxygen minimum zone (OMZ) water column (see Zakem et al., 2019, 2022).

By explicitly representing the metabolisms of these bacterial populations, the model resolves key nitrogen-cycle reactions, including nitrification, denitrification, and anammox, using resource competition arguments rather than prescribing reaction rates as functions of seawater chemistry alone (e.g., `NitrOMZ` from Bianchi et al., 2022).

The ecosystem model is embedded in a 1D advection-diffusion-reaction framework with an imposed sinking flux of organic carbon.

Contact Daniel McCoy (dmccoy@carnegiescience.edu) for assistance.
    
## Table of Contents

- [Updates](#updates)
- [Model Description](#model-description)
- [Getting started](#getting-started)
- [Code structure](#code-structure)
- [Support](#support)
- [How to cite](#how-to-cite)

Requires MATLAB 2013 or above.

## Updates

* 05/05/2026 -- First commit of MicrOMZ_1D 

## Model Description

MicrOMZ_1D models microbial functional types competing for oxidants and reductants in a 1D water column with an imposed sinking flux of organic carbon. Each microbial type $i$

- consumes one or more substrates, such as OM, $\mathrm{NO}_3^-$, $\mathrm{NO}_2^-$, $\mathrm{NH}_4^+$, $\mathrm{O}_2$, and $\mathrm{N}_2\mathrm{O}$,
- produces metabolic byproducts determined by its redox pathway,
- experiences losses from mortality and grazing,
- grows according to Liebig’s law of the minimum.

Functional types interact through substrate competition and metabolic cross-feeding. Resource competition shapes local substrate concentrations through the balance of physical supply, biological consumption, and biological production. All chemical tracers evolve through 1D transport, microbial uptake, and microbial byproduct formation.

The model contains 10 functional types that fit under two broad metabolic categories.

Chemoautotrophs derive energy from the oxidation of inorganic substrates and include:

- **aoa** = aerobic ammonium oxidizers ($\mathrm{NH}_4^+ + \mathrm{O}_2 \rightarrow \mathrm{NO}_2^-$)
- **nob** = aerobic nitrite oxidizers ($\mathrm{NO}_2^- + \mathrm{O}_2 \rightarrow \mathrm{NO}_3^-$)
- **aox** = anaerobic ammonium oxidizers ($\mathrm{NH}_4^+ + \mathrm{NO}_2^- \rightarrow \mathrm{N}_2$)

Heterotrophs respire organic matter using either $\mathrm{O}_2$ or oxidized nitrogen species as electron acceptors and include:

- **aer** = obligately aerobic heterotrophs 
- **nar** = $\mathrm{NO}_3^-$ to $\mathrm{NO}_2^-$ reducing heterotrophs 
- **nai** = $\mathrm{NO}_3^-$ to $\mathrm{N}_2\mathrm{O}$ reducing heterotrophs 
- **nao** = $\mathrm{NO}_3^-$ to $\mathrm{N}_2$ reducing heterotrophs 
- **nir** = $\mathrm{NO}_2^-$ to $\mathrm{N}_2\mathrm{O}$ reducing heterotrophs 
- **nio** = $\mathrm{NO}_2^-$ to $\mathrm{N}_2$ reducing heterotrophs 
- **nos** = $\mathrm{N}_2\mathrm{O}$ to $\mathrm{N}_2$ reducing heterotrophs 

---

### Model Equations

#### 1. Biomass dynamics

Microbial biomass $B_i$ (mmol $\mathrm{m}^{-3}$) evolves through the balance of physical transport, growth, mortality, and grazing:

$$
\frac{\partial B_i}{\partial t}
= \mathcal{T}(B_i) + \left(\mu_i - m_i - g_i\right)B_i,
$$

where $\mathcal{T}(B_i)$ represents 1D physical transport, $\mu_i$ is the realized growth rate of type $i$, $m_i$ is the mortality rate, and $g_i$ is the grazing loss rate. A population is favored where its required substrates are sufficiently available relative to its local loss terms and uptake traits.

---

#### 2. Growth limitation and Monod uptake

Potential growth of microbial type $i$ on resource $j$ follows Monod uptake kinetics:

$$
\mu_{ij}
= y_{ij} V^{\max}_{ij} \frac{R_j}{R_j + K_{ij}},
$$

where:

- $K_{ij}$ — half-saturation coefficient,
- $V^{\max}_{ij}$ — maximum uptake rate (mol resource per mol biomass per day),
- $y_{ij}$ — biomass yield (mol biomass per mol resource),
- $R_j$ — environmental concentration of resource $j$.

For functional types requiring more than one resource, realized growth follows Liebig’s law of the minimum:

$$
\mu_i = \min_j \left(\mu_{ij}\right).
$$

A useful local subsistence concentration for resource $j$ can be defined relative to the total biological loss rate $L_i = m_i + g_i$:

$$
R^*_{i,j} = \frac{K_{ij} L_i}{y_{ij} V^{\max}_{ij} - L_i}.
$$

This $R^*$ expression is used as a diagnostic of resource competition among functional types, while the prognostic model evolves biomass and tracers in the 1D water column.

---

#### 3. Tracer dynamics

All dissolved tracers obey a 1D transport-reaction mass balance. For tracer $X$,

$$
\frac{\partial [\mathrm{X}]}{\partial t} = \mathcal{T}([\mathrm{X}]) + \sum_i e_{i,\mathrm{X}} \mu_i B_i - \sum_i \frac{1}{y_{i,\mathrm{X}}} \mu_i B_i,
$$

where $\mathcal{T}([\mathrm{X}])$ represents 1D physical transport, $e_{i,\mathrm{X}}$ is the production coefficient for tracer $X$ by functional type $i$, and $1/y_{i,\mathrm{X}}$ represents the resource consumed per biomass produced when tracer $X$ limits growth.

The generic transport operator includes vertical diffusion and any imposed vertical advection:

$$
\mathcal{T}([\mathrm{X}])
= \frac{\partial}{\partial z}
\left(K_z \frac{\partial [\mathrm{X}]}{\partial z}\right)
- w \frac{\partial [\mathrm{X}]}{\partial z},
$$

with tracer- or particle-specific modifications where appropriate, such as sinking organic matter fluxes.

## Getting started

#### Set run options via `options.m`

Set model options for variable organic matter input:

```matlab
amplitude  % amplitude of organic matter oscillations
period     % period of oscillations
```

`options.m` also controls which tracers and functional types are included in a given run.
        
#### Run the model 

```matlab
model = run_model;
```

## Code structure 

```text
options.m    == Script to toggle main model settings
run_model.m  == Call this to run the model, e.g. >> run_model
```
                   
#### functions/

Folder where utility functions are stored, mostly for plotting.

#### src/

Folder where core model functions are stored:

```text
init_model     == Initializes model based on options.m
timestepping   == Runs the model forward in time, saving each timestep
sources_sinks  == Calculates biogeochemical sources and sinks called by timestepping.m
calc_rstar     == Calculates Tilman's R* diagnostics for each functional type
ftype_params   == Stores functional type parameters
```

#### plotting/ 

Folder where useful plotting functions for MicrOMZ output are stored:

```text
make_plots  == Generates model spinup figures for each tracer
make_zngi   == Generates zero-net-growth-isolines for oxidant/reductant pairs
```

#### output/

Folder to store output.
      
## Support

Contact Daniel McCoy at the Carnegie Institution for Science (dmccoy@carnegiescience.edu).

## How to cite

Citation information will be added here.
