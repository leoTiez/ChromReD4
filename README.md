# Reaction-Diffusion (ReD) System for Rad4 and Chromatin

This is a simple simulation of a reaction-diffusion system which models how Rad4--involved
in the damage detection in the DNA repair pathway _GG-NER_ in yeast--could potentially probe 
for lesions in the genome. In wildtype strains, the Rad4-homologue in humans XPC seems to be unevenly distributed in the nucleus and colocalises with heterochromatin [1]. Hoogstraten et al. proposed that XPC can only associate transiently to the DNA molecule. After radiation and through the induction of DNA damage, XPC binds more stably to the damage site [1]. Therefore, Rad4 randomly probes for lesions rather than scanning for it. This is supported by this simple model of reaction-diffusion system. Rad4 randomly associates and dissociates to the DNA. After radiation, the dissociation probability is reduced, resulting in a stable DNA:protein interaction. As it can be verified by the model, this is decreasing the concentration gradient, resulting in the dissolution of the Rad4 patches. This allows the distribution of Rad4 from heterochromatin to euchromatin.

## Example
![ChromReD4 Example](ChromReD4_example.gif)

In this simulation, I used reflecting boundaries for creating the chromatin and wrapping boundaries for the reaction-diffusion system. Red represents high concentration of Rad4, blue shades represent low concentrations. The values do not represent real biological values and rather should exemplify how the system works.

## Installation and Usage of Julia scripts
The repo provides two implementations. The Julia model assumes continuous protein levels and uses the Einstein diffusion,
which is therefore a deterministic diffusion model. Install dependencies via

```commandline
julia --project=. -e 'using Pkg; Pkg.instantiate()
```

Run from the project directory
```commandline
julia --project=. src/main.jl 
```
to run the example. If you want to use the model components for your own implementation, include it via

```julia
using ChromRed4
```

and use the functions

```julia
len_chromatin = 10000
size_n = [100, 100]
# Pass chromatin length and size of the 2-dim nucleus
nucleus = init_nucleus(len_chromatin, size_n)  
# Run one update step with the association probability of 0.1,
# the dissociation probability of 0.7
# and the diffusion constant of 0.245
update!(nucleus, asso_prob=.1, disso_prob=.7, diffu=.245)  
# Plot the nucleus state with a string that can be passed to add more information to the title.
# The string is appended to the title.
plot(nucleus, "Before Radiation")
```

## Installation and Usage of Python Script
The Python script is exclusively based on random particle dynamics, including diffusion. It only uses a smoothed map
of the particle distribution to determine which proteins can move most. There can be several proteins at one place,
but a protein can only be moved to a free spot through diffusion. An example is given here:

![diffusion](figures/simulations/diffusion_nucleus.gif)

When the protein is now interacting with the DNA, we obtain a reaction diffusion pattern. After 150 minutes, we simulate
that dynamics have changed (for example through UV irradiation), and proteins can bind less stably. We see the pattern 
pattern changing from then onwards. The DNA is created through a self-avoiding random walk.

![visited DNA](figures/simulations/reaction_diffusion_koff_lowtohigh_nucleus.gif)

Interestingly, if the DNA is coloured w.r.t. whether it has been visited or not, we 
observe a local patches behaviour shortly after the change, which becomes subsequently quickly random again.

![reaction diffusion](figures/simulations/reaction_diffusion_koff_lowtohigh_dna.gif)

Code was implemented in `Python3` and tested on `Python3.8`. To install the requirements, run

```commandline
python3.8 -m pip install -r requirements.txt
```

To run the main program, run

```commandline
python3.8 ReD4.py
```

Change the behaviour by changing the parameters in the file.

## Results
Whilst both simulate the same underlying process, they differ in their implementation. The deterministic diffusion equation
calculates difference to a mean when everything would be perfectly mixed. Indeed, a visualisation of such an approach would
indicate diffusion away from the DNA when binding is more stable. However, the Python implementation reveals that actually 
accumulation becomes stronger in genetic regions, when considering also bound proteins. Consequently, the phenomenon of 
apparently diffusing proteins occur simply the majority of proteins that can move freely are stably bound, and are therefore 
taken out of the system. This slightly contradicts [1]. However, it should be emphasised that this implementation was only 
used for exemplifying different protein distributions in the nucleus, and what they mean for DNA-protein interaction processes.
**An in depth analysis is missing, and the results should be taken w/ a pinch of salt.**

## References
[1] Hoogstraten, Deborah, et al. "Versatile DNA damage detection by the global genome nucleotide excision repair protein XPC." Journal of cell science 121.17 (2008): 2850-2859.
