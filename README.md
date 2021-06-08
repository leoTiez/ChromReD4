# Reaction-Diffusion (ReD) System for Rad4 and Chromatin

This is a simple simulation of a reaction-diffusion system which models how Rad4--involved
in the damage detection in the DNA repair pathway _GG-NER_ in yeast--could potentially probe 
for lesions in the genome. In wildtype strains, Rad4 seems to be unevenly distributed in the nucleus and colocalises with heterochromatin (for human cells, see [1]). Hoogstraten et al. proposed that XPC (the human homologue of Rad4) can only associate transiently to the DNA molecule. After radiation and through the induction of DNA damage, XPC binds more stably to the damage site [1]. Therefore, Rad4 randomly probes for lesions rather than scanning for it. This is supported by this simple model of reaction-diffusion system. Rad4 randomly associates and dissociates to the DNA. After radiation, the dissociation probability is reduced, resulting in a stable DNA:protein interaction. As it can be verified by the model, this is decreasing the concentration gradient, resulting in the dissolution of the Rad4 patches. This allows the distribution of Rad4 from heterochromatin to euchromatin.

## Example
![ChromReD4 Example](ChromReD4_example.gif)
In this simulation, I used reflecting boundaries for creating the chromatin and wrapping boundaries for the reaction-diffusion system. Red represents high concentration of Rad4, blue shades represent low concentrations. The values do not represent real biological values and rather should exemplify how the system works.

## Installation and Usage
The model is implemented in Julia. Run from the project directory
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

## References
[1] Hoogstraten, Deborah, et al. "Versatile DNA damage detection by the global genome nucleotide excision repair protein XPC." Journal of cell science 121.17 (2008): 2850-2859.