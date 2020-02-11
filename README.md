# Fitness Waves: Interactive Evolution Simulations

## About
This is a Streamlit web app which enables users to run simulations of rapidly evolving, asexual populations (e.g populations of *E. coli* or influenza viruses). These are populations in which multiple mutations are often present at once, competing to eventually take over the population. 

In these simulations, individuals produce a random number of offspring based on their **fitness** value. The fitness of an individual is incremented (or decremented) with the acquisition of random **mutations**, which are also inherited by the individual's offspring. Users can specify a population size, a rate at which mutations occur, as well as a probability distribution from which the fitness effect of a new mutation is drawn.

After a simulation runs, the user can interactively visualize two aspects of the population's evolution:

1. **Mutational Trajectories**: Plots the number of individuals carrying each mutation as a function of time. Each mutation is assigned a different color. 
2. **Fitness Distribution**: Plots the distribution of fitnesses in the population as a function of time. As a user mouses over the mutational trajectories plot, the fitness distribution at that time point is plotted as a bar chart. The bar chart is colored according to the mutations carried by the individuals represented (using the same color scheme as in the Mutational Trajectories plot.

## How to Run
The app can be accessed at <http://fitness-waves.herokuapp.com>. Alternatively, the app can be run locally by [installing streamlit](https://docs.streamlit.io/), cloning this repository and running ```streamlit run app.py```.

## Simulation Parameters
Users specify the following simulation parameters using the sidebar:
1. ```N```: Total number of individuals in the population
2. ```U_b```: Rate of beneficial mutations (per individual, per generation)
3. ```U_d```: Rate of deleterious mutations (per individual, per generation)
4. ```s_b```: Average fitness effect size of beneficial mutations. Fitness effects add to the log-fitness of an individual.
5. ```s_d```: Average fitness effect size of deleterious mutations
6. ```rho_b```: Distribution of fitness effects of beneficial mutations. Choices include a single fitness effect, an exponential distribution, or a gamma distribution.
7. ```rho_d```: Distribution of fitness effects of deleterious mutations.
8. ```num_gen```: Number of total generations simulated.
9. ```assay_interval```: The state of the population is recorded for visualization every ```assay_interval``` generations.

## Simulation Steps 

We implement *lineage-based* simulations, which track each of the genetically distinct lineages in the population. Each time a mutation occurs, it seeds a new lineage. We assign each lineage a unique ID, record its fitness, and track the size (number of individuals) within the lineage over time. In large populations, this can be much more computationally feasible than tracking each individual within the population separately. 

Simulations consist of two steps repeated each generation:

1. **Reproduction**: Given lineage sizes ```{n_k(t)}``` at generation ```t``` with corresponding fitnesses ```{X_k}```, lineage sizes at the following generation ```t+1``` are drawn from a Poisson distribution with mean ```n_k(t)\exp\left[X_k - \bar{X(t)} \right] \exp \left[1 - \sum_k{n_k(t)}/N \right]```. The first exponential factor represents the action of selection, which amplifies the size of a lineage based on the difference between its fitness and average fitness of the population. The second exponential factor ensures that the total number of individuals in the population fluctuates around ```N```.
2. **Mutation**: Each lineage acquires a number of beneficial mutations drawn from a binomial distribution with probability ```U_b``` of a mutation per individual. The fitness effect ```s``` of each beneficial mutation is drawn randomly from the probability distribution ```rho_b```. Each beneficial mutation then seeds a new lineage composed of 1 individual with fitness ```X\exp(s)```, where ```X``` is the fitness of the genetic background the mutation landed on. The same process is then repeated for deleterious mutations.

