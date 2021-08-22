# Simulated dataset

* [sim.R](https://github.com/Wancen/airpartpaper/blob/main/sim/sim.R): To assess airpart’s partitioning of cell types by allelic ratio. Half of the total counts were drawn from a NB with a low mean count of 2 while half of the total counts were drawn from a NB with a higher mean count that ranged across different simulations among values of $cnt \in \{5, 10, 20\}$. The number of gene within a gene cluster was varied across $g \in \{5,10,20\}$.
* [comparison.R](https://github.com/Wancen/airpartpaper/blob/main/sim/comparison.R): For evaluating estimation accuracy as well as assessing the effect of similar cell types partition by comparing  with group step and  without group step
* [DAI_test.R](https://github.com/Wancen/airpartpaper/blob/main/sim/DAI_test.R): For testing DAI statistical significance, we simulated one without DAI \{0.5,0.5,0.5,0.5,0.5,0.5\} and with DAI \{0.5,0.5,0.6,0.6,0.7,0.7\} on 6 cell types
