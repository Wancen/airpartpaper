# airpartpaper
Code for analyses in airpart paper

* [Larsson2019](https://htmlpreview.github.io/?https://github.com/Wancen/airpartpaper/blob/main/Larsson2019/Larsson2019.html): [Larsson et al.(2019)](https://www.nature.com/articles/s41586-018-0836-1)  contains  224mouse embryo stem cells (C57BL/6×CAST/EiJ) and 188 mouse embryofibroblasts (CAST/EiJ×C57BL/6J) including across states of cell cycle(G1, S, G2M), as identified by the author
  
* [Deng2014](https://htmlpreview.github.io/?https://github.com/Wancen/airpartpaper/blob/main/Deng2014/Deng2014.html): [Deng et al.(2014)](https://www.nature.com/articles/s41586-018-0836-1) includes 286 pre-implantation mouse embryo cells composed of10 cell types from an F1 cross of female CAST/EiJ and male C57BL/6Jmice. Cells were sampled along a time course from the zygote and early2-cell stages through the late blastocyst stage of development.

* [spatialDmelxsim](https://github.com/Wancen/airpartpaper/blob/main/spatialDmelxsim/spatialDmelxsim.R): [Combsand Fraser (2018)](https://doi.org/10.1371/journal.pgen.1007631) performed RNA-seq of five hybrid Drosophila species D.melanogaster×D.simulans  embryos  sliced  along  their  anterior-posterior axis to identify genes with spatially varying allelic imbalance.

* [sim](https://github.com/Wancen/airpartpaper/blob/main/sim): 
   * [sim.R](https://github.com/Wancen/airpartpaper/blob/main/sim/sim.R): To assess airpart’s partitioning of cell types by allelic ratio. Half of the total counts were drawn from a NB with a low mean count of 2 while half of the total counts were drawn from a NB with a higher mean count that ranged across different simulations among values of $cnt \in \{5, 10, 20\}$. The number of gene within a gene cluster was varied across $g \in \{5,10,20\}$.
   * [comparison.R](https://github.com/Wancen/airpartpaper/blob/main/sim/comparison.R): For evaluating estimation accuracy as well as assessing the effect of similar cell types partition by comparing  with group step and  without group step
   * [DAI_test.R](https://github.com/Wancen/airpartpaper/blob/main/sim/DAI_test.R): For testing DAI statistical significance, we simulated one without DAI \{0.5,0.5,0.5,0.5,0.5,0.5\} and with DAI \{0.5,0.5,0.6,0.6,0.7,0.7\} on 6 cell types
