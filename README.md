# DEGBOE: Discrete time Evolution modeling of Gene mutation through Bayesian inference using qualitative Observation of mutation Events

 
** Many direct and indirect methods have been developed to help estimate occurrences and progression of somatic mutations in various organisms. The main difficulty in estimating rates of mutation involves the fact that DNA changes are extremely rare events which can only be detected on a background of identical DNA. On the Basis of the probabilistic graphical model of Bayesian network, this work presents a computational framework which describes the variability in the evolution of the lung cancer gene mutations as a long-term dependent process. The relation between a targeted gene mutation with other gene mutations is examined.**
 

# A quick look to the DEGBOE computational framework
![](graphic/ Fig_1.png)


## (A) Population of cellular dynamics for genetic mutation.

** s1:** Rapid growth to tumor
The mutation spreads through the population and replaces the resistant (normal) cells. 

** s2:** An additional driver mutation increases the population of mutated cells, which coexist with normal cells. 

** s3:** Equilibrium of cell mutation with no novel cell mutation. At this stage, neither the normal cell nor the mutated cell populations vary. 

** s4:** Decreasing mutant cell population.

 ** s5:** Mutants go extinct. When the delay before the emergence of the next mutant is long, the population may attain one of the outcomes presented at time T4 (**s1**, **s2**, **s3**, **s4**, and **s5**). 

## (B) Preprocessing lung cancer mutation binary dataset. 

Each row represents a gene, and each column represents a mutational sample. The lung tumor mutations are represented as a discrete timeline of mutational events such that the appearance of a mutation in the timeline corresponds to its estimated place in the tumor evolution.

## (C) Computational model. 

Graphical model with observation z_k={z_k }_(k=1)^N: lung cancer mutation binary dataset.

Local hidden variables x_k={x_k }_(k=1)^N. It illustrates the temporal evolution of mutation occurrence and evolution in the 100 mutational samples

Global hidden variables φ={φ_m }_(m=1)^(M-1). It portrays the presence or absence of other mutations as the total number of corresponding lung cancer mutations found.


**DEGBOE is applied to model occurrences and the evolution of somatic mutations as well as the mutation rates of the EGFR, KRAS, and TP53 lung cancer gene mutations.**

# Implementation
Use R code lung cancer binary dataset preprocessing
Use Matlab for DEGBOE implementation
![image](https://user-images.githubusercontent.com/105561377/168443351-c92c9725-586b-4407-873a-124834e4658b.png)
