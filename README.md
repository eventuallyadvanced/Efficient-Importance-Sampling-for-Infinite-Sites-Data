# Efficient-Importance-Sampling-for-Infinite-Sites-Data
Computing the likelihood of real biological data-sets consistent with the infinite sites model commonly requires
intervention of MCMC or impor- tance sampling based methods. This is due to the computational challenge posed by 
evaluating the exact likelihood of such data.  In this project we focus on importance sampling based methods for 
likelihood approxima- tion  of  infinite  sites  data  (in  the  absence  of  recombination).   A  newly developed
counting scheme allows us to compute the number of histories of a given data set quickly.  We develop the
’Combinatorial Importance Sampler’.  The proposal distribution of the new sampling scheme assigns probabilities to
events directly proportional to the number of histories of the parental configuration corresponding to the event. 
This follows from our conjecture that it is likely that the probability of an observed sample to  have  descended
from  a  particular  ancestral  configuration  is  strongly correlated to the number of histories of that configuration.
We implement the counting scheme and the ’Combinatorial Importance Sampler’ as well as  three  previously  implemented
schemes,  namely  the  Griffiths-Tavar ́e, Stephens-Donnelly and Hobolth importance samplers.  The four schemes are analysed
on data sets of different sizes, including simulated data.  A series  of  experiments  are  conducted  to  evaluate  the
performance  of  the samplers.  Additionally, we attempt to understand the conditions under which the Combinatorial Importance
Sampler fails to achieve desired results.
