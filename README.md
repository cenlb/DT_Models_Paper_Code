# DT Models Pape Code
Code accompanying "When and why direct transmission models can be used for environmentally persistent pathogens"

## Data
Python pickles and .txt files containing simulated data sets used in Section 3 of paper.

## MC cubed sampler 

Python3 implementation of MC cubed (*Parallel Metropolis coupled Markov chain Monte Carlo for Bayesian phylogenetic inference, Altekar, Dwarkadas, Huelsenbeck & Ronquist, 2003*)

to sample from posterior

<img src="https://render.githubusercontent.com/render/math?math=p( \beta, \delta, \mathbf{t}^E | \mathbf{t}^I, \mathbf{t}^R )">.  

We can write down the marginal posterior posterior density for the mortality rate, <img src="https://render.githubusercontent.com/render/math?math=\gamma">, so there is no need to draw MCMC samples for this.  We get posterior samples for the direct transmission rate, <img src="https://render.githubusercontent.com/render/math?math=\beta">, first by calculcating

<img src="https://render.githubusercontent.com/render/math?math=A = \int_{t^E_\kappa}^\infty S_t I_t \, dt = \sum_{j=1}^m \sum_{i=1}^N \{ \min(t_j^R, t^E_i) - \min(t_j^I, t_i^E)\}">

for each MCMC draw of <img src="https://render.githubusercontent.com/render/math?math=\mathbf{t}^E">

then drawing 

<img src="https://render.githubusercontent.com/render/math?math=\beta \sim \Gamma( m \text{ plus } \nu_{\beta} \text{ minus } 1, A \text{ plus } \lambda_{\beta} )">

Set parameters in `sample_parallel.py` to adjust block size (number of iterations between attempted chain switches) and blocks to burn-in, etc.

While burning in, stdout shows scale parameter for proposal density.  After burning in, stdout shows information on chain switches.

Note: this is not an ideal implementation of MC_cubed, for example, all chains wait at switch points, regardless of whether they are switch candidates!
  
## SEIR / SEIR-P simulation
  
`simulate.py` implements full Doob-Gillespie to simulate SIR / SIR-P / SEIR / SEIR-P processes.  Set `E_` and `I_` definitions to draw required exposed and infectious lifetime distributions.
  
`wssv.cc` does the same (and more) but much faster. Approximates the hazard function by discretising and forward simulateing pathogen load exactly using one draw from a binomial and one draw from a Poission distribution at each time step.  For approximation to be accurate, set time step to a suitably small value (one order of magnitude small than pathogen's mean lifetime should do it).  
 
compile with
  
g++ -std=c++11 -o wssv wssv.cc
  
g++ -std=c++11 -o wssv_SEIR wssv_SEIR.cc  
