# DT_Models_Paper_Code
Code accompanying "When and why direct transmission models can be used for environmentally persistent pathogens"

## data
Python pickles of simulated data sets used in Section 3 of paper.

## MC_cubed_sampler 

Contains Python3 code to sample from posterior $p( \beta, \deta, \mathbf{t}^E | \mathbf{t}^I, \mathbf{t}^R ) $.  
$\delta$'s marginal posterior can be written down, so no need to sample.

Marginal posterior for $\beta$ obtained after running sampler then drawing from $\Gamma( m + \nu_\beta - 1, A - \lambda_\beta)$, (A in output file "A_<suffix>.txt").  Set parameters in sample_parallel.py to adjust block size (number of iterations between attempted chain switches) and blocks to burn-in, etc.

While burning in, stdout shows scale parameter for proposal density.  After burning in, stdout shows information on chain switches.

Note: this is not an ideal implementation of MC_cubed!
  
## wssv_simulation
 
compile with
  
g++ -std=c++11 -o wssv wssv.cc
  
g++ -std=c++11 -o wssv_SEIR wssv_SEIR.cc  
