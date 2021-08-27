import numpy as np
import multiprocessing as m_proc
import time
import sys
import SEIR_Sample as samp
import SEIR_Parallel as pll
import pickle

if len(sys.argv) != 3:
    print("Usage: python sample_parallel.py <path> <output_suffix>")
    sys.exit()

path = sys.argv[1]
sufx = sys.argv[2]

'''
main thread RNG
'''
main_ranD = np.random.RandomState(int(time.time()) - 1)

'''
Load data
'''
N = 300 ### Change this ###

fI = open(path + 'I_','rb')
fR = open(path + 'R_','rb')


I = np.load(fI)
R = np.load(fR)

print(R.size)

'''
Priors
Pattern = _beta, _gamma, _delta
'''
#prior_nu = 1.0, 1.0, 12*12/10000
#prior_lambda = 0.001, 0.001, 12/10000

prior_nu = 1.0, 1.0, 1.0
prior_lambda = 0.001, 0.001, 0.001


'''
Initial parameters
'''
parameter_inits = (1.0, 1.0, 1.0)
E_0 = 1.1

'''
MH parameter
'''
sigma_delta_init = 0.01

'''
Number of iterations, block size and length of burn-in
'''
block_size = 400
blocks_to_burn = 20
num_blocks = 4000
thin = 40

temps = [1.00, 1.02, 1.04, 1.06, 1.08, 1.10]
num_chains = len(temps)
write = [True]
write = write + [False for j in range(num_chains - 1)]
seeds = [int(time.time()) + j for j in range(num_chains)]

print(temps, num_chains)
print(write)
print(seeds)

beta_ = open(path + 'beta_' + sufx + '.txt','w')
gamma_ = open(path + 'gamma_' + sufx + '.txt','w')
delta_ = open(path + 'delta_' + sufx + '.txt','w')
K_ = open(path + 'K_' + sufx + '.txt','w')
A_ = open(path + 'A_' + sufx + '.txt','w')
E_kappa_ = open(path + 'E_kappa_' + sufx + '.txt','w')


file_handles = beta_, delta_, gamma_, K_, A_, E_kappa_
burning = True

'''
stuff required for synchronisation
'''
parent_connections = []
child_connections = []

sh_temperatures = m_proc.Array('d',temps)
sh_cross_temperatures = m_proc.Array('d',temps)

sh_likelihoods = m_proc.Array('d',num_chains)
sh_cross_likelihoods = m_proc.Array('d',num_chains)

sh_sigma_delta = m_proc.Array('d',[sigma_delta_init for j in range(num_chains)])

sh_write = m_proc.Array('H',write)


for j in range(num_chains):
    p_conn, c_conn = m_proc.Pipe()
    parent_connections.append(p_conn)
    child_connections.append(c_conn)

print(parent_connections)
print(child_connections)

reusable_barrier = pll.ReusableBarrier( num_chains + 1 ) #+1 for _main_ process
mutex = m_proc.Semaphore(1)

'''
initialise SamplerThreads
'''

sampler_threads = [
    pll.SamplerThread(
        index,
        samp.SEIR_Sampler(np.random.RandomState(seed),
                          prior_nu,
                          prior_lambda, I, R, N,
                          parameter_inits,
                          sigma_delta_init, T, wr),
        reusable_barrier,
        mutex,
        child_conn,
        block_size,
        num_blocks,
        blocks_to_burn,
        thin,
        file_handles,
        sh_temperatures,
        sh_cross_temperatures,
        sh_likelihoods,
        sh_cross_likelihoods,
        sh_sigma_delta,
        sh_write) for index, seed, T, wr, child_conn in zip(range(num_chains), seeds, temps, write, child_connections)]


for th in sampler_threads:
    print(th.seir_sampler.T, th.seir_sampler.write, th.seir_sampler.ranD.uniform(), th.connection)
    th.seir_sampler.initialise(-2.0)
    #print th.seir_sampler.E

'''
Start chains in separate processes
'''    
for p in sampler_threads:
    p.start()

    
'''
Main thread loop - burn-in phase.
Acquire
Reports sigma_delta for each chain after each chain has completed block
'''    
for j in range(blocks_to_burn):
    reusable_barrier.wait_phase_one()
    print(j, [x for x in sh_sigma_delta])
    reusable_barrier.wait_phase_two()
    

'''
Main thread loop - sampling phase
'''

for j in range(num_blocks):
    '''
    1. Chain threads WAIT
    '''
    reusable_barrier.wait_phase_one()
    
    k, l = samp.random_transpose( num_chains, main_ranD )
    
    sh_cross_temperatures[k] = sh_temperatures[l]
    sh_cross_temperatures[l] = sh_temperatures[k]


    '''
    2.
    chain threads:
    read T, sigma_delta and write from sh_<>
    perform <block_size> iterations
    write likelihood to sh_likelihoods
    write cross likelihood to sh_cross_likelihoods
    '''
    reusable_barrier.wait_phase_two()
    
    '''
    3.
    chain threads wait
    '''
    reusable_barrier.wait_phase_one()
    
    r = sh_cross_likelihoods[k] + sh_cross_likelihoods[l]
    r = r - sh_likelihoods[k] - sh_likelihoods[l]

    alpha = np.exp( min(0.0, r) )

    u = main_ranD.uniform()

    if u < alpha:
        print(j, '(',k,',',l,')',[x for x in sh_temperatures], sh_temperatures[k], '<->', sh_temperatures[l])
        #print [x for x in sh_write]
        sh_temperatures[k] = sh_cross_temperatures[k]
        sh_temperatures[l] = sh_cross_temperatures[l]

        hold = sh_sigma_delta[k]
        sh_sigma_delta[k] = sh_sigma_delta[l]
        sh_sigma_delta[l] = hold

        hold = sh_write[k]
        sh_write[k] = sh_write[l]
        sh_write[l] = hold

    else:
        print(j, '(',k,',',l,')',[x for x in sh_temperatures])
        #print [x for x in sh_write]

    reusable_barrier.wait_phase_two()
        
    
    





