import multiprocessing as m_proc
import numpy as np
import time

'''
Little Book of Semaphores
'''
class Barrier:
    def __init__(self, n):
        self.n = n
        self.count = m_proc.Value('i',0)
        self.mutex = m_proc.Semaphore(1)
        self.barrier = m_proc.Semaphore(0)

    def wait(self):
        self.mutex.acquire()
        self.count.value += 1
        self.mutex.release()

        if self.count.value == self.n:
            self.barrier.release()

        self.barrier.acquire()
        self.barrier.release()

class ReusableBarrier:
    def __init__(self, n):
        self.n = n
        self.count = m_proc.Value('i',0)
        self.mutex = m_proc.Semaphore(1)
        self.barrier_one = m_proc.Semaphore(0)
        self.barrier_two = m_proc.Semaphore(1)

    def wait_phase_one(self):
        self.mutex.acquire()
        self.count.value += 1
        if self.count.value == self.n:
            self.barrier_two.acquire() #lock second barrier
            self.barrier_one.release() #unlock first barrier
        self.mutex.release()

        self.barrier_one.acquire()
        self.barrier_one.release()

    def wait_phase_two(self):
        self.mutex.acquire()
        self.count.value -= 1
        if self.count.value == 0:
            self.barrier_one.acquire() #lock first barrier
            self.barrier_two.release() #unlock second barrier
        self.mutex.release()

        self.barrier_two.acquire()
        self.barrier_two.release()

    

class SamplerThread(m_proc.Process):
    def __init__(self,
                 index,
                 sampler,
                 reusable_barrier,
                 mutex,
                 connection,
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
                 sh_write):
        m_proc.Process.__init__(self)
        self.index = index
        self.reusable_barrier = reusable_barrier
        self.mutex = mutex
        self.seir_sampler = sampler
        self.connection = connection
        self.block_size = block_size
        self.num_blocks = num_blocks
        self.blocks_to_burn = blocks_to_burn
        self.thin = thin
        self.file_handles = file_handles
        self.sh_temperatures = sh_temperatures
        self.sh_cross_temperatures = sh_cross_temperatures
        self.sh_likelihoods = sh_likelihoods
        self.sh_cross_likelihoods = sh_cross_likelihoods
        self.sh_sigma_delta = sh_sigma_delta
        self.sh_write = sh_write
        
    def run(self):
        '''
        burn-in phase.  Iterate <block_size> number of times and then report current sigma_delta
        '''
        self.seir_sampler.T = self.sh_temperatures[self.index]
        self.seir_sampler.sigma_delta = self.sh_sigma_delta[self.index]
        self.seir_sampler.write = self.sh_write[self.index]

        for i in range(self.blocks_to_burn):
            self.seir_sampler.run( self.block_size, self.thin, self.file_handles, True )
            self.sh_sigma_delta[self.index] = self.seir_sampler.sigma_delta
                    
            self.reusable_barrier.wait_phase_one()
            self.reusable_barrier.wait_phase_two()

        '''
        sampling phase.  Iterate <block_size> number of times
        '''    
        for i in range(self.num_blocks):
            '''
            1.
            main thread:
            draw random transposition (k,l)
            update sh_cross_temperatures (kth and lth positions transposed)
            '''
            self.reusable_barrier.wait_phase_one()
            
            '''
            2.
            main thread waits
            '''
            self.reusable_barrier.wait_phase_two()

            
            self.seir_sampler.T = self.sh_temperatures[self.index]
            self.seir_sampler.sigma_delta = self.sh_sigma_delta[self.index]
            self.seir_sampler.write = self.sh_write[self.index]

            self.seir_sampler.run( self.block_size, self.thin, self.file_handles, False )
            self.sh_likelihoods[self.index] = self.seir_sampler.log_likelihood_partial_integrated( self.seir_sampler.T )
            self.sh_cross_likelihoods[self.index] = self.seir_sampler.log_likelihood_partial_integrated( self.sh_cross_temperatures[self.index] )

           

            '''
            3.
            main thread:
            calculate acceptance probability
            accept / reject
            if accept:
            transpose sh_<>
            '''
            self.reusable_barrier.wait_phase_one()
            self.reusable_barrier.wait_phase_two()
          
    
            
