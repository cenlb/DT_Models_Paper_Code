import numpy as np
import time
import scipy.stats as stats
import sys

## @brief Sample from piecewise constant hazard function, \f$ h(t) \f$
#
# @param left_end_points ndarray sorted in increasing order of left end points \f$ t_0, t_1, \ldots ,t_p \f$ (right end points of intervals are \f$ t_1, \ldots, t_p, +\infty \f$)
#
# @param hazard numpy ndarray of values \f$ h(t) \f$ for \f$ t \in [t_i, t_{i+1} )\f$
#
# @param ranD instance of numpy.random.RandomState
#
# @usage hazard and left_end_points assumed to have equal length
#
# @return random sample of R.V. with hazard function \f$h\f$

def sample_from_step_hazard(hazard, left_end_points, ranD):
    res = float('inf')

    r_end_points = left_end_points[1::]
    r_end_points = np.hstack( (r_end_points, [float('inf')]) )

    for h,t0,t1 in zip(hazard, left_end_points, r_end_points):
        #print h, t0, t1
        
        if h == 0:
            continue
        tau = ranD.exponential(scale = 1/h)
        if tau > t1 - t0:
            continue
        else:
            res = t0 + tau
            return res
        
    return res

## @brief Helper function to convert infection and removal times to a step function
#
# @param[in] I,R unsorted numpy ndarrays, of EQUAL LENGTH, of INFECTION and REMOVAL times.  Indexed by HOST.
#
# @param[in] current_time Before this, hazard is zero
#
# @param[out] infected, left_end_points format of input to sample_from_step_hazard
def convert_to_step_function(I, R, current_time):
    infected = np.hstack((np.ones(I.size),-np.ones(R.size)))

    tx = np.hstack((I, R))
    inds = np.argsort(tx)

    tx = tx[inds]
    infected = np.cumsum(infected[inds])

    #now need pop. size at current_time
    ax = np.argwhere(tx <= current_time)
    if ax.size == 0:
        current_infected = 0
    else:
        i = np.amax(ax)
        current_infected = infected[i]

    #not interested in anything before current_time
    infected = infected[tx > current_time]
    tx = tx[tx > current_time]
    
    infected = np.hstack(( [current_infected], infected ))

    left_end_points = np.hstack(([current_time],tx))
    
    return infected, left_end_points
    

## @brief helper function for sample_indirect_exposure()
#
# @param[in] ranD instance of numpy.random.RandomState
#
# @param[in] c multiplier of hazard function
#
# @param[out] p, t size of pathogen population \f$ P_t \f$ and time of sampled exposure or \f$+\infty\f$ if no exposure occurs
#
# @return a random sample from hazard function \f$ h(t) = c P_t\f$ for \f$ t \in [t_0, t_1) \f$, zero otherwise. \f$ P_t\f$ is immigration-death process \f$ P_{t_0} = p_0\f$ with immigration and death rates \f$\epsilon\f$ and \f$\rho\f$, respectively.
def sample_indirect_exposure_uniform(epsilon, rho, t0, t1, c, p0, ranD):
    p = p0
    t = t0

    while(True):
        #print 'current t,p = ', t, p

        if t >= t1:
            return p, float('inf')

        q = epsilon + rho * p

        if q == 0:
            return p, float('inf')

        dt = ranD.exponential( scale = 1/q )
        t_prime = min(t1, t + dt)

        v = ranD.uniform() * q
        if v < epsilon:
            p_prime = p + 1
        else:
            p_prime = p - 1

        #print 't_prime, p_prime = ',t_prime, p_prime

        if p == 0:
            tau = float('inf')
        else:
            tau = ranD.exponential( scale = 1 /( c * p ) )
            
        #print t + tau
        
        if t + tau >= t_prime:
            t = t_prime
            p = p_prime
            continue

        else:
            return p, t + tau

## @brief steps through stages of emission rate step function (usually \f$ \epsilon I_t \f$
#
# @param[in] emission_rates, left_end_points step function giving emission rates
#
# @param[in] rho, p0, c as in sample_indirect_exposure_uniform
#
# @param[in] break_at not current used
#
# @param[in] ranD instance of numpy.random.RandomState
#
# @return \f$ P_{t'}, t' \f$, where \f$t'\f$ is time of indirect exposure sampled and inf if no exposure occurs

def sample_indirect_exposure(emission_rates, left_end_points, rho, p0, c, break_at, ranD):
    ret = float('inf')
    p = p0
    right_end_points = left_end_points[1::]
    right_end_points = np.hstack( ( right_end_points, [float('inf')] ))
    #for epsilon, t0, t1 in zip( emission_rates,  left_end_points, right_end_points ):
        #print epsilon, t0, t1
    for epsilon, t0, t1 in zip( emission_rates,  left_end_points, right_end_points ):
        #print 'epsilon, rho, t0, t1, p = ', epsilon, rho, t0, t1, p
        p, t = sample_indirect_exposure_uniform(epsilon, rho, t0, t1, c, p, ranD)
        if t == float('inf'):
            #print epsilon, t, p
            continue
        else:
            return p, t

    return p, float('inf')


## @brief Simulates sets of host event times arising from S(E)IR(-P) processes
class SEIRP_Simulate:
    ## Constructor
    #
    # @param ranD instance of numpy.RandomState
    def __init__(self, ranD):
        self.E = np.array([])
        self.I = np.array([])
        self.R = np.array([])
        self.ranD = ranD
       
    ## @brief Simulate host event times.
    #
    # @param[in] alpha, beta, epsilon, rho indirect and direct transmission rates, pathogen emission and pathogen decay rates
    #
    # @param[in] E_sampler, I_sampler function parameters to sample from lifetime distributions for durations in E and I states
    #
    # @param[in] E_sampler_params, I_sampler_params tuples of auxillary parameters to be unpacked in _sampler bodies (should include self.ranD)
    #
    # @param[in] N, t0 total host population size and time of first exposure
    
    def run(self, alpha, beta, epsilon, rho, E_sampler, E_sampler_params, I_sampler, I_sampler_params, N, t0):
        t = t0
        nS = N-1
        nP = 0

        self.E = np.append( self.E, [t] )
        tau1 = E_sampler(E_sampler_params)
        tau2 = I_sampler(I_sampler_params)
        
        self.I = np.append( self.I, [t + tau1] )
        self.R = np.append( self.R, [t + tau1 + tau2] )

        while(True):
            #print 'top of loop. Susceptibles remaining:', nS
            #sys.stdout.flush()
            if nS <= 0:
                break
            t = self.E[-1]

            infected, l_end_points = convert_to_step_function(self.I, self.R, t)

            t_direct = float('inf')
            if beta > 0.0:
                t_direct = sample_from_step_hazard(beta * nS * infected, l_end_points, self.ranD)
            
            t_indirect = float('inf')
            if alpha > 0.0:
                emission_rates = epsilon * infected
                nP, t_indirect = sample_indirect_exposure(emission_rates, l_end_points, rho, nP, nS * alpha, 0.0, self.ranD)

            t = min( t_direct, t_indirect )

            if t == float('inf'):
                break

            else:
                self.E = np.append( self.E, [t] )
                tau1 = E_sampler(E_sampler_params)
                tau2 = I_sampler(I_sampler_params)
        

                self.I = np.append(self.I, [t + tau1])
                self.R = np.append(self.R, [t + tau1 + tau2])

                nS = nS - 1

        r = np.argsort(self.R)
        self.E = self.E[r]
        self.I = self.I[r]
        self.R = self.R[r]

    def reset(self):
        self.E = np.array([])
        self.I = np.array([])
        self.R = np.array([])

        
                
## @brief Get host susceptible, exposed, infected and removed host numbers at given times
#
# @param N total host population size
#
# @param E, I, R unsorted numpy ndarrays of exposure, infection and removal times
#
# @param t numpy ndarray of times
#
# @warning Assumes: \f$m = |E| = |I| = |R|\f$, index identifies HOST, \f$N \geq m \f$
def hosts_in_state(N, E, I, R, t ):
    tx, Ey = np.ix_( t, E )
    _E = (Ey <= tx)

    tx, Iy = np.ix_( t, I )
    _I = (Iy <= tx)

    tx, Ry = np.ix_( t, R )
    _R = (Ry <= tx)

    s = N - np.sum( _E , axis = 1 )
    e = np.sum( np.logical_and( _E, np.logical_not(_I) ), axis = 1 )
    i = np.sum( np.logical_and( _I, np.logical_not(_R) ), axis = 1 )
    r = np.sum( _R, axis = 1 )

    return s,e,i,r


## @brief How many individuals were infected
def final_outbreak_size(E,I,R):
    return R.size

## @brief Time and size of peak outbreak (when number of infected peaks)
def time_and_size_of_outbreak_peak(E,I,R):
    #sort by I
    ind = np.argsort(I)
    E = E[ind]
    I = I[ind]
    R = R[ind]
    _s,_e,_i,_r = hosts_in_state( R.size, E, I, R, I )

    #time and size of outbreak at its peak
    ind = np.argmax(_i)
    return I[ind]-I[0], _i[ind]
    


'''
def simulate_pathogen_uniform(epsilon ,rho, t0, t1, P0, ranD, break_at):    
        self.pathogen = np.append(self.pathogen,[P0])
        self.pathogen_time = np.append(self.pathogen_time,[t0])

        p = P0
        t = t0

        while(True):
            q = epsilon + (rho * p)
            if(q == 0):
                return

            #if(q < 0):
            #    print epsilon, rho, p, q
            #    sys.exit()

            dt = self.ranD.exponential( scale = 1/q )    
            t = t + dt

            if t > t1:
                return
            elif t0 <= break_at and break_at < t1 and t > break_at:
                return
            else:    
                v = self.ranD.uniform() * q
                if(v <= epsilon):            
                    p = p + 1
                else:
                    p = p - 1

                self.pathogen = np.append(self.pathogen,[p])
                self.pathogen_time = np.append(self.pathogen_time,[t])
'''
