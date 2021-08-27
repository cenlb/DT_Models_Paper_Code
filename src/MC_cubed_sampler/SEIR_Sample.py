import numpy as np
import multiprocessing as m_proc
import time
import scipy.stats as stats
import sys

ranD_lock = m_proc.Lock()
file_lock = m_proc.Lock()

## @brief random uniform sample from union of collection of intervals
#
# @author https://stackoverflow.com/questions/21444597/random-number-in-an-union-of-intervals-in-python 30th May 2019
#
# @param[in] intervals Sequence of start,end tuples
def random_from_intervals(intervals, ranD): # intervals is a sequence of start,end tuples
    total_size = sum(end-start for start,end in intervals)

    ranD_lock.acquire()
    n = ranD.uniform(0.0, total_size)
    ranD_lock.release()

    for start, end in intervals:
        if n < end-start:
            return start + n
        n -= end-start

## @brief randomly selects an initial array @c E of exposure times
#
# @param[in] I,R numpy ndarrays of infection and removal times, index identifies HOST
#
# @param[in] init first exposure time
#
# @param[in] ranD instance of numpy.random.RandomState
#
# @warning (I,R) are assumed to be sorted with I in increasing order!
#
# @return list of exposure times, index identifies HOST, with first element @c init
def set_initial_E( I, R, init, ranD ):
    E = []
    intervals = []

    intervals.append( (I[0], R[0]) )
    E.append( I[0] - init )

    #print intervals[-1]
    #print E[-1]

    for i in range( 1, I.size):
        current_i, current_r = intervals[-1]
        if I[i] > current_r:
            intervals.append( (I[i], R[i]) )
            current_i, current_r = intervals[-1]
        elif R[i] > current_r:
            intervals.pop()
            intervals.append( ( current_i, R[i] ) )
            current_i, current_r = intervals[-1]
        intervals_prime = list(intervals)
        l, u = intervals_prime.pop()
        
        intervals_prime.append( (l, min(u,I[i])) )
        E.append( random_from_intervals( intervals_prime, ranD ) )
        #print intervals[-1]
        #print E[-1]
    #print intervals
    return E
    


## @brief helper function for SEIR likelihood
#
# @param V1, V2 Vectors of start and end times
#
# @param gam Exponential dist. rate, \f$\gamma\f$
# @usage Assumes V1 < V2
#
# @return \f$ \log{\prod_{j=1}^m \gamma e^{-\gamma (V^2_j - V^1_j)}} = m \log \gamma - \gamma \sum_{j=1}^m (V^2_j - V^1_j)\f$
def lifetime(V1,V2,gam):
    m = R.size
    return m*np.log(gam) - gam*np.sum(V2 - V1)



## @brief helper function for SEIR likelihood 
#
# @param E, I, R unsorted arrays of exposure, infection and removal times
#
# @warning Assumes: \f$m = |E| = |I| = |R| \f$, index identifies HOST
#
# @usage For S/E/I/R, K(E,I,R). For S/I/R, K(I,I,R) - note, \f$I_j\f$ does not fall in \f$(I_j,R_j)\f$
#
# @return \f$\log {\prod_{i \neq \kappa} I_{E_i-} }\f$, where \f$E_\kappa = \min_{E_1, \ldots ,E_m}\f$ and -1 when \f$(E,I,R)\f$ not legal - i.e. when some \f$ E_i \f$ does not fall within some \f$ (I_j, R_j) \f$
#
# @warning Always test for -1 first!
def K(E,I,R):
    kappa = np.argmin(E)
    E_short = np.delete(E,kappa) #Product is over all exposures EXCEPT the first, indexed as \kappa

    x1 , y1 = np.ix_(E_short,I)
    x1 , y2 = np.ix_(E_short,R)

    a = (x1 > y1)*(x1 < y2) #For each j (infection event) Does exposure event fall in interval (I_j,R_j) ?
    aSum = np.sum(a.astype(bool),axis=1) #For each exposure (not including \kappa) count number of (I_j,R_j) it falls between
    #print aSum

    if not np.all(aSum):
        return -1  #returns -1 if E,D,R not legal, i.e. no dead prior to an exposure - there is a zero (no I's at time of exposure)
    else:
        return np.sum(np.log(aSum))

## As K(E,I,R) but prints aSum for inspection
def K_display_sum(E,I,R):
    kappa = np.argmin(E)
    E_short = np.delete(E,kappa)

    x1 , y1 = np.ix_(E_short,I)
    x1 , y2 = np.ix_(E_short,R)

    a = (x1 >= y1)*(x1 < y2)
    aSum = np.sum(a.astype(bool),axis=1)
    print(aSum)

    if not np.all(aSum):
        return -1  #returns -1 if E,D,R not legal, i.e. no dead prior to an exposure
    else:
        return np.sum(np.log(aSum))

## @brief helper function: factor of SEIR likelihood 
#
# \f$S_t, I_t \f$ numbers of susceptible and infectious hosts
#
# @param E, I, R unsorted arrays of exposure, infection and removal times
#
# @param N Total host population size
#
# @usage For SEIR, A(E,I,R,N).  For SIR, A(I,I,R,N)
#
# @warning Assumes: \f$m = |E| = |I| = |R|\f$, index identifies HOST, \f$N \geq m \f$
# @return \f$ \int_{E_\kappa}^\infty S_t I_t \, dt\f$ 
def A(E,I,R,N):
    m = R.size
    A=(N - m) * np.sum(R - I)

    Ei, Rj = np.ix_(E, R)
    s1 = np.minimum(Ei, Rj)

    Ei, Ij = np.ix_(E, I)
    s2 = np.minimum(Ei, Ij)

    return A + np.sum(s1 - s2)

class SEIR_Sampler:
    ## @brief Constructor
    #
    def __init__(self, ranD, prior_nu, prior_lambda, I, R, N, parameter_inits, sigma_delta, T, write):

        self.nu_beta, self.nu_gamma, self.nu_delta = prior_nu
        self.lambda_beta, self.lambda_gamma, self.lambda_delta = prior_lambda

        self.beta, self.gamma, self.delta = parameter_inits
        self.E_hold = -1.0
        self.num_hold = -1.0
        self.like_hold = -1.0
        self.E_proposals_accepted = 0
        

        self.E = np.zeros(R.size)
        self.I = I
        self.R = R
        self.N = N        
        self.ranD = ranD
        self.T = T
        self.write = write
        self.sigma_delta = sigma_delta
        print('SEIR_Sampler constructor called')

    ## @brief Initialise E via particle filter using parameter_inits
    #
    # (however, currently using crude method - delete E0 input param)
    def initialise(self, E0):
        ind = np.argsort(self.I)
        self.E = set_initial_E( self.I[ind], self.R[ind], E0, self.ranD )
        
        

    def run(self, iterations, thin, file_handles, burning):
        if self.T != 1.0 and self.write:
            print('Something horrible as happend...\n', self.T, self.write)
            sys.exit()
        beta_out, delta_out, gamma_out, K_out, A_out, E_kappa_out = file_handles
        for i in range(iterations):
            self.delta_MH( self.T, burning )
            self.propose_E( self.T )

            if self.write and i % thin == 0 and not burning:
                file_lock.acquire()
                beta_out.write('%s\n' %self.beta )
                delta_out.write('%s\n' %self.delta )
                gamma_out.write('%s\n' %self.gamma )
                K_out.write('%s\n' %K( self.E, self.I, self.R ))
                A_out.write('%s\n' %A( self.E, self.I, self.R, self.N ))
                E_kappa_out.write('%s\n' %np.amin( self.E ))
                file_lock.release()
        
        beta_out.flush()
        delta_out.flush()
        gamma_out.flush()
        K_out.flush()
        A_out.flush()
        E_kappa_out.flush()

                
    
        
        
    ## @brief heated log-likelihood, \f$ \log \pi(\beta, \delta, \gamma \,|\,E,I,R)^{\frac{1}{T}} \f$ for exponentially-distributed E and I lifetimes
    #
    # @return \f$ \log \left (\prod_{i \neq \kappa} \{\beta I_{E_i -}\} e^{-\beta \int_{E_\kappa}^\infty S_t I_t\, dt}  \prod_{j=1}^m \delta e^{-\delta (I_j - E_j)} \prod_{j=1}^m \gamma e^{-\gamma (R_j - I_j)}\right)^{\frac{1}{T}}\f$
    #
    # @return \f$ = \frac{1}{T}\left(\mathtt{K}(E,I,R) - \beta \mathtt{A}(E,I,R,N) + \mathtt{lifetime}(E,I,\delta) + \mathtt{lifetime}(I,R,\gamma)\right)\f$ 
    def log_likelihood(self, T):
        m = self.R.size
        res = (m - 1) * np.log(self.beta) + K(self.E, self.I, self.R)
        res = res - self.beta * A(self.E,self.I,self.R,self.N)
        res = res + lifetime(self.E, self.I, self.delta)
        res = res + lifetime(self.I, self.R, self.gamma)
        return res / T


    ## @brief \f$ \log \int \pi(E,I,R \,|\,\beta, \delta, \gamma)\pi(\beta)\,\pi(\delta)\,\pi(\gamma)\,d\beta\, d\delta\, d\gamma \f$ . Independent gamma priors: (shape, rate) = (nu_, lambda_)
    #
    # Can lose Gamma(m+nu_beta -1) as constant since m is fixed by this problem (what does this mean???)
    #
    # @return \f$ \log \left( \prod_{i \neq \kappa} \{ I_{E_i -}\} \left( \int_{E_\kappa}^\infty S_t I_t\, dt + \lambda_\beta\right)^{-(m + \nu_\beta -1)}  \left(\sum_{j=1}^m (I_j-E_j) + \lambda_\delta\right)^{-(m + \nu_\delta)} \left(\sum_{j=1}^m (R_j-I_j) + \lambda_\gamma\right)^{- (m + \nu_\gamma)}\right)\f$
    #
    # @return \f$ = \mathtt{K}(E,I,R) - (m + \nu_\beta - 1)\log(\mathtt{A}(E,I,R,N) + \lambda_\beta)\f$
    #
    # @return \f$- (m + \nu_\delta)\log(\sum_{j=1}^m (I_j-E_j) + \lambda_\delta) - (m + \nu_\gamma)\log(\sum_{j=1}^m (R_j-I_j) + \lambda_\gamma)\f$
    def log_likelihood_integrated(self, T):
        m = self.R.size
        res = K(self.E, self.I, self.R)
        res = res - (m + self.nu_beta - 1) * np.log( A(self.E,self.I,self.R,self.N) + self.lambda_beta )
        res = res - (m + self.nu_delta) * np.log( np.sum(self.I - self.E) + self.lambda_delta )
        res = res - (m + self.nu_gamma) * np.log( np.sum(self.R - self.I) + self.lambda_gamma )
        return res / T
        '''
    ## @brief \f$ \log \int \pi(E,I,R \,|\,\beta, \delta, \gamma)\pi(\beta)\,\pi(\gamma)\,d\beta\, d\gamma + \log p(\delta)\f$ .
    def log_likelihood_partial_integrated(self, T):
        m = self.R.size
        res = K(self.E, self.I, self.R)
        res = res - (m + self.nu_beta - 1) * np.log( A(self.E,self.I,self.R,self.N) + self.lambda_beta )
        res = res - (m + self.nu_gamma) * np.log( np.sum(self.R - self.I) + self.lambda_gamma )
        res = res + m*np.log(self.delta) - self.delta * np.sum(self.I - self.E)
        res = res + (self.nu_delta - 1) * np.log( self.delta ) - self.lambda_delta * self.delta
        return res / T
        '''
    ## @brief \f$ \log \int \pi(E,I,R \,|\,\beta, \delta, \gamma)\pi(\beta)\,\pi(\gamma)\,d\beta\, d\gamma + \log p(\delta)\f$. \f$ p(\delta) \sim U(0,10.0) \f$
    #def log_likelihood_partial_integrated_uniform_delta_prior(self, T):
    def log_likelihood_partial_integrated(self, T):
        if self.delta < 0.0 or self.delta > 10.0:
            return -float('inf')
        else:
            m = self.R.size
            res = K(self.E, self.I, self.R)
            res = res - (m + self.nu_beta - 1) * np.log( A(self.E,self.I,self.R,self.N) + self.lambda_beta )
            res = res - (m + self.nu_gamma) * np.log( np.sum(self.R - self.I) + self.lambda_gamma )
            res = res + m*np.log(self.delta) - self.delta * np.sum(self.I - self.E)
            #res = res + (self.nu_delta - 1) * np.log( self.delta ) - self.lambda_delta * self.delta
            return res / T
        
    ## @brief Gibbs sample parameter \f$\beta\f$
    #
    # Full conditional dist. \f$ p(\,\beta \,|\, \mathbf{E}, \mathbf{I}, \mathbf{R} ,\delta, \gamma \,) \sim \Gamma(m + \nu_\beta -1,   \int_{E_\kappa}^\infty S_t I_t\, dt + \lambda_\beta) \f$ (shape, rate)
    def beta_gibbs( self, T ):
        m = self.R.size
        a = (m + self.nu_beta - 1) / T
        sc = T / (self.lambda_beta + A(self.E, self.I, self.R, self.N))
        ranD_lock.acquire()
        ret = self.ranD.gamma( a , scale = sc)
        ranD_lock.release()
        self.beta = ret

    ## @brief Gibbs sample parameter \f$\gamma\f$
    #
    # Full conditional dist. \f$ p(\,\gamma \,|\, \mathbf{E}, \mathbf{I}, \mathbf{R} ,\beta, \delta \,) \sim \Gamma(m + \nu_\gamma, \sum_{j=1}^m (R_j - I_j) + \lambda_\gamma) \f$ (shape, rate)
    def gamma_gibbs( self, T ):
        m = self.R.size
        a = (m + self.nu_gamma) / T
        sc = T / (self.lambda_gamma + np.sum(self.R - self.I))
        ranD_lock.acquire()
        ret = self.ranD.gamma( a, scale = sc)
        ranD_lock.release()
        self.gamma = ret

    ## @brief Gibbs sample parameter \f$\delta\f$
    #
    # Full conditional dist. \f$ p(\,\delta \,|\, \mathbf{E}, \mathbf{I}, \mathbf{R} ,\beta, \gamma \,) \sim \Gamma(m + \nu_\delta, \sum_{j=1}^m (I_j - E_j) + \lambda_\delta) \f$ (shape, rate)
    def delta_gibbs( self, T ):
        m = self.R.size
        a = (m + self.nu_delta) / T
        sc = T / (self.lambda_delta + np.sum(self.I - self.E))
        ranD_lock.acquire()
        ret = self.ranD.gamma( a, scale = sc)
        ranD_lock.release()
        self.delta = ret

    def delta_MH( self, T , burning):
        self.num_hold = self.delta
        self.like_hold = self.log_likelihood_partial_integrated( T )

        ranD_lock.acquire()
        self.delta = self.delta + self.ranD.normal( scale = self.sigma_delta )
        ranD_lock.release()
        
        if self.delta <= 0:
            al = 0
        else:
            likeNew = self.log_likelihood_partial_integrated( T )        
            r = 0 if likeNew > self.like_hold else (likeNew - self.like_hold)/T
            al = np.exp(r)

        ranD_lock.acquire()
        u = self.ranD.uniform()
        ranD_lock.release()

        if u < al:
            if burning:
                self.sigma_delta = self.sigma_delta * 1.002

        else:
            self.delta = self.num_hold
            if burning:
                self.sigma_delta = self.sigma_delta * 0.999
        

    ## @brief Metropolis-Hastings step to propose and accept/reject change to one \f$E_i\f$ chosen at random
    #
    # Uses log_likelihood_integrated() rather than log_likelihood() for accept/reject
    #
    # \f$ E_j' = I_j - \tau \f$  with \f$ \tau \sim \exp(\delta)\f$ (rate)
    #
    # ACCEPT with probability \f$ \min\{1,\frac{\pi(E',I,R)}{\pi(E,I,R)}\frac{\delta e^{-\delta (I_j - E_j)}}{\delta e^{-\delta (I_j - E_j)}} \} = \min\{1,\frac{\pi(E',I,R)}{\pi(E,I,R)} e^{-\delta (E'_j - E_j)} \}\f$ 
    def propose_E( self, T ):
        m = self.R.size
        ranD_lock.acquire()
        j = self.ranD.randint(0, m)
        ranD_lock.release()
        self.like_hold = self.log_likelihood_partial_integrated( T )
        self.num_hold = self.E[j]

        ranD_lock.acquire()
        self.E[j] = self.I[j] - self.ranD.exponential( scale = 1 / self.delta )
        ranD_lock.release()

        fac = self.delta * (self.E[j] - self.num_hold)
        r = min(0, (self.log_likelihood_partial_integrated( T ) - self.like_hold) - fac)
        if K(self.E, self.I, self.R) < 0:
            al = 0
        else:
            al = np.exp(r)
            
        ranD_lock.acquire()    
        u = self.ranD.uniform()
        ranD_lock.release()
        
        if u < al:
            self.E_proposals_accepted = self.E_proposals_accepted + 1
        else:
            self.E[j] = self.num_hold

    def propose_E_uniform( self, T):
        m = self.R.size
        m = self.R.size
        ranD_lock.acquire()
        j = self.ranD.randint(0, m)
        ranD_lock.release()
        self.like_hold = self.log_likelihood_partial_integrated( T )
        self.num_hold = self.E[j]

        I__ = np.minimum( self.I, self.I[j] )
        R__ = np.minimum( self.R, self.I[j] )

        print('----------------------------')
        print(self.I[j], self.R[j])
        print('----------------------------')
        for i,r in zip( I__, R__ ):
            print(i,r)
        
            

def exchange( samplerA, samplerB):
    T = samplerA.T
    w = samplerA.write
    sig = samplerA.sigma_delta
    
    samplerA.T = samplerB.T
    samplerA.write = samplerB.write
    samplerA.sigma_delta = samplerB.sigma_delta

    samplerB.T = T
    samplerB.write = w
    samplerB.sigma_delta = sig

def random_transpose( N, ranD):
    ranD_lock.acquire()
    i = ranD.randint( 0, N )
    j = ranD.randint( 0, N - 1 )
    ranD_lock.release()
    j = j if j < i else j + 1
    return i,j



        
        
        
        


