import numpy as np
import pickle
import time
import SEIR_Simulate as sim
from scipy import stats
import sys

if len(sys.argv) != 9:
    print('<exec> path file_no alp bet gam eps rho N')
    sys.exit()

path = sys.argv[1]
file_no = sys.argv[2]
alp = float(sys.argv[3])
bet = float(sys.argv[4])
gam = float(sys.argv[5])
eps = float(sys.argv[6])
rho = float(sys.argv[7])
N = int(sys.argv[8])

s = int(time.time())
ranD = np.random.default_rng(s)

def E_(args):
    return 0.0

def I_(args):
    gamma, rD = args
    return rD.gamma( shape = 1.10, scale = 1/gamma )

mySimulator = sim.SEIRP_Simulate(ranD)

mySimulator.run(alp, bet, eps, rho, E_, (mySimulator.ranD), I_, (gam, mySimulator.ranD), N, 0.0)

with open(path + 'E_' + file_no, 'wb') as fE:
    with open(path + 'I_' + file_no, 'wb') as fI:
        with open(path + 'R_' + file_no,'wb') as fR:
            pickle.dump(mySimulator.E, fE)
            pickle.dump(mySimulator.I, fI)
            pickle.dump(mySimulator.R, fR)
            

if int(file_no) % 10 == 0:
    print(file_no)
