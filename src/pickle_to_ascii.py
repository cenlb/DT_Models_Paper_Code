import numpy as np
import sys

path = sys.argv[1]

E = np.load(open(path + 'E_', 'rb'))
I = np.load(open(path + 'I_', 'rb'))
R = np.load(open(path + 'R_', 'rb'))


outfile = open(path + 'EIR.txt', 'w')
outfile.write( 'E\tI\tR\n' )
for e, i, r in zip(E, I, R):
    outfile.write( '%s\t%s\t%s\n' %(e, i, r) )

outfile.close()

                   


