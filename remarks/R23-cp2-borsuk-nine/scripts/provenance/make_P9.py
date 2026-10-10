# PROVENANCE ONLY: writes ../P9.txt (exact hex floats) from audit_capcover_best_4.npy. Not run by run_all.sh.
import numpy as np
P = np.load('audit_capcover_best_4.npy')
with open('../P9.txt', 'w') as f:
    f.write('# 9 points of C^3 (row i = p_i), each entry: Re Im as exact IEEE-754 hex floats\n')
    for row in P:
        f.write(' '.join(f'{z.real.hex()} {z.imag.hex()}' for z in row) + '\n')
