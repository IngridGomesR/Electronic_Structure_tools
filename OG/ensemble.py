import os
from subprocess import call
import subprocess
import random
import numpy as np

sistema = 4
grupos = np.round(np.linspace(0.017,0.28,10),2) 
total = 1

for n in range(total):
    for group in grupos:
    #group = random.choice(grupos)
        subprocess.call(['python','gen.py',str(sistema), str(group), str(n)])
