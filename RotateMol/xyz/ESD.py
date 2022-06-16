import numpy as np

def to_float(num):
    try:
        num = float(num)
    except:
        num = -np.inf
    return num 

def pega_homo(file):
    HOMOS = []
    with open(file, 'r') as f:
        for line in f:
            if ('OPT' in file and "Optimized Parameters" in line) :
                HOMOS = [] 
            if "occ. eigenvalues" in line:
                line = line.split()
                homos = line[4:]
                HOMOS.extend(homos)
    if len(HOMOS) == 0:
        with open(file, 'r') as f:
            for line in f:
                if "occ. eigenvalues" in line:
                    line = line.split()
                    homos = line[4:]
                    HOMOS.extend(homos)
    HOMOS = list(map(to_float,HOMOS))
    return HOMOS

def pega_lumo(file):
    LUMOS = []
    with open(file, 'r') as f:
        for line in f:
            if ('OPT' in file and "Optimized Parameters" in line) :
                LUMOS = [] 
            if "virt. eigenvalues" in line:
                line = line.split()
                lumos = line[4:]
                LUMOS.extend(lumos)
    if len(LUMOS) == 0:
        with open(file, 'r') as f:
            for line in f:
                if "virt. eigenvalues" in line:
                    line = line.split()
                    lumos = line[4:]
                    LUMOS.extend(lumos)
    LUMOS = list(map(to_float,LUMOS))
    return LUMOS


# Energy Splitting in Dimer
def ESD(file):
    HOMO = pega_homo(file) 
    LUMO = pega_lumo(file)
    homo1 = HOMO[-2]
    homo = HOMO[-1]
    lumo1 = LUMO[1]
    lumo = LUMO[0]
    print(homo1,homo,lumo1,lumo)
    t_e = (lumo1 - lumo)/2
    print("t_e: {0:.5f}".format(t_e*27.2114)) #eV
    t_b = (homo - homo1)/2
    print("t_b: {0:.5f}".format(t_b*27.2114)) #eV

#print('Para o dímero com distância 2A:')
#ESD('BenzDist2.log')

#print('Para o dímero com distância 4A:')
#ESD('BenzDist4.log')

#print('Para o dímero com distância 6A:')
#ESD('BenzDist6.log')

#print('Para o dímero com distância 4A e rotação de 12 graus:')
#ESD('benz1.log')

#print('Para o dímero com distância 4A e rotação de 24 graus:')
#ESD('benz2.log')

#print('Para o dímero com distância 4A e rotação de 36 graus:')
#ESD('benz3.log')

#print('Para o dímero com distância 4A e rotação de 48 graus:')
#ESD('benz4.log')

#print('Para o dímero virado 90º:')
#ESD('Benzeno90.log')

print('Para o dímero paralelo:')
ESD('Naphthalene.log')

print('Para o dímero virado 90º em x:')
ESD('Naphthalene90x.log')

print('Para o dímero virado 90º em y:')
ESD('Naphthalene90y.log')

print('Para o dímero virado 90º em z:')
ESD('Naphthalene90z.log')