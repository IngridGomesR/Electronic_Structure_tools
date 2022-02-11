import numpy as np

massas={'H' : 1.008,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
                 'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
                 'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,\
                 'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,\
                 'Cl' : 35.453, 'AR' : 39.948, 'K' : 39.098, 'Ca' : 40.078,\
                 'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,\
                 'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,\
                 'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,\
                 'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,\
                 'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,\
                 'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,\
                 'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,\
                 'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,\
                 'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,\
                 'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,\
                 'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'Gd' : 157.25,\
                 'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,\
                 'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,\
                 'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,\
                 'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,\
                 'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,\
                 'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,\
                 'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,\
                 'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247, 'Bk' : 247,\
                 'Ct' : 251, 'Es' : 252, 'Fm' : 257, 'Md' : 258, 'No' : 259,\
                 'Lr' : 262, 'Rf' : 261, 'Db' : 262, 'Sg' : 266, 'Bh' : 264,\
                 'Hs' : 269, 'Mt' : 268, 'Ds' : 271, 'Rg' : 272, 'Cn' : 285,\
                 'Nh' : 284, 'Fl' : 289, 'Mc' : 288, 'Lv' : 292, 'Ts' : 294,\
                 'Og' : 294}

# PEGA O MÓDULO E O VETOR MOMENTO DE DIPOLO DA MOLÉCULA
def get_moment(fileS1):
    busca = "transition electric dipole moments" 
    n = -1
    with open(fileS1 ,'r') as f:
        for line in f:
            if busca in line:
                n = 0
                dip_sqrd = []
            elif n >= 0 and n < 1:
                n += 1
            elif n == 1 and "transition velocity" not in line:    
                line = line.split()
                dip_vec = []
                for j in range(1,4):
                    dip_vec.append(float(line[j]))
                dip_sqrd.append(float(line[4]))
                n +=1
            elif n >= 3:
                n = -1
        mu=np.sqrt(float(dip_sqrd[-1])) 
    return mu, dip_vec

def gera_base(vetores,atomos):
    a = vetores[atomos[1]-1,:] - vetores[atomos[0]-1,:]
    b = vetores[atomos[2]-1,:] - vetores[atomos[0]-1,:]
    c = np.cross(a,b)
    a = a/np.sqrt(np.inner(a,a))
    b = b/np.sqrt(np.inner(b,b))
    c = c/np.sqrt(np.inner(c,c))
    matriz = np.vstack((a,b))
    matriz = np.vstack((matriz,c))
    return matriz

# PEGA OS ATOMOS E VETORES DOS DÍMEROS
def coords(fileDim):
    vetores = np.array([0,0,0])
    atomos  = np.array([0])
    with open(fileDim, 'r') as f:
        for line in f:
            line = line.split()
            if len(line) < 4:
                pass
            else:
                vetor   = np.array([float(line[1]),float(line[2]),float(line[3])])
                vetores = np.vstack((vetores,vetor))
                atomo   = np.array([line[0]])
                atomos  = np.vstack((atomos,atomo))
    atomos  = np.delete(atomos,0,0)
    vetores = np.delete(vetores,0,0)
    return atomos, vetores

# CALCULA O VETOR TRANSLAÇÃO E A DISTÂNCIA DE CENTRO DE MASSAS ENTRE OS DÍMEROS
def trans_dim(fileDim): 
    atomos, vetores = coords(fileDim)
    n_atoms = np.shape(atomos)[0]
    mol1    = vetores[:int(n_atoms/2),:]
    mol2    = vetores[int((n_atoms/2)):,:]
    atomos1 = atomos[:int(n_atoms/2),:]
    atomos2 = atomos[int((n_atoms/2)):,:]
    massas1 = np.array([massas.get(str(i[0])) for i in atomos1])
    massas2 = np.array([massas.get(str(i[0])) for i in atomos2])
    cm1     = [np.sum(massas1*mol1[:,0].flatten()),np.sum(massas1*mol1[:,1].flatten()),np.sum(massas1*mol1[:,2].flatten())]/np.sum(massas1)
    cm2     = [np.sum(massas2*mol2[:,0].flatten()),np.sum(massas2*mol2[:,1].flatten()),np.sum(massas2*mol2[:,2].flatten())]/np.sum(massas2)
    trans = cm2 - cm1
    dcm = np.sqrt(np.inner(trans,trans))*0.1 #nm
    return trans, dcm

def pega_geom(file):
    vetores = np.array([0,0,0])
    atomos  = np.array([0])
    arquivo = open(file, 'r')
    lista = arquivo.readlines()
    i=1
    for linha in lista:
        i += 1
        if 'Standard orientation:' in linha:
            j=i+3
        if "NAtoms=" in linha:
            aux=linha.split()
            num_a=int(aux[1])
    
    geom=lista[j:j+num_a]
    for linha in geom:
        aux=linha.split()
        vetor   = np.array([float(aux[3]),float(aux[4]),float(aux[5])])
        vetores = np.vstack((vetores,vetor))
        atomo   = np.array([aux[1]])
        atomos  = np.vstack((atomos,atomo))
    atomos  = np.delete(atomos,0,0)
    vetores = np.delete(vetores,0,0)    
    return atomos, vetores

def kappa2(fileS1,fileDim,atomos_opt,atomos_dim1,atomos_dim2):
    _, vetores = pega_geom(fileS1)
    matriz  = gera_base(vetores,atomos_opt)
    inversa = np.linalg.inv(matriz)
    _, dip_vec = get_moment(fileS1)
    nmu     = np.matmul(dip_vec,inversa)
    _ , vetores = coords(fileDim)
    matriz_dim1 = gera_base(vetores,atomos_dim1)
    mu_final1   = np.matmul(nmu,matriz_dim1)
    matriz_dim2 = gera_base(vetores,atomos_dim2)
    mu_final2   = np.matmul(nmu,matriz_dim2)
    mu_d = mu_final1/np.sqrt(np.inner(mu_final1,mu_final1))
    mu_a = mu_final2/np.sqrt(np.inner(mu_final2,mu_final2))
    trans, _  = trans_dim(fileDim)
    r = trans/np.sqrt(np.inner(trans,trans))
    kappa  = np.inner(mu_d,mu_a) - 3*np.inner(r,mu_d)*np.inner(r,mu_a)
    kappa2 = kappa**2
    return kappa2

#kappa2 = kappa2('3bS1.log','Dimero3b1.xyz',[9,18,16],[32,15,1],[92,75,61])
kappa2 = kappa2('3bS1.log','Dimero3b2.xyz',[9,18,16],[32,15,1],[92,75,61])
#kappa2 = kappa2('3bS1.log','Dimero3b3.xyz',[9,18,16],[2,45,31],[62,105,91])
print(kappa2)
