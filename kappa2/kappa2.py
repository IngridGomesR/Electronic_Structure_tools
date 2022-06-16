from scipy.spatial.transform import Rotation as R
import numpy as np
from lx import tools

weights = {
    '1': 1.00797,    '2': 4.00260,    '3': 6.941,    '4': 9.01218,    '5': 10.81,    '6': 12.011,    '7': 14.0067,    '8': 15.9994,
    '9': 18.998403,    '10': 20.179,    '11': 22.98977,    '12': 24.305,    '13': 26.98154,    '14': 28.0855,    '15': 30.97376,
    '16': 32.06,    '17': 35.453,    '19': 39.0983,    '18': 39.948,    '20': 40.08,    '21': 44.9559,    '22': 47.90,    '23': 50.9415,
    '24': 51.996,    '25': 54.9380,    '26': 55.847,    '28': 58.70,    '27': 58.9332,    '29': 63.546,    '30': 65.38,    '31': 69.72,
    '32': 72.59,    '33': 74.9216,    '34': 78.96,    '35': 79.904,    '36': 83.80,    '37': 85.4678,    '38': 87.62,    '39': 88.9059,
    '40': 91.22,    '41': 92.9064,    '42': 95.94,    '43': 98,    '44': 101.07,    '45': 102.9055,    '46': 106.4,    '47': 107.868,
    '48': 112.41,    '49': 114.82,    '50': 118.69,    '51': 121.75,    '53': 126.9045,    '52': 127.60,    '54': 131.30,    '55': 132.9054,
    '56': 137.33,    '57': 138.9055,    '58': 140.12,    '59': 140.9077,    '60': 144.24,    '61': 145,    '62': 150.4,    '63': 151.96,
    '64': 157.25,    '65': 158.9254,    '66': 162.50,    '67': 164.9304,    '68': 167.26,    '69': 168.9342,    '70': 173.04,    '71': 174.967,
    '72': 178.49,    '73': 180.9479,    '74': 183.85,    '75': 186.207,    '76': 190.2,    '77': 192.22,    '78': 195.09,    '79': 196.9665,
    '80': 200.59,    '81': 204.37,    '82': 207.2,    '83': 208.9804,    '84': 209,    '85': 210,    '86': 222,    '87': 223,    '88': 226.0254,
    '89': 227.0278,    '91': 231.0359,    '90': 232.0381,    '93': 237.0482,    '92': 238.029,    '94': 242,    '95': 243,    '97': 247,
    '96': 247,    '102': 250,    '98': 251,    '99': 252,    '108': 255,    '109': 256,    '100': 257,    '101': 258,    '103': 260,
    '104': 261,    '107': 262,    '105': 262,    '106': 263,    '110': 269,    '111': 272,    '112': 277,
}

mass = {'H' : 1.008,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
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
    return dip_vec, mu 

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
#def trans_dim(fileDim): 
    #vec_dim, atom_dim  = tools.pega_geom(fileDim)
    #n_atoms = np.shape(atomos)[0]
    #for i in vec_dim:
    #    if vec_dim[i+1] - vec_dim[i] > 4:

    #mol1    = vetores[:int(n_atoms/2),:]
    #mol2    = vetores[int((n_atoms/2)):,:]
    #atomos1 = atomos[:int(n_atoms/2),:]
    #atomos2 = atomos[int((n_atoms/2)):,:]
    #avaliar de o atomo é numero ou letra, se numero usar weights
    #massas1 = np.array([weights.get(str(i[0])) for i in atomos1])
    #massas2 = np.array([weights.get(str(i[0])) for i in atomos2])
    #cm1     = [np.sum(massas1*mol1[:,0].flatten()),np.sum(massas1*mol1[:,1].flatten()),np.sum(massas1*mol1[:,2].flatten())]/np.sum(massas1)
    #cm2     = [np.sum(massas2*mol2[:,0].flatten()),np.sum(massas2*mol2[:,1].flatten()),np.sum(massas2*mol2[:,2].flatten())]/np.sum(massas2)
    #trans = cm2 - cm1
    #dcm = np.sqrt(np.inner(trans,trans))*0.1 #nm
    #return trans, dcm


#def old_kappa2(fileS1,fileDim,atomos_opt,atomos_dim1,atomos_dim2):
    #vetores = tools.pega_geom(fileS1)
    #matriz  = gera_base(vetores,atomos_opt)
    #inversa = np.linalg.inv(matriz)
    #dip_vec, _ = get_moment(fileS1)
    #print(np.sqrt(np.inner(dip_vec,dip_vec)))
    #nmu     = np.matmul(dip_vec,inversa)
    #_ , vetores = coords(fileDim)
    #matriz_dim1 = gera_base(vetores,atomos_dim1)
    #mu_final1   = np.matmul(nmu,matriz_dim1)
    #print(np.sqrt(np.inner(mu_final1,mu_final1)))
    #atomos_dim2 = atomos_dim + numero de atomos da molécula 
    #matriz_dim2 = gera_base(vetores,atomos_dim2)
    #mu_final2   = np.matmul(nmu,matriz_dim2)
    #mu_d = mu_final1/np.sqrt(np.inner(mu_final1,mu_final1))
    #mu_a = mu_final2/np.sqrt(np.inner(mu_final2,mu_final2))
    #trans, dcm  = trans_dim(fileDim)
    #print(dcm)
    #r = trans/np.sqrt(np.inner(trans,trans))
    #kappa  = np.inner(mu_d,mu_a) - 3*np.inner(r,mu_d)*np.inner(r,mu_a)
    #kappa2 = kappa**2
    #return kappa2

#old_kappa2('dim1S1.log','Dimero.xyz',[27,8,26],[27,8,26],[63,44,62])

def translation(atomos1, vetores1, atomos2, vetores2):
    #avaliar se atomos é numero ou letra
    massas1 = np.array([weights.get(str(i[0])) for i in atomos1])
    massas2 = np.array([weights.get(str(i[0])) for i in atomos2])
    cm1     = [np.sum(massas1*vetores1[:,0].flatten()),np.sum(massas1*vetores1[:,1].flatten()),np.sum(massas1*vetores1[:,2].flatten())]/np.sum(massas1)
    cm2     = [np.sum(massas2*vetores2[:,0].flatten()),np.sum(massas2*vetores2[:,1].flatten()),np.sum(massas2*vetores2[:,2].flatten())]/np.sum(massas2)
    trans = cm2 - cm1
    dcm = np.sqrt(np.inner(trans,trans))*0.1 #nm
    return trans, dcm

def kappa2(file1S1, file2S1, fileDim):
    vec_opt1, atom_opt1 = tools.pega_geom(file1S1)
    vec_opt2, atom_opt2 = tools.pega_geom(file2S1)

    vec_dim, atom_dim  = tools.pega_geom(fileDim) 
    n_atoms = np.shape(atom_dim)[0] #Pega o numero de atomos total
    atom_dim1 = atom_dim[:int(n_atoms/2)] #Pega a primeira metade
    atom_dim2 = atom_dim[int((n_atoms/2)):] #Pega a segunda metade
    vec_dim1    = vec_dim[:int(n_atoms/2),:]
    vec_dim2    = vec_dim[int((n_atoms/2)):,:]
    
    dip_vec1, mu1 = get_moment(file1S1)
    dip_vec2, mu2 = get_moment(file2S1)

    trans_dim1, dcm_dim1 = translation(atom_opt1, vec_opt1, atom_dim1, vec_dim1)
    dim1_trans = vec_dim1 - trans_dim1 #translada o dim pro centro de massa do opt
    r1 =  R.align_vectors(dim1_trans, vec_opt1) #matriz que alinha o opt na posição do dim

    trans_dim2, dcm_dim2 = translation(atom_opt2, vec_opt2, atom_dim2, vec_dim2)
    dim2_trans = vec_dim2 - trans_dim2
    r2 =  R.align_vectors(dim2_trans, vec_opt2)

    mu_final1 = r1[0].apply(dip_vec1) 
    mu_final2 = r2[0].apply(dip_vec2) 
    
    mu_d = mu_final1/np.sqrt(np.inner(mu_final1,mu_final1))
    mu_a = mu_final2/np.sqrt(np.inner(mu_final2,mu_final2))
        
    trans = trans_dim2 -  trans_dim1
    r = trans/np.sqrt(np.inner(trans,trans))
    kappa  = np.inner(mu_d,mu_a) - 3*np.inner(r,mu_d)*np.inner(r,mu_a)
    kappa2 = kappa**2
    return kappa2


def run_kappa():
    file1S1 = tools.fetch_file("log",['.log'])#input('log file from mol1: ')
    file2S1 = tools.fetch_file("log",['.log'])#input('log file from mol2: ')
    fileDim = tools.fetch_file("dimer",['.log','.xyz'])#input('xyz file from dimer: ')
    kappa = kappa2(file1S1, file2S1, fileDim)
    print("Orientation Factor (k^2): {0:.2f}".format(kappa))

run_kappa()
