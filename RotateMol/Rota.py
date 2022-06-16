from scipy.spatial.transform import Rotation as R
import numpy as np


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

def gerafile(namefile, atomos1, atomos2, vetores1, vetores2):
    with open(namefile,'w') as f:
        f.write(str(np.shape(vetores1)[0]+np.shape(vetores2)[0])+"\n"+"\n")
        for i in range(len(atomos1)):
            f.write(str(atomos1[i][0])+"\t"+str(vetores1[i,0])+"\t"+str(vetores1[i,1])+"\t"+str(vetores1[i,2])+"\n")
        for i in range(len(atomos2)):
            f.write(str(atomos2[i][0])+"\t"+str(vetores2[i,0])+"\t"+str(vetores2[i,1])+"\t"+str(vetores2[i,2])+"\n")
    return 


def centro(atomos1, vetores1):
    massas = np.array([mass.get(str(i[0])) for i in atomos1])
    cm     = [np.sum(massas*vetores1[:,0].flatten()),np.sum(massas*vetores1[:,1].flatten()),np.sum(massas*vetores1[:,2].flatten())]/np.sum(massas)
    return cm


def rota(file, trans, thetax, thetay, thetaz, newfile):
    atomos1, vetores1 = coords(file)
    r = R.from_rotvec([thetax, thetay, thetaz])
    vecrota = r.apply(vetores1)
    cm1 = centro(atomos1, vetores1)
    cm2 = centro(atomos1, vecrota)
    cm = cm2 - cm1
    vecrota = vecrota - cm
    vectrans = vecrota + trans
    gerafile(newfile, atomos1, atomos1, vetores1, vectrans)
    

rota('naphthalene.xyz', np.array([0,0,4]), 0, np.pi/2, 0, 'Naphthalene90y.xyz')




    
        