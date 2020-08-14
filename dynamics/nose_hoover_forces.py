#! /usr/bin/python


## SIMULATION OF A SYSTEM AT CONSTANT TEMPERATURE

##### PLEASE READ ALL THE COMMENTS!!!!


#### REMEMBER THAT IN PYTHON INDEX BEGINS IN 0 !!!
 
###################################################################
## If initial forces are set to zero because the initial geometry # 
## is a stationary point                                          #
##       !! IT IS DONE IN THE INSTANTATION !!                     #
##################################################################


import Mytools as my
import numpy as np
import os
import random as rnd
from math import *
import sys

__metaclass__= type

############### CLASS DEFINITION ###########################

class Hessian_problem(Exception):
    pass

class Point:
    """ Represent a point in a space of three dimensions.
    attributes: x,y,z."""

    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return'%g, %g, %g' % (self.x, self.y, self.z)

class Atoms:
    """Represent an atom in a molecule.
    attributes: mass, cart, vel, force"""

    def __init__(self, mass = 0.0):
        self.mass = mass
        self.cart = Point()
        self.vel  = Point()
        self.force = Point()

    def momenta(self,mv):
        t1 = 0.5 * mv[0] * mv[0] /self.mass
        t2 = 0.5 * mv[1] * mv[1] /self.mass
        t3 = 0.5 * mv[2] * mv[2] /self.mass
        return t1+t2+t3
  
    def scaling_A(self,factor,mv):
        self.vel.x += factor*mv[0]/self.mass
        self.vel.y += factor*mv[1]/self.mass
        self.vel.z += factor*mv[2]/self.mass
  
    def scaling_B(self,factor):  
        self.vel.x = factor*self.vel.x
        self.vel.y = factor*self.vel.y
        self.vel.z = factor*self.vel.z

    def gradiente(self,vector):
        self.force.x = -vector[0]
        self.force.y = -vector[1]
        self.force.z = -vector[2]
   
    def cinetica(self):
        k = 0.5*my.cuad(self.vel,self.vel)*self.mass
        return k

class Molecule:
    """Contains all the internal and cartesian coordinates of the atoms.
    attributes: atoms, bonds, angles, dihedrals, b0, a0, d0"""


    def __init__(self,all_atoms_cart = []):
        self.cart = all_atoms_cart    

    def internas_bonds(self,conex,nbond):
        self.bonds = [my.enlaces(self.cart[conex[i][0]-1], self.cart[conex[i][1]-1]) \
        for i in range(nbond)]		
	
    def internas_ang(self,conex,nang):	
        self.angles = [my.angulos(self.cart[conex[i][0]-1], self.cart[conex[i][1]-1], \
        self.cart[conex[i][2]-1]) for i in range(nang)]		

    def internas_dihe(self,conex,ndih):	
        try:
            self.qt_1 = self.dihedrals 
            self.dihedrals = [my.dihedros(self.cart[conex[i][0]-1], self.cart[conex[i][1]-1], \
            self.cart[conex[i][2]-1], self.cart[conex[i][3]-1],self.qt_1[i]) for i in range(ndih)]
        except AttributeError:
            self.qt_1 = [0.0 for w in range(ndih)]
            self.dihedrals = [my.dihedros(self.cart[conex[i][0]-1], self.cart[conex[i][1]-1], \
            self.cart[conex[i][2]-1], self.cart[conex[i][3]-1],self.qt_1[i]) for i in range(ndih)]

    def potential(self,mtx_transp,grad_int,mtx_Hess):
        qo = self.b0 + self.a0 + self.d0
        q1 = self.bonds + self.angles + self.dihedrals
        Q = np.array(([q1[i] - qo[i] for i in range(len(qo))]), dtype = float)
        
        pot = grad_int+np.dot(mtx_Hess,Q)
        Gx =  np.dot(mtx_transp,pot)
        return Gx,pot
    
    def second_state_pot(self,vector_grad,mtx_Hess,gap):
        qo = self.b0 + self.a0 + self.d0
        q1 = self.bonds + self.angles + self.dihedrals
        Q = np.array(([q1[i] - qo[i] for i in range(len(qo))]), dtype = float)
        product1 = np.dot(Q,vector_grad)
        product2 = np.dot(Q,np.dot(mtx_Hess,Q))
   
        gap_dynamic = gap + product1 + 0.5 * product2
        
        return gap_dynamic
 
         
    def E_potencial(self,potencial):
        qo = self.b0 + self.a0 + self.d0
        q1 = self.bonds + self.angles + self.dihedrals
        Q = np.array(([q1[i] - qo[i] for i in range(len(qo))]), dtype = float)
        E = 0.5 * np.dot(Q,potencial)
        return E

    def output(self,file,counter,number_atoms,element,EA,ET):
        if counter % 20 == 0: # geometries saved 
            file.write('%g \n' % number_atoms)
	    file.write(' EA %g  ET:1 %g  ET:2  %g \n' % (EA,ET[0],ET[1]))
            for i in range(number_atoms):
                if i ==0:
                    file.write('%s  0.   0.   0. \n' % (element[i]))
                else:
                    x = (self.cart[i].x-self.cart[0].x)*my.a0
                    y = (self.cart[i].y-self.cart[0].y)*my.a0
                    z = (self.cart[i].z-self.cart[0].z)*my.a0
                    file.write('%s  %g  %g  %g \n' % (element[i],x,y,z))


    def rangos(self,Hessii,Etot,n_interv):
        """ 0.5*k*DQ^2 = E_tot"""
        my.hess_check(Hessii)
        d  = [dict() for n in range(len(Hessii))]
        DQ = [sqrt(2.0 * Etot/hii) for hii in Hessii]
        n1 = 0
        for bond in self.bonds:
            n2 = 0  
            low = self.b0[n1] - DQ[n1]; interv = 2.*DQ[n1] /float(n_interv)
            for i in range(n_interv):
                rango = (low + i*interv, low + (i+1)*interv)
                d[n1][rango] = 0 
            n1 += 1
         
        n1 = len(self.b0) 
        for ang in self.angles:
            low = self.a0[n1-len(self.b0)] - DQ[n1]; interv = 2.*DQ[n1] /float(n_interv)
            for i in range(n_interv):
                rango = (low + i*interv, low + (i+1)*interv)
                d[n1][rango] = 0
            n1 += 1

        n1 = len(self.b0) + len(self.a0)
        cte =len(self.b0) + len(self.a0)
        for dieh in self.dihedrals:
            low = self.d0[n1-cte] - DQ[n1]; interv = 2.*DQ[n1] /float(n_interv)
            for i in range(n_interv):
                rango = (low + i*interv, low + (i+1)*interv)
                d[n1][rango] = 0
            n1 += 1

        return d     

    def histogram(self,freq): 
        n1 = 0
        for bond in self.bonds:
            for key in freq[n1]: 
                if key[0] <= bond and bond < key[1]:
                    freq[n1][key] += 1
            n1 += 1

        n1 = len(self.b0) 
        for ang in self.angles:
            for key in freq[n1]: 
                if key[0] <= ang and ang < key[1]:
                    freq[n1][key] += 1
            n1 += 1

        n1 = len(self.b0) + len(self.a0)
        for dieh in self.dihedrals:
            for key in freq[n1]: 
                if key[0] <= dieh and dieh < key[1]:
                    freq[n1][key] += 1
            n1 += 1
        
        return freq

class Thermostat:
    def __init__(self, Q1 = 0., Q2 = 0., x1 = 0., vx1 = 0., x2 = 0., vx2 = 0.):
        self.Q1 = Q1
        self.Q2 = Q2
        self.x1 = x1
        self.vx1 = vx1
        self.x2 = x2
        self.vx2 = vx2
        self.scale = 1.
     
def Freq(ET,D,FC,factor,Low,High):
    if ET > Low and ET < High: 
        x = (ET-Low) / factor 	  
        D[int(x)] += 1
 
    return D
    
def calc_grad_ext(F_ext,atom1,atom2):
    v1 = my.attr_val(atom1); v2 = my.attr_val(atom2)
#   vector = [v1[i]-v2[i] for i in range(3)]
    vector = [v2[i]-v1[i] for i in range(3)]
    norm = np.linalg.norm(vector)
    U = [x/norm for x in vector]
    g_ext = [F_ext*x for x in U]

    return g_ext

def check_anchor(F,atom1,atom2,numat):
    L = []
    for x in range(numat):
        if atom1 == x:
            L.append(Point(F[0],F[1],F[2]))
        elif atom2 == x:
            L.append(Point(-F[0],-F[1],-F[2]))
        else:
            L.append(Point())

    return L
def flat_points(lists_points):
    L = []
    for val in lists_points:
        x,y,z = my.attr_val(val)
        L.append(x)
        L.append(y)
        L.append(z)

    return L


###########  READING DATA ########################

# READING DATA FROM THE INPUT FILE

input = open('input.dat')
files = [line for line in input.read().split()]

#reading of the Cartesian Cordinates, number of Atoms, type of Atoms,
#Masses and Energy
file_st1 = files[0]

cord, numat, typ, symb, mass, st1_energy, st1_grad_cart, st1_hess_cart = \
my.initial_data(file_st1,'state1')


# READING CONECTIVITY 

#my.mtx_conect('conex.dat')
conexion =list(open('internas.dat'))
nbond = int(conexion[0])
nang  = int(conexion[nbond+1])
ndih  = int(conexion[nbond+nang+2])
bond = [tuple([int(m) for m in w.split()]) for w in conexion[1:nbond+1]]
ang  = [tuple([int(m) for m in w.split()]) for w in conexion[nbond+2:nbond+nang+2]]
dih  = [tuple([int(m) for m in w.split()]) for w in conexion[nbond+nang+3:nbond+nang+ndih+3]]
redun = bond + ang + dih
ndim = len(redun)

# The final temperature which you want to heat the molecule, in Kelvin'
T = float(files[3])

#The time in femtoseconds in which you want to heat up the system'
time = float(files[4])

#The time in Femtoseconds in which you want to run the dynamics'
run_time = float(files[5])
input.close()

# External Forces
# atomic unit of force  a.u. 8.238722e-8
mod_Fext = float(files[6])/my.auN

# Anchor points
anchor1 = int(files[7]) - 1
anchor2 = int(files[8]) - 1


#########  INSTANTIATION OF ATOMS ###############################

for i in xrange(len(symb)):
    symb[i] = Atoms()

# Assignation of mass and cartesian Coordinates

for m in xrange(numat):
    j = 3*m
    symb[m].cart.x = cord[j]
    symb[m].cart.y = cord[j+1]
    symb[m].cart.z = cord[j+2]
    symb[m].mass = mass[m]

######### INSTANTIATION OF MOLECULES #########################

# Instantiation of the molecule
mol = Molecule([var.cart for var in symb] )

# Transformation to internal coordinates
mol.internas_bonds(bond,nbond)
mol.internas_ang(ang,nang)			
mol.internas_dihe(dih,ndih)

# Internal coordinates of the minimun
mol.b0 = mol.bonds
mol.a0 = mol.angles
mol.d0 = mol.dihedrals
q0 = mol.bonds + mol.angles + mol.dihedrals


## CALCULATING THE GRADIENT AND HESSIAN MATRIX IN INTERNAL COORDINATES 
## FOR THE FIRST STATE

second_derv = my.segunda_wilson(symb,ndim,numat,bond,ang,dih)
derv_trans = np.transpose(second_derv)
Bwilson,transp = my.matrix_transf(symb,bond,ang,dih)
G_mtx = np.dot(Bwilson,transp)
G_inv = my.invertir_mtx(G_mtx)

Grad_st1 =np.dot(G_inv,np.dot(Bwilson,st1_grad_cart))
mtx_B_G1 = np.dot(derv_trans,Grad_st1)
mtx_resta1 = st1_hess_cart - mtx_B_G1
Hess_st1 = np.dot(np.dot(np.dot(G_inv,Bwilson),mtx_resta1),np.dot(transp,G_inv))


########### DATA OF THE SECOND AND THIRD STATES ############
#############################################################
file_st2 = files[1]
st2_energy, st2_grad_cart, st2_hess_cart = my.initial_data(file_st2,'state2')

# CALCULATING THE GRADIENT AND HESSIAN MATRIX IN INTERNAL COORDINATES 
# FOR THE SECOND STATE
Grad_st2 =np.dot(G_inv,np.dot(Bwilson,st2_grad_cart))

mtx_B_G2 = np.dot(derv_trans,Grad_st2)
mtx_resta2 = st2_hess_cart - mtx_B_G2
Hess_st2 =np.dot(np.dot(np.dot(G_inv,Bwilson),mtx_resta2),np.dot(transp,G_inv))

########### DATA OF THE SECOND AND THIRD STATES ############
###################################################################
file_st3 = files[2]
st3_energy, st3_grad_cart, st3_hess_cart = my.initial_data(file_st3,'state2')

# CALCULATING THE GRADIENT AND HESSIAN MATRIX IN INTERNAL COORDINATES 
# FOR THE SECOND STATE
Grad_st3 =np.dot(G_inv,np.dot(Bwilson,st3_grad_cart))

mtx_B_G3 = np.dot(derv_trans,Grad_st3)
mtx_resta3 = st3_hess_cart - mtx_B_G3
Hess_st3 =np.dot(np.dot(np.dot(G_inv,Bwilson),mtx_resta3),np.dot(transp,G_inv))


############# HEATING USING MAXWELL-BOLTZMANN DISTRIBUTION ###################
# " From Statistical thermodynamics it is known that the velocities,         #
#  of the atoms in a classical system are distributed according              #
# to the Maxwell-Boltzmann distribution.This says that if the temperature    #
# of the system is T, the probability of each component of the velocity      #
# of the ith atom  having a value between v and v + dv is                    #
#                                                                            # 
# f(V) dV = sqrt(Massi /2 pi Kb T) * exp(-Massi V^2 / 2 Kb T) * dV           #
#                                                                            #  
#The values of the velocities of the atoms can be assigned by treating them  #
# as independent Gaussian random variables drawm from the above distribution"#
# with mean value of 0 and standart deviation of sqrt(KbT/Massi)             #
##############################################################################

## OUTPUT FILES
total_energy = open('energy.dat', 'w')
kinetic = open('kinetic_energy.out','w')
geometry = open('geometry.dat','w')
temperature = open( 'temperature.out','w')
Eact = open('Activation_energy.out','w')
probability = open('Probability_Maxwell','w')
Normal = open('Normal_dist.out','w')

ET_gap1 = open('ET_GAP1.out','w')
ET_gap2 = open('ET_GAP2.out','w')

########### HEATING #################
d0,ndih_inicial,dih_inicial,Hess_inicial = mol.d0,ndih,dih,Hess_st1

# time of integration and heating time per unit of temperature
dt = 0.1/my.au_time; slice = (time/my.au_time)/T

# Initial heating temperature
Temper = T
ncounter2 = 1
t_heat = 0.

######### Properties of the Thermostat#########################
# The relaxation factor range between 0.5 and 2 Ps generally  #
thermo = Thermostat()                                         #
freq = 1./(2.2e1/my.au_time)                                                                 
thermo.Q1 = 3*numat*Temper*my.kb/(freq**2)                    #  
thermo.Q2 = Temper*(freq**2)/my.kb                            #
###############################################################


# Initial random velocities

for n in range(numat):
    symb[n].vel.x = my.random_maxwell(symb[n],2*T)
    symb[n].vel.y = my.random_maxwell(symb[n],2*T)
    symb[n].vel.z = my.random_maxwell(symb[n],2*T)


while t_heat <= time/my.au_time:

    Mass_total = sum(symb[k].mass for k in range(numat))
    MC_vel_X = sum(symb[k].mass*symb[k].vel.x for k in xrange(numat))/Mass_total
    MC_vel_Y = sum(symb[k].mass*symb[k].vel.y for k in xrange(numat))/Mass_total
    MC_vel_Z = sum(symb[k].mass*symb[k].vel.z for k in xrange(numat))/Mass_total

    angular_vel = my.Inertia_matrix(symb)
    rot_comp = [my.cross(angular_vel,[symb[k].cart.x,symb[k].cart.y,symb[k].cart.z]) for k in range(numat)]
    
    for i in xrange(numat):
        symb[i].vel.x = symb[i].vel.x - MC_vel_X - rot_comp[i][0]
        symb[i].vel.y = symb[i].vel.y - MC_vel_Y - rot_comp[i][1]
        symb[i].vel.z = symb[i].vel.z - MC_vel_Z - rot_comp[i][2]

    my.nose_hoover_1(symb,dt,thermo,T,numat)

# Transformation to internal coordinates
    mol.cart =[var.cart for var in symb]
    mol.internas_bonds(bond,nbond)
    mol.internas_ang(ang,nang)
    mol.internas_dihe(dih,ndih)
    
# Calculation of the transpose Wilson Matrix
    wilson, transp = my.matrix_transf(symb,bond,ang,dih)

# force field calculation
    grad,pot = mol.potential(transp,Grad_st1,Hess_st1)
    vect_ext = calc_grad_ext(mod_Fext, symb[anchor1].cart, symb[anchor2].cart)
    grad_ext = flat_points(check_anchor(vect_ext,anchor1,anchor2,numat))
    new_grad = grad  + grad_ext
    for i in range(numat):
        symb[i].gradiente((new_grad[3*i],new_grad[1+3*i],new_grad[2+3*i]))


# last step
    my.nose_hoover_2(symb,dt,thermo,T,numat)

    Ek =  sum([symb[i].cinetica() for i in range(numat)])
    Ep = mol.E_potencial(pot)
    Etot = Ek + Ep
    t_heat += dt

#Output 
#   total_energy.write('%g \n' % Etot)
#   kinetic.write('%g \n' % Ek)
#   mol.output(geometry,ncounter2,numat,typ, 0.0, 0.0)
        
###### RUNNING OF THE DYNAMICS #######################
freq = 1./(2.2e1/my.au_time)                                                                 
thermo.Q1 = 3*numat*Temper*my.kb/(freq**2)                      
thermo.Q2 = Temper*my.kb/(freq**2)
ndym = 0 
#####################################################

# Sigma = Its the maximum allow variation of the energetic gap

## Averages
average_T = list(); vel_prom = list()
ET1_mean = list(); ET2_mean = list()

FC_ET1 = (st2_energy - st1_energy )*my.eh_kcal
FC_ET2 = ((st3_energy - st1_energy))*my.eh_kcal 
sigma = 40. ; boxes = 200
Low1 = FC_ET1 - sigma; High1 = FC_ET1 + sigma
Low2 = FC_ET2 - sigma; High2 = FC_ET2 + sigma
histo1, histo2 = dict(), dict()
for k in range(boxes):
    histo1[k], histo2[k] = 0, 0

#############DYNAMICS ######################

while run_time/my.au_time >= dt * float(ndym):

    Mass_total = sum(symb[k].mass for k in range(numat))
    MC_vel_X = sum(symb[k].mass*symb[k].vel.x for k in xrange(numat))/Mass_total
    MC_vel_Y = sum(symb[k].mass*symb[k].vel.y for k in xrange(numat))/Mass_total
    MC_vel_Z = sum(symb[k].mass*symb[k].vel.z for k in xrange(numat))/Mass_total

    angular_vel = my.Inertia_matrix(symb)
    rot_comp = [my.cross(angular_vel,[symb[k].cart.x,symb[k].cart.y,symb[k].cart.z]) for k in range(numat)]
 

    for i in xrange(numat):
        symb[i].vel.x = symb[i].vel.x - MC_vel_X - rot_comp[i][0]
        symb[i].vel.y = symb[i].vel.y - MC_vel_Y - rot_comp[i][1]
        symb[i].vel.z = symb[i].vel.z - MC_vel_Z - rot_comp[i][2]

    my.nose_hoover_1(symb,dt,thermo,T,numat)

# Transformation to internal coordinates
    mol.cart =[var.cart for var in symb]
    mol.internas_bonds(bond,nbond)
    mol.internas_ang(ang,nang)
    mol.internas_dihe(dih,ndih)

# Calculation of the transpose Wilson Matrix
    wilson, transp = my.matrix_transf(symb,bond,ang,dih)

# force field calculation
    grad,pot = mol.potential(transp,Grad_st1,Hess_st1)
    vect_ext = calc_grad_ext(mod_Fext, symb[anchor1].cart, symb[anchor2].cart)
    grad_ext = flat_points(check_anchor(vect_ext,anchor1,anchor2,numat))
    new_grad = grad  + grad_ext
    for i in range(numat):
        symb[i].gradiente((new_grad[3*i],new_grad[1+3*i],new_grad[2+3*i]))

# last step of the dynamic
    my.nose_hoover_2(symb,dt,thermo,T,numat)

    ndym += 1
    
# COMPUTING THE ENERGY AND THE GAP BETWEEN S0 AND DIFFERENT ENERGETIC LEVELS
    Ek = sum(symb[i].cinetica() for i in range(numat))
    Ep = mol.E_potencial(pot)
    Etot = Ek + Ep
    prom_T = 2.0 * Ek /(3*numat*my.kb) 
    average_T.append(prom_T)
    temperature.write(' %g \n' % prom_T)
   
####Maxwell-Boltzmann distribution
    prob_maxwell = my.maxwell_vel(symb[6],T)
    vel_prom.append(symb[6].vel.x)
    probability.write('%g %g \n' % (prob_maxwell[1],prob_maxwell[0]))
  
	
#Output 
    total_energy.write('%g \n' % Etot)
    kinetic.write('%g \n' % Ek)

    vertical1 =   st2_energy - st1_energy 
    vertical2 =   (st3_energy - st1_energy)
    second_Epot = mol.second_state_pot(Grad_st2,Hess_st2,vertical1)
    third_Epot  = mol.second_state_pot(Grad_st3,Hess_st3,vertical2)
    gap1 = (second_Epot - Ep)* my.eh_kcal
    gap2 = (third_Epot  - Ep)* my.eh_kcal

    ET1_mean.append(gap1); ET2_mean.append(gap2)
    histo1 = Freq(gap1,histo1,FC_ET1,(2*sigma/boxes),Low1,High1)
    histo2 = Freq(gap2,histo2,FC_ET2,(2*sigma/boxes),Low2,High2)

    act_energy = Ep*my.eh_kcal
    Eact.write('%g \n' % act_energy)
    ET_gap1.write( '%g \n' % gap1)
    ET_gap2.write( '%g \n' % gap2)
    mol.output(geometry,ndym,numat,typ,act_energy,(gap1,gap2))

Eact.close(); total_energy.close(); kinetic.close()


########Writting of the averages######################
T_array = np.array(average_T)
std = np.std(T_array); mean = np.mean(T_array)
print std, mean, 'Stadistics Temperature'


f1 = open('distribution_transition_S0_S1.out','w')
f2 = open('distribution_transition_S0_S2.out','w')
for k in range(boxes) :
    x1 = ((Low1 + 2*k*sigma/boxes) + (Low1 + 2*(k+1)*sigma/boxes))*0.5
    x2 = ((Low2 + 2*k*sigma/boxes) + (Low2 + 2*(k+1)*sigma/boxes))*0.5
    f1.write('{0:>8f} {1:4d} \n'.format(x1,histo1[k]))
    f2.write('{0:>8f} {1:4d} \n'.format(x2,histo2[k]))
f1.close(); f2.close()

