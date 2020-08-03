"""Run molecular dynamics."""

from typing import NamedTuple, Optional

import numpy as np
from scipy.constants import physical_constants

from .molecule import Molecule

# from femtoseconds to au
AU_TIME = 1e15 * physical_constants["atomic unit of time"][0]


class Configuration(NamedTuple):
    """NamedTuple that contains the configuration values to run a simulation."""

    molecule: Molecule
    connectivity: np.ndarray
    hessian: np.ndarray
    gradient: Optional[np.ndarray]
    time: int
    dt: float


def run_simulation(config: Configuration) -> None:
    """Run a MD simulation using a given `config`."""
    mol = config.molecule
    # Initialize the velocities
    mol.generate_random_velocities()

    # convert time delta to AU
    dt = config.dt * AU_TIME


# def dif(a, b):
#     q = (a.x - b.x)**2 + (a.y - b.y)**2 + (a.z - b.z)**2
#     return q


# def dif2(a, b, c):
#     q = (a.x-b.x)*(c.x-b.x)+(a.y-b.y)*(c.y-b.y) + (a.z-b.z)*(c.z-b.z)
#     return q


# def atom_type(molecule):
#     w = list()
#     element = {1: 'H', 6.: 'C', 7.: 'N', 8.: 'O'}
#     for val in molecule:
#         if val in element:
#             w = w + list(element[val])
#     return w


# def cuad(a, b):
#     q = a.x*b.x + a.y*b.y + a.z*b.z
#     return q


# def rad(q):
#     w = q*180./pi
#     return w


# def attr_val(A):
#     lista = [A.x, A.y, A.z]

#     return lista


# def kronecker(m, n):
#     if m == n:
#         kron = 1.
#     else:
#         kron = 0.

#     return kron


# def cross(A, B):
#     prod = list(range(3))
#     try:
#         prod[0] = A[1]*B[2]-A[2]*B[1]
#         prod[1] = A[2]*B[0]-A[0]*B[2]
#         prod[2] = A[0]*B[1]-A[1]*B[0]

#     except TypeError:
#         prod[0] = A.y * B.z - A.z * B.y
#         prod[1] = A.z * B.x - A.x * B.z
#         prod[2] = A.x * B.y - A.y * B.x

#     return prod


# def invertir_mtx(A):
#     ndim = len(A)
#     eig_val, eig_vec = np.linalg.eig(A)
#     unitaria = np.eye(ndim, dtype=float)
#     eig_trans = np.transpose(eig_vec)
#     for i in range(ndim):
#         if abs(eig_val[i]) > 1.e-5:
#             unitaria[i, i] = 1./eig_val[i]
#         else:
#             unitaria[i, i] = 0.

#     mtx_inversa = np.dot(eig_vec, np.dot(unitaria, eig_trans))

#     return mtx_inversa


# def scalar_vector(s, V):
#     R = Point()
#     R.x = s*V.x
#     R.y = s*V.y
#     R.z = s*V.z

#     return R


###################### VELOCITY VERLET ################################################

# def move_first(mol: Molecule, dt: float):
#     """AT THE START OF A TIMESTEP, MOVEA IS CALLED TO ADVANCE THE
#     POSITIONS AND 'HALF-ADVANCE' THE VELOCITIES.  THEN THE FORCE
#     ROUTINE IS CALLED, AND THIS IS FOLLOWED BY MOVEB WHICH
#     COMPLETES THE ADVANCEMENT OF VELOCITIES. WHERE mol DEFINE A LIST
#     THAT CONTAINS ALL THE ATOMS IN THE MOLECULE"""
#     dt2 = dt*.5
#     dtsq2 = dt * dt2

#     for i in xrange(len(mol)):
#         mol[i].cart.x = mol[i].cart.x + dt * mol[i].vel.x + \
#             dtsq2 * mol[i].force.x / mol[i].mass
#         mol[i].cart.y = mol[i].cart.y + dt * mol[i].vel.y + \
#             dtsq2 * mol[i].force.y / mol[i].mass
#         mol[i].cart.z = mol[i].cart.z + dt * mol[i].vel.z + \
#             dtsq2 * mol[i].force.z / mol[i].mass
#         mol[i].vel.x = mol[i].vel.x + dt2 * mol[i].force.x / mol[i].mass
#         mol[i].vel.y = mol[i].vel.y + dt2 * mol[i].force.y / mol[i].mass
#         mol[i].vel.z = mol[i].vel.z + dt2 * mol[i].force.z / mol[i].mass


# def MOVEB(mol, dt):
#     """ Second part of the velocity verlet"""
#     dt2 = dt*.5
#     for i in xrange(len(mol)):
#         mol[i].vel.x = mol[i].vel.x + dt2 * mol[i].force.x / mol[i].mass
#         mol[i].vel.y = mol[i].vel.y + dt2 * mol[i].force.y / mol[i].mass
#         mol[i].vel.z = mol[i].vel.z + dt2 * mol[i].force.z / mol[i].mass

# ################# INTEGRATION OF NOSE-HOOVER THERMOSTAT ###################################
# # Integration of a molecular system coupled to a bath, through the Liouville Operator
# # Using the trotter factorization


# def bath(mol, dt, thermo, Ek, T, numat):
#     dt2 = dt*.5
#     dt4 = dt*0.25
#     dt8 = dt*0.125

#     G2 = (thermo.Q1*(thermo.vx1**2) - T*kb)/thermo.Q2
#     thermo.vx2 += G2*dt4
#     thermo.vx1 *= exp(-thermo.vx2*dt8)
#     G1 = (2*Ek - 3*numat*T*kb)/thermo.Q1
#     thermo.vx1 += G1*dt4
#     thermo.vx1 *= exp(-thermo.vx2*dt8)
#     thermo.scale = exp(-thermo.vx1*dt2)

#     for i in range(numat):
#         mol[i].vel.x *= thermo.scale
#         mol[i].vel.y *= thermo.scale
#         mol[i].vel.z *= thermo.scale

#     Ek *= thermo.scale**2.

#     thermo.x1 += thermo.vx1*dt2
#     thermo.x2 += thermo.vx2*dt2
#     thermo.vx1 *= exp(-thermo.vx2*dt8)
#     G1 = (2*Ek - 3*numat*T*kb)/thermo.Q1
#     thermo.vx1 += G1*dt4
#     thermo.vx1 *= exp(-thermo.vx2*dt8)
#     G2 = (thermo.Q1*(thermo.vx1**2) - T*kb)/thermo.Q2
#     thermo.vx2 += G2*dt4

#     return mol


# def nose_hoover_1(mol, dt, thermo, T, numat):

#     Ekinetic = sum([0.5*mol[j].mass*cuad(mol[j].vel, mol[j].vel)
#                     for j in range(numat)])

#     mol = bath(mol, dt, thermo, Ekinetic, T, numat)

#     MOVEA(mol, dt)


# def nose_hoover_2(mol, dt, thermo, T, numat):

#     MOVEB(mol, dt)

#     Ekinetic = sum([0.5*mol[j].mass*cuad(mol[j].vel, mol[j].vel)
#                     for j in range(numat)])

#     mol = bath(mol, dt, thermo, Ekinetic, T, numat)


# # ROTATIONAL MOTION

# def Inertia_matrix(atom):
#     numat = len(atom)
#     I_matrix = np.zeros((3, 3), dtype=float)
#     L = np.zeros((3), dtype=float)

#     I_matrix[0, 0] = sum(
#         [atom[n].mass * (atom[n].cart.y**2 + atom[n].cart.z**2) for n in range(numat)])
#     I_matrix[0, 1] = -sum([atom[n].mass*atom[n].cart.x *
#                            atom[n].cart.y for n in range(numat)])
#     I_matrix[0, 2] = -sum([atom[n].mass*atom[n].cart.x *
#                            atom[n].cart.z for n in range(numat)])

#     I_matrix[1, 0] = I_matrix[0, 1]
#     I_matrix[1, 1] = sum(
#         [atom[n].mass * (atom[n].cart.x**2 + atom[n].cart.z**2) for n in range(numat)])
#     I_matrix[1, 2] = -sum([atom[n].mass*atom[n].cart.y *
#                            atom[n].cart.z for n in range(numat)])

#     I_matrix[2, 0] = I_matrix[0, 2]
#     I_matrix[2, 1] = I_matrix[1, 2]
#     I_matrix[2, 2] = sum(
#         [atom[n].mass * (atom[n].cart.x**2 + atom[n].cart.y**2) for n in range(numat)])

#     # Angular Momentum
#     lin_mom = [scalar_vector(atom[k].mass, atom[k].vel) for k in range(numat)]

#     lista = [cross(atom[n].cart, lin_mom[n]) for n in range(numat)]
#     L[0] = sum([lista[n][0] for n in range(numat)])
#     L[1] = sum([lista[n][1] for n in range(numat)])
#     L[2] = sum([lista[n][2] for n in range(numat)])

#     I_inverse = np.linalg.inv(I_matrix)
#     ang_vel = np.dot(I_inverse, L)

#     return ang_vel

# ################## Coordinates transformation ############


# def enlaces(A_1, A_2):
#     d = sqrt(dif(A_1, A_2))
#     return d


# def angulos(A_1, A_2, A_3):
#     w1 = dif2(A_1, A_2, A_3)
#     r1 = sqrt(dif(A_1, A_2))
#     r2 = sqrt(dif(A_3, A_2))
#     return acos(w1/(r1*r2))


# def dihedros(A_1, A_2, A_3, A_4, qt_1):
#     xba = A_2.x - A_1.x
#     xca = A_3.x - A_1.x
#     yba = A_2.y - A_1.y
#     yca = A_3.y - A_1.y
#     zba = A_2.z - A_1.z
#     zca = A_3.z - A_1.z
#     xcb = A_3.x - A_2.x
#     xdb = A_4.x - A_2.x
#     ycb = A_3.y - A_2.y
#     ydb = A_4.y - A_2.y
#     zcb = A_3.z - A_2.z
#     zdb = A_4.z - A_2.z
#     nx1 = (yba*zca)-(yca*zba)
#     ny1 = (xca*zba)-(xba*zca)
#     nz1 = (xba*yca)-(xca*yba)
#     nx2 = (ycb*zdb)-(ydb*zcb)
#     ny2 = (xdb*zcb)-(xcb*zdb)
#     nz2 = (xcb*ydb)-(xdb*ycb)
#     pdot = (nx1*nx2)+(ny1*ny2)+(nz1*nz2)
#     w1 = sqrt((nx1*nx1)+(ny1*ny1)+(nz1*nz1))
#     w2 = sqrt((nx2*nx2)+(ny2*ny2)+(nz2*nz2))
#     signo = (nx2*xba)+(ny2*yba)+(nz2*zba)
#     u = pdot/(w1*w2)
#     if abs(u) >= 1.0:
#         u = (1.-1.e-20)*(u/abs(u))
#     if signo >= 0.:
#         dieh = acos(u)
#         qa = dieh - 2*pi
#     if signo < 0.:
#         dieh = -acos(u)
#         qa = dieh + 2*pi

#     wa = dieh - qt_1
#     wb = qa - qt_1
#     if abs(wa) <= abs(wb):
#         ang_dih = dieh
#     else:
#         ang_dih = qa

#     return ang_dih


# ############## Wilson's Matrix Calculation ################

# # Distance derivative

# def derv_bond(A, B):
#     r = sqrt(dif(A, B))
#     q = [(A.x-B.x)/r, (A.y-B.y)/r, (A.z-B.z)/r]
#     return q


# def derv_ang(A, B, C):
#     t1 = list()
#     t2 = list()
#     t3 = list()
#     r1 = sqrt(dif(A, B))
#     r2 = sqrt(dif(B, C))
#     teta = angulos(A, B, C)
#     d1 = [B.x - C.x, B.y - C.y, B.z - C.z]
#     d2 = [A.x - B.x, A.y - B.y, A.z - B.z]
#     d3 = [2.0*B.x-A.x-C.x, 2.0*B.y-A.y-C.y, 2.0*B.z-A.z-C.z]

#     for i in range(3):
#         a1 = d1[i]/(r1*r2*sin(teta))
#         a2 = (1.0/tan(teta))*(d2[i]/(r1*r1))
#         t1.append(a1+a2)
#         b1 = (1.0/tan(teta))*(((-1*d2[i])/(r1*r1))+(d1[i]/(r2*r2)))
#         b2 = d3[i]/(r1*r2*sin(teta))
#         t2.append(b1-b2)
#         c1 = (-1*d2[i])/(r1*r2*sin(teta))
#         c2 = (1.0/tan(teta))*((-1*d1[i])/(r2*r2))
#         t3.append(c1+c2)
#     return (t1, t2, t3)


# def derv_dih(A, B, C, D):
#     r2 = sqrt(dif(B, C))
#     uxw = np.zeros((3), dtype=float)
#     vxw = np.zeros((3), dtype=float)
#     v1 = [A.x - B.x, A.y - B.y, A.z - B.z]
#     v2 = [C.x - B.x, C.y - B.y, C.z - B.z]
#     v3 = [C.x - D.x, C.y - D.y, C.z - D.z]

#     uxw = cross(v1, v2)
#     vxw = cross(v2, v3)

#     omega1 = sum([uxw[i]*uxw[i] for i in range(3)])
#     omega2 = sum([vxw[i]*vxw[i] for i in range(3)])

#     DA = [(r2/omega1)*uxw[i] for i in range(3)]
#     DD = [-1*(r2/omega2)*vxw[i] for i in range(3)]

#     v1pv2 = sum([v1[i]*v2[i] for i in range(3)])
#     v3pv2 = sum([v3[i]*v2[i] for i in range(3)])
#     p1 = v1pv2/(r2*r2) - 1.0
#     p2 = v3pv2/(r2*r2)
#     DB = [(p1*DA[i])-(p2*DD[i]) for i in range(3)]

#     p3 = v3pv2/(r2*r2) - 1.0
#     p4 = v1pv2/(r2*r2)
#     DC = [(p3 * DD[i])-(p4*DA[i]) for i in range(3)]

#     return (DA, DB, DC, DD)


# def matrix_transf(mol, bond, ang, dih):
#     """ Function for calculating the transpose Wilson Matrix"""
#     numat = len(mol)
#     ndim = len(bond) + len(ang) + len(dih)
#     wilson = np.zeros((ndim, 3*numat), dtype=float)

#     for i in range(len(bond)):
#         m1 = bond[i][0]
#         n1 = 3*(m1-1)
#         m2 = bond[i][1]
#         n2 = 3*(m2-1)
#         wilson[i, n1:n1+3] = derv_bond(mol[m1-1].cart, mol[m2-1].cart)
#         wilson[i, n2:n2+3] = derv_bond(mol[m2-1].cart, mol[m1-1].cart)

#     for i in range(len(ang)):
#         m1 = ang[i][0]
#         m2 = ang[i][1]
#         m3 = ang[i][2]
#         j = i + len(bond)
#         n1 = 3*(m1-1)
#         n2 = 3*(m2-1)
#         n3 = 3*(m3-1)
#         wilson[j, n1:n1+3], wilson[j, n2:n2+3], wilson[j, n3:n3+3] = \
#             derv_ang(mol[m1-1].cart, mol[m2-1].cart, mol[m3-1].cart)

#     for i in range(len(dih)):
#         m1 = dih[i][0]
#         m2 = dih[i][1]
#         m3 = dih[i][2]
#         m4 = dih[i][3]
#         j = i+len(bond)+len(ang)
#         n1 = 3*(m1-1)
#         n2 = 3*(m2-1)
#         n3 = 3*(m3-1)
#         n4 = 3*(m4-1)
#         wilson[j, n1:n1+3], wilson[j, n2:n2+3], wilson[j, n3:n3+3], wilson[j, n4:n4+3] = \
#             derv_dih(mol[m1-1].cart, mol[m2-1].cart,
#                      mol[m3-1].cart, mol[m4-1].cart)

#     Btrans = wilson.transpose()
#     return wilson, Btrans

# # Second derivatives


# def seg_derv_bond(A, B, m, n):
#     r = sqrt(dif(A, B))
#     d_ij = kronecker(m, n)
#     Ax = attr_val(A)
#     Bx = attr_val(B)
#     derv = (-1./r)*(((Ax[m]-Bx[m])*(Ax[n]-Bx[n])/(r*r)) - d_ij)

#     return derv


# def seg_derv_ang(A, B, C, DA, DB, DC, i, j):
#     Ax = attr_val(A)
#     Bx = attr_val(B)
#     Cx = attr_val(C)
#     r1 = sqrt(dif(A, B))
#     r2 = sqrt(dif(C, B))
#     d_ij = kronecker(i, j)
#     teta = angulos(A, B, C)
#     a1 = cos(teta)
#     a2 = sin(teta)
#     a3 = tan(teta)

#     p1 = (Ax[i]-Bx[i])*(Cx[j]-Bx[j])/(r1*r2)
#     p2 = (Ax[i]-Bx[i])*(Ax[j]-Bx[j])/(r1*r1)
#     p3 = (Ax[j]-Bx[j])*(Cx[i]-Bx[i])/(r1*r2)
#     p4 = (Cx[j]-Bx[j])*(Cx[i]-Bx[i])/(r2*r2)

#     t1 = (p1 + p3 - 3*p2*a1 + d_ij*a1)/(a2*r1*r1)
#     t2 = (p3 + p1 - 3*p4*a1 + d_ij*a1)/(a2*r2*r2)
#     t3 = (p2 + p4 - p1*a1 - d_ij)/(a2*r1*r2)
#     t4 = (p4 + p2 - p3*a1 - d_ij)/(a2*r1*r2)

#     dos_AA = t1 - DA[i]*DA[j]/a3
#     dos_AB = -t1 - t3 - DA[i]*DB[j]/a3
#     dos_AC = t3 - DA[i]*DC[j]/a3

#     dos_BA = -t1 - t4 - DB[i]*DA[j]/a3
#     dos_BB = t1 + t2 + t3 + t4 - DB[i]*DB[j]/a3
#     dos_BC = -t2 - t3 - DB[i]*DC[j]/a3

#     dos_CA = t4 - DC[i]*DA[j]/a3
#     dos_CB = -t2 - t4 - DC[i]*DB[j]/a3
#     dos_CC = t2 - DC[i]*DC[j]/a3

#     dervs_A = dos_AA, dos_AB, dos_AC
#     dervs_B = dos_BA, dos_BB, dos_BC
#     dervs_C = dos_CA, dos_CB, dos_CC

#     return dervs_A, dervs_B, dervs_C


# def numerical_derv_AA(i, j, A, B, C, dx):
#     RCB = sqrt(dif(C, B))
#     rcb = [C.x - B.x, C.y - B.y, C.z - B.z]

#     ra = [A.x, A.y, A.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [ra[n] + 2*DX[n] for n in range(3)]
#     D2 = [ra[n] + DX[n] for n in range(3)]
#     D3 = [ra[n] - DX[n] for n in range(3)]
#     D4 = [ra[n] - 2*DX[n] for n in range(3)]

#     A1 = Point(D1[0], D1[1], D1[2])
#     A2 = Point(D2[0], D2[1], D2[2])
#     A3 = Point(D3[0], D3[1], D3[2])
#     A4 = Point(D4[0], D4[1], D4[2])

#     rab_D1 = [A1.x - B.x, A1.y - B.y, A1.z - B.z]
#     rab_D2 = [A2.x - B.x, A2.y - B.y, A2.z - B.z]
#     rab_D3 = [A3.x - B.x, A3.y - B.y, A3.z - B.z]
#     rab_D4 = [A4.x - B.x, A4.y - B.y, A4.z - B.z]

#     rmb_D1 = cross(rab_D1, rcb)
#     RMB_D1 = sqrt(sum([r*r for r in rmb_D1]))
#     rmb_D2 = cross(rab_D2, rcb)
#     RMB_D2 = sqrt(sum([r*r for r in rmb_D2]))
#     rmb_D3 = cross(rab_D3, rcb)
#     RMB_D3 = sqrt(sum([r*r for r in rmb_D3]))
#     rmb_D4 = cross(rab_D4, rcb)
#     RMB_D4 = sqrt(sum([r*r for r in rmb_D4]))

#     f1 = RCB*rmb_D1[i]/(RMB_D1**2)
#     f2 = RCB*rmb_D2[i]/(RMB_D2**2)
#     f3 = RCB*rmb_D3[i]/(RMB_D3**2)
#     f4 = RCB*rmb_D4[i]/(RMB_D4**2)

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_AB(i, j, A, B, C, dx):

#     rb = [B.x, B.y, B.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [rb[n] + 2*DX[n] for n in range(3)]
#     D2 = [rb[n] + DX[n] for n in range(3)]
#     D3 = [rb[n] - DX[n] for n in range(3)]
#     D4 = [rb[n] - 2*DX[n] for n in range(3)]

#     B1 = Point(D1[0], D1[1], D1[2])
#     B2 = Point(D2[0], D2[1], D2[2])
#     B3 = Point(D3[0], D3[1], D3[2])
#     B4 = Point(D4[0], D4[1], D4[2])

#     rab_D1 = [A.x - B1.x, A.y - B1.y, A.z - B1.z]
#     rab_D2 = [A.x - B2.x, A.y - B2.y, A.z - B2.z]
#     rab_D3 = [A.x - B3.x, A.y - B3.y, A.z - B3.z]
#     rab_D4 = [A.x - B4.x, A.y - B4.y, A.z - B4.z]

#     rcb_D1 = [C.x - B1.x, C.y - B1.y, C.z - B1.z]
#     RCB_D1 = sqrt(sum([r*r for r in rcb_D1]))
#     rcb_D2 = [C.x - B2.x, C.y - B2.y, C.z - B2.z]
#     RCB_D2 = sqrt(sum([r*r for r in rcb_D2]))
#     rcb_D3 = [C.x - B3.x, C.y - B3.y, C.z - B3.z]
#     RCB_D3 = sqrt(sum([r*r for r in rcb_D3]))
#     rcb_D4 = [C.x - B4.x, C.y - B4.y, C.z - B4.z]
#     RCB_D4 = sqrt(sum([r*r for r in rcb_D4]))

#     rmb_D1 = cross(rab_D1, rcb_D1)
#     RMB_D1 = sqrt(sum([r*r for r in rmb_D1]))
#     rmb_D2 = cross(rab_D2, rcb_D2)
#     RMB_D2 = sqrt(sum([r*r for r in rmb_D2]))
#     rmb_D3 = cross(rab_D3, rcb_D3)
#     RMB_D3 = sqrt(sum([r*r for r in rmb_D3]))
#     rmb_D4 = cross(rab_D4, rcb_D4)
#     RMB_D4 = sqrt(sum([r*r for r in rmb_D4]))

#     f1 = RCB_D1*rmb_D1[i]/(RMB_D1**2)
#     f2 = RCB_D2*rmb_D2[i]/(RMB_D2**2)
#     f3 = RCB_D3*rmb_D3[i]/(RMB_D3**2)
#     f4 = RCB_D4*rmb_D4[i]/(RMB_D4**2)

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_AC(i, j, A, B, C, dx):
#     rab = [A.x - B.x, A.y - B.y, A.z - B.z]

#     rc = [C.x, C.y, C.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [rc[n] + 2*DX[n] for n in range(3)]
#     D2 = [rc[n] + DX[n] for n in range(3)]
#     D3 = [rc[n] - DX[n] for n in range(3)]
#     D4 = [rc[n] - 2*DX[n] for n in range(3)]

#     C1 = Point(D1[0], D1[1], D1[2])
#     C2 = Point(D2[0], D2[1], D2[2])
#     C3 = Point(D3[0], D3[1], D3[2])
#     C4 = Point(D4[0], D4[1], D4[2])

#     rcb_D1 = [C1.x - B.x, C1.y - B.y, C1.z - B.z]
#     RCB_D1 = sqrt(sum([r*r for r in rcb_D1]))
#     rcb_D2 = [C2.x - B.x, C2.y - B.y, C2.z - B.z]
#     RCB_D2 = sqrt(sum([r*r for r in rcb_D2]))
#     rcb_D3 = [C3.x - B.x, C3.y - B.y, C3.z - B.z]
#     RCB_D3 = sqrt(sum([r*r for r in rcb_D3]))
#     rcb_D4 = [C4.x - B.x, C4.y - B.y, C4.z - B.z]
#     RCB_D4 = sqrt(sum([r*r for r in rcb_D4]))

#     rmb_D1 = cross(rab, rcb_D1)
#     RMB_D1 = sqrt(sum([r*r for r in rmb_D1]))
#     rmb_D2 = cross(rab, rcb_D2)
#     RMB_D2 = sqrt(sum([r*r for r in rmb_D2]))
#     rmb_D3 = cross(rab, rcb_D3)
#     RMB_D3 = sqrt(sum([r*r for r in rmb_D3]))
#     rmb_D4 = cross(rab, rcb_D4)
#     RMB_D4 = sqrt(sum([r*r for r in rmb_D4]))

#     f1 = RCB_D1*rmb_D1[i]/(RMB_D1**2)
#     f2 = RCB_D2*rmb_D2[i]/(RMB_D2**2)
#     f3 = RCB_D3*rmb_D3[i]/(RMB_D3**2)
#     f4 = RCB_D4*rmb_D4[i]/(RMB_D4**2)

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_BA(i, j, A, B, C, D, dx):

#     RCB = sqrt(dif(C, B))
#     rcb = [C.x - B.x, C.y - B.y, C.z - B.z]
#     rcd = [C.x - D.x, C.y - D.y, C.z - D.z]

#     ra = [A.x, A.y, A.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [ra[n] + 2*DX[n] for n in range(3)]
#     D2 = [ra[n] + DX[n] for n in range(3)]
#     D3 = [ra[n] - DX[n] for n in range(3)]
#     D4 = [ra[n] - 2*DX[n] for n in range(3)]

#     A1 = Point(D1[0], D1[1], D1[2])
#     A2 = Point(D2[0], D2[1], D2[2])
#     A3 = Point(D3[0], D3[1], D3[2])
#     A4 = Point(D4[0], D4[1], D4[2])

#     rab_D1 = [A1.x - B.x, A1.y - B.y, A1.z - B.z]
#     rab_D2 = [A2.x - B.x, A2.y - B.y, A2.z - B.z]
#     rab_D3 = [A3.x - B.x, A3.y - B.y, A3.z - B.z]
#     rab_D4 = [A4.x - B.x, A4.y - B.y, A4.z - B.z]

#     a1 = sum([rab_D1[n]*rcb[n] for n in range(3)])/(RCB**2) - 1.
#     a2 = sum([rab_D2[n]*rcb[n] for n in range(3)])/(RCB**2) - 1.
#     a3 = sum([rab_D3[n]*rcb[n] for n in range(3)])/(RCB**2) - 1.
#     a4 = sum([rab_D4[n]*rcb[n] for n in range(3)])/(RCB**2) - 1.

#     DA_D1, dB, dC, DD_D1 = derv_dih(A1, B, C, D)
#     DA_D2, dB, dC, DD_D2 = derv_dih(A2, B, C, D)
#     DA_D3, dB, dC, DD_D3 = derv_dih(A3, B, C, D)
#     DA_D4, dB, dC, DD_D4 = derv_dih(A4, B, C, D)

#     p1 = sum([rcd[n]*rcb[n] for n in range(3)])/(RCB**2)
#     f1 = a1*DA_D1[i] - p1*DD_D1[i]
#     f2 = a2*DA_D2[i] - p1*DD_D2[i]
#     f3 = a3*DA_D3[i] - p1*DD_D3[i]
#     f4 = a4*DA_D4[i] - p1*DD_D4[i]

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_BB(i, j, A, B, C, D, dx):
#     rcd = [C.x - D.x, C.y - D.y, C.z - D.z]

#     rb = [B.x, B.y, B.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [rb[n] + 2*DX[n] for n in range(3)]
#     D2 = [rb[n] + DX[n] for n in range(3)]
#     D3 = [rb[n] - DX[n] for n in range(3)]
#     D4 = [rb[n] - 2*DX[n] for n in range(3)]

#     B1 = Point(D1[0], D1[1], D1[2])
#     B2 = Point(D2[0], D2[1], D2[2])
#     B3 = Point(D3[0], D3[1], D3[2])
#     B4 = Point(D4[0], D4[1], D4[2])

#     rab_D1 = [A.x - B1.x, A.y - B1.y, A.z - B1.z]
#     rab_D2 = [A.x - B2.x, A.y - B2.y, A.z - B2.z]
#     rab_D3 = [A.x - B3.x, A.y - B3.y, A.z - B3.z]
#     rab_D4 = [A.x - B4.x, A.y - B4.y, A.z - B4.z]

#     rcb_D1 = [C.x - B1.x, C.y - B1.y, C.z - B1.z]
#     RCB_D1 = sqrt(sum([r*r for r in rcb_D1]))
#     rcb_D2 = [C.x - B2.x, C.y - B2.y, C.z - B2.z]
#     RCB_D2 = sqrt(sum([r*r for r in rcb_D2]))
#     rcb_D3 = [C.x - B3.x, C.y - B3.y, C.z - B3.z]
#     RCB_D3 = sqrt(sum([r*r for r in rcb_D3]))
#     rcb_D4 = [C.x - B4.x, C.y - B4.y, C.z - B4.z]
#     RCB_D4 = sqrt(sum([r*r for r in rcb_D4]))

#     a1 = sum([rab_D1[n]*rcb_D1[n] for n in range(3)])/(RCB_D1**2) - 1.
#     a2 = sum([rab_D2[n]*rcb_D2[n] for n in range(3)])/(RCB_D2**2) - 1.
#     a3 = sum([rab_D3[n]*rcb_D3[n] for n in range(3)])/(RCB_D3**2) - 1.
#     a4 = sum([rab_D4[n]*rcb_D4[n] for n in range(3)])/(RCB_D4**2) - 1.

#     a5 = sum([rcb_D1[n]*rcd[n] for n in range(3)])/(RCB_D1**2)
#     a6 = sum([rcb_D2[n]*rcd[n] for n in range(3)])/(RCB_D2**2)
#     a7 = sum([rcb_D3[n]*rcd[n] for n in range(3)])/(RCB_D3**2)
#     a8 = sum([rcb_D4[n]*rcd[n] for n in range(3)])/(RCB_D4**2)

#     DA_D1, dB, dC, DD_D1 = derv_dih(A, B1, C, D)
#     DA_D2, dB, dC, DD_D2 = derv_dih(A, B2, C, D)
#     DA_D3, dB, dC, DD_D3 = derv_dih(A, B3, C, D)
#     DA_D4, dB, dC, DD_D4 = derv_dih(A, B4, C, D)

#     f1 = a1*DA_D1[i] - a5*DD_D1[i]
#     f2 = a2*DA_D2[i] - a6*DD_D2[i]
#     f3 = a3*DA_D3[i] - a7*DD_D3[i]
#     f4 = a4*DA_D4[i] - a8*DD_D4[i]

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_BC(i, j, A, B, C, D, dx):

#     rab = [A.x - B.x, A.y - B.y, A.z - B.z]

#     rc = [C.x, C.y, C.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [rc[n] + 2*DX[n] for n in range(3)]
#     D2 = [rc[n] + DX[n] for n in range(3)]
#     D3 = [rc[n] - DX[n] for n in range(3)]
#     D4 = [rc[n] - 2*DX[n] for n in range(3)]

#     C1 = Point(D1[0], D1[1], D1[2])
#     C2 = Point(D2[0], D2[1], D2[2])
#     C3 = Point(D3[0], D3[1], D3[2])
#     C4 = Point(D4[0], D4[1], D4[2])

#     rcb_D1 = [C1.x - B.x, C1.y - B.y, C1.z - B.z]
#     RCB_D1 = sqrt(sum([r*r for r in rcb_D1]))
#     rcb_D2 = [C2.x - B.x, C2.y - B.y, C2.z - B.z]
#     RCB_D2 = sqrt(sum([r*r for r in rcb_D2]))
#     rcb_D3 = [C3.x - B.x, C3.y - B.y, C3.z - B.z]
#     RCB_D3 = sqrt(sum([r*r for r in rcb_D3]))
#     rcb_D4 = [C4.x - B.x, C4.y - B.y, C4.z - B.z]
#     RCB_D4 = sqrt(sum([r*r for r in rcb_D4]))

#     rcd_D1 = [C1.x - D.x, C1.y - D.y, C1.z - D.z]
#     RCD_D1 = sqrt(sum([r*r for r in rcd_D1]))
#     rcd_D2 = [C2.x - D.x, C2.y - D.y, C2.z - D.z]
#     RCD_D2 = sqrt(sum([r*r for r in rcd_D2]))
#     rcd_D3 = [C3.x - D.x, C3.y - D.y, C3.z - D.z]
#     RCD_D3 = sqrt(sum([r*r for r in rcd_D3]))
#     rcd_D4 = [C4.x - D.x, C4.y - D.y, C4.z - D.z]
#     RCD_D4 = sqrt(sum([r*r for r in rcd_D4]))

#     a1 = sum([rab[n]*rcb_D1[n] for n in range(3)])/(RCB_D1**2) - 1.
#     a2 = sum([rab[n]*rcb_D2[n] for n in range(3)])/(RCB_D2**2) - 1.
#     a3 = sum([rab[n]*rcb_D3[n] for n in range(3)])/(RCB_D3**2) - 1.
#     a4 = sum([rab[n]*rcb_D4[n] for n in range(3)])/(RCB_D4**2) - 1.

#     a5 = sum([rcb_D1[n]*rcd_D1[n] for n in range(3)])/(RCB_D1**2)
#     a6 = sum([rcb_D2[n]*rcd_D2[n] for n in range(3)])/(RCB_D2**2)
#     a7 = sum([rcb_D3[n]*rcd_D3[n] for n in range(3)])/(RCB_D3**2)
#     a8 = sum([rcb_D4[n]*rcd_D4[n] for n in range(3)])/(RCB_D4**2)

#     DA_D1, dB, dC, DD_D1 = derv_dih(A, B, C1, D)
#     DA_D2, dB, dC, DD_D2 = derv_dih(A, B, C2, D)
#     DA_D3, dB, dC, DD_D3 = derv_dih(A, B, C3, D)
#     DA_D4, dB, dC, DD_D4 = derv_dih(A, B, C4, D)

#     f1 = a1*DA_D1[i] - a5*DD_D1[i]
#     f2 = a2*DA_D2[i] - a6*DD_D2[i]
#     f3 = a3*DA_D3[i] - a7*DD_D3[i]
#     f4 = a4*DA_D4[i] - a8*DD_D4[i]

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_BD(i, j, A, B, C, D, dx):

#     RCB = sqrt(dif(C, B))
#     rcb = [C.x - B.x, C.y - B.y, C.z - B.z]
#     rab = [A.x - B.x, A.y - B.y, A.z - B.z]

#     rd = [D.x, D.y, D.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [rd[n] + 2*DX[n] for n in range(3)]
#     D2 = [rd[n] + DX[n] for n in range(3)]
#     D3 = [rd[n] - DX[n] for n in range(3)]
#     D4 = [rd[n] - 2*DX[n] for n in range(3)]

#     d1 = Point(D1[0], D1[1], D1[2])
#     d2 = Point(D2[0], D2[1], D2[2])
#     d3 = Point(D3[0], D3[1], D3[2])
#     d4 = Point(D4[0], D4[1], D4[2])

#     rcd_D1 = [C.x - d1.x, C.y - d1.y, C.z - d1.z]
#     RCD_D1 = sqrt(sum([r*r for r in rcd_D1]))
#     rcd_D2 = [C.x - d2.x, C.y - d2.y, C.z - d2.z]
#     RCD_D2 = sqrt(sum([r*r for r in rcd_D2]))
#     rcd_D3 = [C.x - d3.x, C.y - d3.y, C.z - d3.z]
#     RCD_D3 = sqrt(sum([r*r for r in rcd_D3]))
#     rcd_D4 = [C.x - d4.x, C.y - d4.y, C.z - d4.z]
#     RCD_D4 = sqrt(sum([r*r for r in rcd_D4]))

#     a1 = sum([rcb[n]*rcd_D1[n] for n in range(3)])/(RCB**2)
#     a2 = sum([rcb[n]*rcd_D2[n] for n in range(3)])/(RCB**2)
#     a3 = sum([rcb[n]*rcd_D3[n] for n in range(3)])/(RCB**2)
#     a4 = sum([rcb[n]*rcd_D4[n] for n in range(3)])/(RCB**2)
#     p1 = sum([rab[n]*rcb[n] for n in range(3)])/(RCB**2) - 1.

#     DA_D1, dB, dC, DD_D1 = derv_dih(A, B, C, d1)
#     DA_D2, dB, dC, DD_D2 = derv_dih(A, B, C, d2)
#     DA_D3, dB, dC, DD_D3 = derv_dih(A, B, C, d3)
#     DA_D4, dB, dC, DD_D4 = derv_dih(A, B, C, d4)

#     f1 = p1*DA_D1[i] - a1*DD_D1[i]
#     f2 = p1*DA_D2[i] - a2*DD_D2[i]
#     f3 = p1*DA_D3[i] - a3*DD_D3[i]
#     f4 = p1*DA_D4[i] - a4*DD_D4[i]

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_CA(i, j, A, B, C, D, dx):

#     RCB = sqrt(dif(C, B))
#     rcb = [C.x - B.x, C.y - B.y, C.z - B.z]
#     rcd = [C.x - D.x, C.y - D.y, C.z - D.z]

#     ra = [A.x, A.y, A.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [ra[n] + 2*DX[n] for n in range(3)]
#     D2 = [ra[n] + DX[n] for n in range(3)]
#     D3 = [ra[n] - DX[n] for n in range(3)]
#     D4 = [ra[n] - 2*DX[n] for n in range(3)]

#     A1 = Point(D1[0], D1[1], D1[2])
#     A2 = Point(D2[0], D2[1], D2[2])
#     A3 = Point(D3[0], D3[1], D3[2])
#     A4 = Point(D4[0], D4[1], D4[2])

#     rab_D1 = [A1.x - B.x, A1.y - B.y, A1.z - B.z]
#     rab_D2 = [A2.x - B.x, A2.y - B.y, A2.z - B.z]
#     rab_D3 = [A3.x - B.x, A3.y - B.y, A3.z - B.z]
#     rab_D4 = [A4.x - B.x, A4.y - B.y, A4.z - B.z]

#     a1 = sum([rab_D1[n]*rcb[n] for n in range(3)])/(RCB**2)
#     a2 = sum([rab_D2[n]*rcb[n] for n in range(3)])/(RCB**2)
#     a3 = sum([rab_D3[n]*rcb[n] for n in range(3)])/(RCB**2)
#     a4 = sum([rab_D4[n]*rcb[n] for n in range(3)])/(RCB**2)

#     DA_D1, dB, dC, DD_D1 = derv_dih(A1, B, C, D)
#     DA_D2, dB, dC, DD_D2 = derv_dih(A2, B, C, D)
#     DA_D3, dB, dC, DD_D3 = derv_dih(A3, B, C, D)
#     DA_D4, dB, dC, DD_D4 = derv_dih(A4, B, C, D)

#     p1 = sum([rcd[n]*rcb[n] for n in range(3)])/(RCB**2)
#     f1 = p1*DD_D1[i] - a1*DA_D1[i]
#     f2 = p1*DD_D2[i] - a2*DA_D2[i]
#     f3 = p1*DD_D3[i] - a3*DA_D3[i]
#     f4 = p1*DD_D4[i] - a4*DA_D4[i]

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_CB(i, j, A, B, C, D, dx):
#     rcd = [C.x - D.x, C.y - D.y, C.z - D.z]

#     rb = [B.x, B.y, B.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [rb[n] + 2*DX[n] for n in range(3)]
#     D2 = [rb[n] + DX[n] for n in range(3)]
#     D3 = [rb[n] - DX[n] for n in range(3)]
#     D4 = [rb[n] - 2*DX[n] for n in range(3)]

#     B1 = Point(D1[0], D1[1], D1[2])
#     B2 = Point(D2[0], D2[1], D2[2])
#     B3 = Point(D3[0], D3[1], D3[2])
#     B4 = Point(D4[0], D4[1], D4[2])

#     rab_D1 = [A.x - B1.x, A.y - B1.y, A.z - B1.z]
#     rab_D2 = [A.x - B2.x, A.y - B2.y, A.z - B2.z]
#     rab_D3 = [A.x - B3.x, A.y - B3.y, A.z - B3.z]
#     rab_D4 = [A.x - B4.x, A.y - B4.y, A.z - B4.z]

#     rcb_D1 = [C.x - B1.x, C.y - B1.y, C.z - B1.z]
#     RCB_D1 = sqrt(sum([r*r for r in rcb_D1]))
#     rcb_D2 = [C.x - B2.x, C.y - B2.y, C.z - B2.z]
#     RCB_D2 = sqrt(sum([r*r for r in rcb_D2]))
#     rcb_D3 = [C.x - B3.x, C.y - B3.y, C.z - B3.z]
#     RCB_D3 = sqrt(sum([r*r for r in rcb_D3]))
#     rcb_D4 = [C.x - B4.x, C.y - B4.y, C.z - B4.z]
#     RCB_D4 = sqrt(sum([r*r for r in rcb_D4]))

#     a1 = sum([rab_D1[n]*rcb_D1[n] for n in range(3)])/(RCB_D1**2)
#     a2 = sum([rab_D2[n]*rcb_D2[n] for n in range(3)])/(RCB_D2**2)
#     a3 = sum([rab_D3[n]*rcb_D3[n] for n in range(3)])/(RCB_D3**2)
#     a4 = sum([rab_D4[n]*rcb_D4[n] for n in range(3)])/(RCB_D4**2)

#     a5 = sum([rcb_D1[n]*rcd[n] for n in range(3)])/(RCB_D1**2) - 1.
#     a6 = sum([rcb_D2[n]*rcd[n] for n in range(3)])/(RCB_D2**2) - 1.
#     a7 = sum([rcb_D3[n]*rcd[n] for n in range(3)])/(RCB_D3**2) - 1.
#     a8 = sum([rcb_D4[n]*rcd[n] for n in range(3)])/(RCB_D4**2) - 1.

#     DA_D1, dB, dC, DD_D1 = derv_dih(A, B1, C, D)
#     DA_D2, dB, dC, DD_D2 = derv_dih(A, B2, C, D)
#     DA_D3, dB, dC, DD_D3 = derv_dih(A, B3, C, D)
#     DA_D4, dB, dC, DD_D4 = derv_dih(A, B4, C, D)

#     f1 = a5*DD_D1[i] - a1*DA_D1[i]
#     f2 = a6*DD_D2[i] - a2*DA_D2[i]
#     f3 = a7*DD_D3[i] - a3*DA_D3[i]
#     f4 = a8*DD_D4[i] - a4*DA_D4[i]

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_CC(i, j, A, B, C, D, dx):

#     rab = [A.x - B.x, A.y - B.y, A.z - B.z]
#     rc = [C.x, C.y, C.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [rc[n] + 2*DX[n] for n in range(3)]
#     D2 = [rc[n] + DX[n] for n in range(3)]
#     D3 = [rc[n] - DX[n] for n in range(3)]
#     D4 = [rc[n] - 2*DX[n] for n in range(3)]

#     C1 = Point(D1[0], D1[1], D1[2])
#     C2 = Point(D2[0], D2[1], D2[2])
#     C3 = Point(D3[0], D3[1], D3[2])
#     C4 = Point(D4[0], D4[1], D4[2])

#     rcb_D1 = [C1.x - B.x, C1.y - B.y, C1.z - B.z]
#     RCB_D1 = sqrt(sum([r*r for r in rcb_D1]))
#     rcb_D2 = [C2.x - B.x, C2.y - B.y, C2.z - B.z]
#     RCB_D2 = sqrt(sum([r*r for r in rcb_D2]))
#     rcb_D3 = [C3.x - B.x, C3.y - B.y, C3.z - B.z]
#     RCB_D3 = sqrt(sum([r*r for r in rcb_D3]))
#     rcb_D4 = [C4.x - B.x, C4.y - B.y, C4.z - B.z]
#     RCB_D4 = sqrt(sum([r*r for r in rcb_D4]))

#     rcd_D1 = [C1.x - D.x, C1.y - D.y, C1.z - D.z]
#     rcd_D2 = [C2.x - D.x, C2.y - D.y, C2.z - D.z]
#     rcd_D3 = [C3.x - D.x, C3.y - D.y, C3.z - D.z]
#     rcd_D4 = [C4.x - D.x, C4.y - D.y, C4.z - D.z]

#     a1 = sum([rab[n]*rcb_D1[n] for n in range(3)])/(RCB_D1**2)
#     a2 = sum([rab[n]*rcb_D2[n] for n in range(3)])/(RCB_D2**2)
#     a3 = sum([rab[n]*rcb_D3[n] for n in range(3)])/(RCB_D3**2)
#     a4 = sum([rab[n]*rcb_D4[n] for n in range(3)])/(RCB_D4**2)

#     a5 = sum([rcb_D1[n]*rcd_D1[n] for n in range(3)])/(RCB_D1**2) - 1.
#     a6 = sum([rcb_D2[n]*rcd_D2[n] for n in range(3)])/(RCB_D2**2) - 1.
#     a7 = sum([rcb_D3[n]*rcd_D3[n] for n in range(3)])/(RCB_D3**2) - 1.
#     a8 = sum([rcb_D4[n]*rcd_D4[n] for n in range(3)])/(RCB_D4**2) - 1.

#     DA_D1, dB, dC, DD_D1 = derv_dih(A, B, C1, D)
#     DA_D2, dB, dC, DD_D2 = derv_dih(A, B, C2, D)
#     DA_D3, dB, dC, DD_D3 = derv_dih(A, B, C3, D)
#     DA_D4, dB, dC, DD_D4 = derv_dih(A, B, C4, D)

#     f1 = a5*DD_D1[i] - a1*DA_D1[i]
#     f2 = a6*DD_D2[i] - a2*DA_D2[i]
#     f3 = a7*DD_D3[i] - a3*DA_D3[i]
#     f4 = a8*DD_D4[i] - a4*DA_D4[i]

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_CD(i, j, A, B, C, D, dx):

#     RCB = sqrt(dif(C, B))
#     rcb = [C.x - B.x, C.y - B.y, C.z - B.z]
#     rab = [A.x - B.x, A.y - B.y, A.z - B.z]

#     rd = [D.x, D.y, D.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [rd[n] + 2*DX[n] for n in range(3)]
#     D2 = [rd[n] + DX[n] for n in range(3)]
#     D3 = [rd[n] - DX[n] for n in range(3)]
#     D4 = [rd[n] - 2*DX[n] for n in range(3)]

#     d1 = Point(D1[0], D1[1], D1[2])
#     d2 = Point(D2[0], D2[1], D2[2])
#     d3 = Point(D3[0], D3[1], D3[2])
#     d4 = Point(D4[0], D4[1], D4[2])

#     rcd_D1 = [C.x - d1.x, C.y - d1.y, C.z - d1.z]
#     RCD_D1 = sqrt(sum([r*r for r in rcd_D1]))
#     rcd_D2 = [C.x - d2.x, C.y - d2.y, C.z - d2.z]
#     RCD_D2 = sqrt(sum([r*r for r in rcd_D2]))
#     rcd_D3 = [C.x - d3.x, C.y - d3.y, C.z - d3.z]
#     RCD_D3 = sqrt(sum([r*r for r in rcd_D3]))
#     rcd_D4 = [C.x - d4.x, C.y - d4.y, C.z - d4.z]
#     RCD_D4 = sqrt(sum([r*r for r in rcd_D4]))

#     a1 = sum([rcb[n]*rcd_D1[n] for n in range(3)])/(RCB**2) - 1.
#     a2 = sum([rcb[n]*rcd_D2[n] for n in range(3)])/(RCB**2) - 1.
#     a3 = sum([rcb[n]*rcd_D3[n] for n in range(3)])/(RCB**2) - 1.
#     a4 = sum([rcb[n]*rcd_D4[n] for n in range(3)])/(RCB**2) - 1.
#     p1 = sum([rab[n]*rcb[n] for n in range(3)])/(RCB**2)

#     DA_D1, dB, dC, DD_D1 = derv_dih(A, B, C, d1)
#     DA_D2, dB, dC, DD_D2 = derv_dih(A, B, C, d2)
#     DA_D3, dB, dC, DD_D3 = derv_dih(A, B, C, d3)
#     DA_D4, dB, dC, DD_D4 = derv_dih(A, B, C, d4)

#     f1 = a1*DD_D1[i] - p1*DA_D1[i]
#     f2 = a2*DD_D2[i] - p1*DA_D2[i]
#     f3 = a3*DD_D3[i] - p1*DA_D3[i]
#     f4 = a4*DD_D4[i] - p1*DA_D4[i]

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_DB(i, j, B, C, D, dx):
#     rcd = [C.x - D.x, C.y - D.y, C.z - D.z]

#     rb = [B.x, B.y, B.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [rb[n] + 2*DX[n] for n in range(3)]
#     D2 = [rb[n] + DX[n] for n in range(3)]
#     D3 = [rb[n] - DX[n] for n in range(3)]
#     D4 = [rb[n] - 2*DX[n] for n in range(3)]

#     B1 = Point(D1[0], D1[1], D1[2])
#     B2 = Point(D2[0], D2[1], D2[2])
#     B3 = Point(D3[0], D3[1], D3[2])
#     B4 = Point(D4[0], D4[1], D4[2])

#     rcb_D1 = [C.x - B1.x, C.y - B1.y, C.z - B1.z]
#     RCB_D1 = sqrt(sum([r*r for r in rcb_D1]))
#     rcb_D2 = [C.x - B2.x, C.y - B2.y, C.z - B2.z]
#     RCB_D2 = sqrt(sum([r*r for r in rcb_D2]))
#     rcb_D3 = [C.x - B3.x, C.y - B3.y, C.z - B3.z]
#     RCB_D3 = sqrt(sum([r*r for r in rcb_D3]))
#     rcb_D4 = [C.x - B4.x, C.y - B4.y, C.z - B4.z]
#     RCB_D4 = sqrt(sum([r*r for r in rcb_D4]))

#     rnc_D1 = cross(rcb_D1, rcd)
#     RNC_D1 = sqrt(sum([r*r for r in rnc_D1]))
#     rnc_D2 = cross(rcb_D2, rcd)
#     RNC_D2 = sqrt(sum([r*r for r in rnc_D2]))
#     rnc_D3 = cross(rcb_D3, rcd)
#     RNC_D3 = sqrt(sum([r*r for r in rnc_D3]))
#     rnc_D4 = cross(rcb_D4, rcd)
#     RNC_D4 = sqrt(sum([r*r for r in rnc_D4]))

#     f1 = -RCB_D1*rnc_D1[i]/(RNC_D1**2)
#     f2 = -RCB_D2*rnc_D2[i]/(RNC_D2**2)
#     f3 = -RCB_D3*rnc_D3[i]/(RNC_D3**2)
#     f4 = -RCB_D4*rnc_D4[i]/(RNC_D4**2)

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_DC(i, j, B, C, D, dx):
#     rcd = [C.x - D.x, C.y - D.y, C.z - D.z]

#     rc = [C.x, C.y, C.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [rc[n] + 2*DX[n] for n in range(3)]
#     D2 = [rc[n] + DX[n] for n in range(3)]
#     D3 = [rc[n] - DX[n] for n in range(3)]
#     D4 = [rc[n] - 2*DX[n] for n in range(3)]

#     C1 = Point(D1[0], D1[1], D1[2])
#     C2 = Point(D2[0], D2[1], D2[2])
#     C3 = Point(D3[0], D3[1], D3[2])
#     C4 = Point(D4[0], D4[1], D4[2])

#     rcb_D1 = [C1.x - B.x, C1.y - B.y, C1.z - B.z]
#     RCB_D1 = sqrt(sum([r*r for r in rcb_D1]))
#     rcb_D2 = [C2.x - B.x, C2.y - B.y, C2.z - B.z]
#     RCB_D2 = sqrt(sum([r*r for r in rcb_D2]))
#     rcb_D3 = [C3.x - B.x, C3.y - B.y, C3.z - B.z]
#     RCB_D3 = sqrt(sum([r*r for r in rcb_D3]))
#     rcb_D4 = [C4.x - B.x, C4.y - B.y, C4.z - B.z]
#     RCB_D4 = sqrt(sum([r*r for r in rcb_D4]))

#     rcd_D1 = [C1.x - D.x, C1.y - D.y, C1.z - D.z]
#     rcd_D2 = [C2.x - D.x, C2.y - D.y, C2.z - D.z]
#     rcd_D3 = [C3.x - D.x, C3.y - D.y, C3.z - D.z]
#     rcd_D4 = [C4.x - D.x, C4.y - D.y, C4.z - D.z]

#     rnc_D1 = cross(rcb_D1, rcd_D1)
#     RNC_D1 = sqrt(sum([r*r for r in rnc_D1]))
#     rnc_D2 = cross(rcb_D2, rcd_D2)
#     RNC_D2 = sqrt(sum([r*r for r in rnc_D2]))
#     rnc_D3 = cross(rcb_D3, rcd_D3)
#     RNC_D3 = sqrt(sum([r*r for r in rnc_D3]))
#     rnc_D4 = cross(rcb_D4, rcd_D4)
#     RNC_D4 = sqrt(sum([r*r for r in rnc_D4]))

#     f1 = -RCB_D1*rnc_D1[i]/(RNC_D1**2)
#     f2 = -RCB_D2*rnc_D2[i]/(RNC_D2**2)
#     f3 = -RCB_D3*rnc_D3[i]/(RNC_D3**2)
#     f4 = -RCB_D4*rnc_D4[i]/(RNC_D4**2)

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def numerical_derv_DD(i, j, B, C, D, dx):

#     RCB = sqrt(dif(C, B))
#     rcb = [C.x - B.x, C.y - B.y, C.z - B.z]

#     rd = [D.x, D.y, D.z]
#     DX = [0., 0., 0.]
#     for m in range(3):
#         if m == j:
#             DX[j] = dx

#     D1 = [rd[n] + 2*DX[n] for n in range(3)]
#     D2 = [rd[n] + DX[n] for n in range(3)]
#     D3 = [rd[n] - DX[n] for n in range(3)]
#     D4 = [rd[n] - 2*DX[n] for n in range(3)]

#     d1 = Point(D1[0], D1[1], D1[2])
#     d2 = Point(D2[0], D2[1], D2[2])
#     d3 = Point(D3[0], D3[1], D3[2])
#     d4 = Point(D4[0], D4[1], D4[2])

#     rcd_D1 = [C.x - d1.x, C.y - d1.y, C.z - d1.z]
#     rcd_D2 = [C.x - d2.x, C.y - d2.y, C.z - d2.z]
#     rcd_D3 = [C.x - d3.x, C.y - d3.y, C.z - d3.z]
#     rcd_D4 = [C.x - d4.x, C.y - d4.y, C.z - d4.z]

#     rnc_D1 = cross(rcb, rcd_D1)
#     RNC_D1 = sqrt(sum([r*r for r in rnc_D1]))
#     rnc_D2 = cross(rcb, rcd_D2)
#     RNC_D2 = sqrt(sum([r*r for r in rnc_D2]))
#     rnc_D3 = cross(rcb, rcd_D3)
#     RNC_D3 = sqrt(sum([r*r for r in rnc_D3]))
#     rnc_D4 = cross(rcb, rcd_D4)
#     RNC_D4 = sqrt(sum([r*r for r in rnc_D4]))

#     f1 = -RCB*rnc_D1[i]/(RNC_D1**2)
#     f2 = -RCB*rnc_D2[i]/(RNC_D2**2)
#     f3 = -RCB*rnc_D3[i]/(RNC_D3**2)
#     f4 = -RCB*rnc_D4[i]/(RNC_D4**2)

#     derv = (-f1 + 8*f2 - 8*f3 + f4)/(12.*dx)

#     return derv


# def segunda_wilson(mol, ndim, numat, bond, ang, dih):
#     seg_wilson = np.zeros((ndim, 3*numat, 3*numat), dtype=float)

#     for k in range(len(bond)):
#         m1 = bond[k][0]
#         n1 = 3*(m1-1)
#         m2 = bond[k][1]
#         n2 = 3*(m2-1)
#         for i in range(3):
#             for j in range(3):
#                 seg_wilson[k, n1+i, n1 +
#                            j] = seg_derv_bond(mol[m1-1].cart, mol[m2-1].cart, i, j)
#                 seg_wilson[k, n1+i, n2+j] = -1 * \
#                     seg_derv_bond(mol[m1-1].cart, mol[m2-1].cart, i, j)
#                 seg_wilson[k, n2+j, n1+i] = seg_wilson[k, n1+i, n2+j]
#                 seg_wilson[k, n2+j, n2+i] = seg_wilson[k, n1+i, n1+j]

#     for k in range(len(ang)):
#         m1 = ang[k][0]
#         m2 = ang[k][1]
#         m3 = ang[k][2]
#         nq = k + len(bond)
#         n1 = 3*(m1-1)
#         n2 = 3*(m2-1)
#         n3 = 3*(m3-1)
#         DA, DB, DC = derv_ang(mol[m1-1].cart, mol[m2-1].cart, mol[m3-1].cart)
#         for i in range(3):
#             for j in range(3):
#                 part1, part2, part3 = seg_derv_ang(
#                     mol[m1-1].cart, mol[m2-1].cart, mol[m3-1].cart, DA, DB, DC, i, j)

#                 seg_wilson[nq, n1+i, n1+j], seg_wilson[nq, n1 +
#                                                        i, n2+j], seg_wilson[nq, n1+i, n3+j] = part1
#                 seg_wilson[nq, n2+i, n1+j], seg_wilson[nq, n2 +
#                                                        i, n2+j], seg_wilson[nq, n2+i, n3+j] = part2
#                 seg_wilson[nq, n3+i, n1+j], seg_wilson[nq, n3 +
#                                                        i, n2+j], seg_wilson[nq, n3+i, n3+j] = part3

#     for k in range(len(dih)):
#         m1 = dih[k][0]
#         m2 = dih[k][1]
#         m3 = dih[k][2]
#         m4 = dih[k][3]
#         nq = k+len(bond)+len(ang)
#         n1 = 3*(m1-1)
#         n2 = 3*(m2-1)
#         n3 = 3*(m3-1)
#         n4 = 3*(m4-1)

#         A, B, C, D = mol[m1-1].cart, mol[m2 -
#                                          1].cart, mol[m3-1].cart, mol[m4-1].cart
#         RAB = sqrt(dif(A, B))
#         RCB = sqrt(dif(C, B))
#         RDC = sqrt(dif(C, D))

#         rab = [A.x - B.x, A.y - B.y, A.z - B.z]
#         rcb = [C.x - B.x, C.y - B.y, C.z - B.z]
#         rdc = [D.x - C.x, D.y - C.y, D.z - C.z]

#         Nrab = [rab[x]/RAB for x in range(3)]
#         Nrcb = [rcb[x]/RCB for x in range(3)]
#         Nrdc = [rdc[x]/RDC for x in range(3)]

#         uxw = cross(Nrab, Nrcb)
#         RMB = sqrt(sum([r*r for r in uxw]))
#         vxw = cross(Nrdc, Nrcb)
#         RNC = sqrt(sum([r*r for r in vxw]))

#         teta1 = angulos(A, B, C)
#         teta2 = angulos(B, C, D)

#         a1, a2 = sin(teta1), sin(teta2)
#         b1, b2 = cos(teta1), cos(teta2)

#         dx = 0.0001

#         for i in range(3):
#             for j in range(3):
#                 seg_wilson[nq, n1+i, n1 +
#                            j] = numerical_derv_AA(i, j, A, B, C, dx)
#                 seg_wilson[nq, n1+i, n2 +
#                            j] = numerical_derv_AB(i, j, A, B, C, dx)
#                 seg_wilson[nq, n1+i, n3 +
#                            j] = numerical_derv_AC(i, j, A, B, C, dx)
#                 seg_wilson[nq, n1+i, n4+j] = 0.

#                 seg_wilson[nq, n2+j, n1 +
#                            i] = numerical_derv_BA(j, i, A, B, C, D, dx)
#                 seg_wilson[nq, n2+i, n2 +
#                            j] = numerical_derv_BB(i, j, A, B, C, D, dx)
#                 seg_wilson[nq, n2+i, n3 +
#                            j] = numerical_derv_BC(i, j, A, B, C, D, dx)
#                 seg_wilson[nq, n2+i, n4 +
#                            j] = numerical_derv_BD(i, j, A, B, C, D, dx)

#                 seg_wilson[nq, n3+j, n1 +
#                            i] = numerical_derv_CA(j, i, A, B, C, D, dx)
#                 seg_wilson[nq, n3+j, n2 +
#                            i] = numerical_derv_CB(j, i, A, B, C, D, dx)
#                 seg_wilson[nq, n3+i, n3 +
#                            j] = numerical_derv_CC(i, j, A, B, C, D, dx)
#                 seg_wilson[nq, n3+i, n4 +
#                            j] = numerical_derv_CD(i, j, A, B, C, D, dx)

#                 seg_wilson[nq, n4+j, n1+i] = 0.
#                 seg_wilson[nq, n4+j, n2 +
#                            i] = numerical_derv_DB(j, i, B, C, D, dx)
#                 seg_wilson[nq, n4+j, n3 +
#                            i] = numerical_derv_DC(j, i, B, C, D, dx)
#                 seg_wilson[nq, n4+i, n4 +
#                            j] = numerical_derv_DD(i, j, B, C, D, dx)

#     for k in range(ndim):
#         for i in range(3*numat):
#             for j in range(3*numat):
#                 w = seg_wilson[k, i, j] - seg_wilson[k, j, i]
#                 if w > 1.e-2:
#                     print(seg_wilson[k, i, j], seg_wilson[k, j, i], w, i, j, k)

#     return seg_wilson


# def triplet_pot(gap, ro, r1, gradient, Hessian):

#     Q = np.array(([r1[i] - ro[i] for i in range(len(r1))]), dtype=float)
#     product1 = np.dot(Q, gradient)
#     product2 = np.dot(Q, np.dot(Hessian, Q))

#     energy = gap + product1 + product2

#     return energy

# # Functions of the conectivity


# def mtx_conect(file):

#     conex = list(open(file))
#     numat = int(conex[0])
#     A = [tuple([int(m) for m in w.split()]) for w in conex[1:]]
#     M = np.zeros((numat, numat), dtype=int)
#     for val in A:
#         i = val[0]-1
#         j = val[1]-1
#         M[i, j] = 1
#         M[j, i] = 1

#     ang = []
#     dih = []
#     for i in xrange(numat):
#         for j in xrange(numat):
#             r = max(i+1, j+1)
#             for k in xrange(r, numat):
#                 x = M[i, j]
#                 y = M[j, k]
#                 w = x*y
#                 if w != 0 and i != j and j != k and i != k:
#                     tup = (i+1, j+1, k+1)
#                     ang.append(tup)
#                 for m in xrange(k+1, numat):
#                     z1 = w * M[i, m]
#                     z2 = w * M[j, m]
#                     z3 = w * M[k, m]
#                     if z1 == 1:
#                         d = (m+1, i+1, j+1, k+1)
#                         dih.append(d)
#                     elif z2 == 1:
#                         d = (i+1, j+1, m+1, k+1)
#                         dih.append(d)
#                     elif z3 == 1:
#                         d = (i+1, j+1, k+1, m+1)
#                         dih.append(d)

#     internas = open('internas.dat', 'w')
#     internas.write('%d' % len(A))
#     internas.write('\n')
#     for i in xrange(len(A)):
#         internas.write('%d  %d' % (A[i][0], A[i][1]))
#         internas.write('\n')
#     internas.write('%d' % len(ang))
#     internas.write('\n')
#     for j in xrange(len(ang)):
#         internas.write('%d  %d  %d' % (ang[j][0], ang[j][1], ang[j][2]))
#         internas.write('\n')
#     internas.write('%d' % len(dih))
#     internas.write('\n')
#     for k in xrange(len(dih)):
#         internas.write('%d  %d  %d  %d' %
#                        (dih[k][0], dih[k][1], dih[k][2], dih[k][3]))
#         internas.write('\n')


# # new is a list of tuples after removed the dihedral angles that
# # have indefinitions

# def resize(problems, conex_dih):
#     a = [list(n) for n in problems]
#     z = list()
#     for trio in a:
#         for w in conex_dih:
#             val = list(w)
#             if trio[0] in val and trio[1] in val and trio[2] in val and w not in z:
#                 z.append(w)
#     new = [val for val in conex_dih if val not in z]

#     return new


# def internas_resize(conex_bond, conex_ang, conex_dih):

#     internas = open('internas.dat', 'w')
#     internas.write('%d' % len(conex_bond))
#     internas.write('\n')
#     for i in xrange(len(conex_bond)):
#         internas.write('%d  %d' % (conex_bond[i][0], conex_bond[i][1]))
#         internas.write('\n')
#     internas.write('%d' % len(conex_ang))
#     internas.write('\n')
#     for j in xrange(len(conex_ang)):
#         internas.write('%d  %d  %d' %
#                        (conex_ang[j][0], conex_ang[j][1], conex_ang[j][2]))
#         internas.write('\n')
#     internas.write('%d' % len(conex_dih))
#     internas.write('\n')
#     for k in xrange(len(conex_dih)):
#         internas.write('%d  %d  %d  %d' % (
#             conex_dih[k][0], conex_dih[k][1], conex_dih[k][2], conex_dih[k][3]))
#         internas.write('\n')


# def write_freq(freq, nbond, nang, ndih, (floor, roof), m):
#     nombre = 'histograma' + '_' + \
#         str(float(floor)) + '_' + str(float(roof)) + '.out'
#     histo = open(nombre, 'w')
#     histo.write(
#         ' the interval of energy is [%g,%g] in kcal/mol' % (floor, roof))
#     histo.write('\n')
#     histo.write('First is printed the internal coordinate number, then the middle point of the %d slices is printed with \
#     the corresponding frecuencies for that slice' % len(freq[0]))
#     histo.write('\n')
#     cont = 0
#     for d in freq:
#         if cont < nbond:
#             histo.write('bond number %d' % cont)
#             histo.write('\n')
#             lista = [key for key in d]
#             lista.sort()
#             for val in lista:
#                 mid = (val[0]+val[1])/2.
#                 histo.write('%g %d' % (mid, freq[cont][val]))
#                 histo.write('\n')

#         histo.write('\n')
#         if nbond <= cont and cont < nang+nbond:
#             number = 1+cont-nbond
#             histo.write('angle number %d' % (number))
#             histo.write('\n')
#             lista = [key for key in d]
#             lista.sort()
#             for val in lista:
#                 mid = (val[0]+val[1])/2.
#                 histo.write('%g %d' % (mid, freq[cont][val]))
#                 histo.write('\n')

#         histo.write('\n')
#         if nang+nbond <= cont:
#             number = 1+cont-nbond-nang
#             histo.write('dihedral number %d' % (number))
#             histo.write('\n')
#             lista = [key for key in d]
#             lista.sort()
#             for val in lista:
#                 mid = (val[0]+val[1])/2.
#                 histo.write('%g %d' % (mid, freq[cont][val]))
#                 histo.write('\n')

#         cont += 1
