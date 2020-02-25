# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:23:39 2020

@author: Gebruiker
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 15:10:10 2020

@author: Shivan
"""
from matplotlib import pyplot as plt
import numpy as np
import math
#-------------------------------------------------------------------------
# Input paramters aileron
#-------------------------------------------------------------------------

C_a = 0.547              #chord length aileron [m]
l_a = 2.771              #span of the aileron [m]
x1 = 0.153               #x-location of hing 1 [m]
x2 = 1.281               #x-location of hing 2 [m]
x3 = 2.681               #x-location of hing 3 [m]
x_a = 0.28               #distance between actuator 1 and 2 [m]
h = 0.225/2              #aileron height [m]
t_sk = 0.0011            #skin thickness [m]
t_sp = 0.0029            #spar thickness [m]
t_st = 0.0012            #stiffener thickness [m]
h_st = 0.015             #heigt of stiffener [m]
w_st = 0.020             #width of stiffener [m]
n_st = 17                #number of stiffeners equally spaced [-]
d1 = 0.01103             #vertical displacement hinge 1
d3 = 0.01642             #vertical displacement hinge 3
theta = np.radians(26)   #maximum upward deflection
P = 91700                #load in actuator 2


#Calculation of stiffener spacing
"""This function first computes the perimeter of the aileron
    Than the constant stiffener spacing is calculated by dividing the perimeter by the total number of stringers"""
def Stiffener_spacing(h,C_a):
    perimeter_semicircle = np.pi*h
    perimeter_triangle = 2*np.sqrt((C_a-h)**2+h**2)
    stiffener_spacing_semicircle = perimeter_semicircle/6
    stiffener_spacing_triangle = perimeter_triangle/14
    return stiffener_spacing_triangle, stiffener_spacing_semicircle, perimeter_semicircle,perimeter_triangle       #[m]

#calculation of stiffener locations
def Stiffener_locations(h,C_a,stiffener_spacing_semicircle):
    y_lst = []
    z_lst = []

    for i in range(3):
        theta = np.radians(i * stiffener_spacing_semicircle/(np.pi*h)*180)
        delta_y = h*np.sin(theta)
        delta_z = h-h*np.cos(theta)
        y_lst.append(delta_y)
        z_lst.append(-delta_z)
    for i in range(1,7):
        y = h - i*h/7
        z = -h-i*(C_a-h)/7
        y_lst.append(y)
        z_lst.append(z)
    reverse_z = z_lst[::-1]
    reverse_y = y_lst[::-1]
    for i in range(len(reverse_z)-1):
       z_lst.append(reverse_z[i])
       y_lst.append(-reverse_y[i])
    return (z_lst, y_lst)

delta_st_triangle, delta_st_semicircle,per_semicircle,per_triangle = Stiffener_spacing(h,C_a)          # [m]
z_lst, y_lst = Stiffener_locations(h,C_a,delta_st_semicircle)     # [m]

#stiffener area calculation
def Stiffener_A(w_st,h_st,t_st):
    A_st = w_st*t_st+(h_st-t_st)*t_st
    return A_st
A_st=Stiffener_A(w_st,h_st,t_st)

#centroid

def Centroid(z_lst,A_st,C_a,t_sk,t_sp,h):
    yc=0
    stff = sum([i * A_st for i in z_lst])
    triangle = -2*(h+(C_a-h)/2)*np.sqrt((C_a-h)**2+h**2)*t_sk
    semic =  -h*(1-2/np.pi)*np.pi*2*h*t_sk/2
    spar=-t_sp*h*2*h

    zc=(stff+semic+triangle+spar)/(17*A_st+np.sqrt((C_a-h)**2+h**2)*t_sk*2+np.pi*2*h*t_sk/2+t_sp*2*h)
    return (yc, zc)

yc,zc=Centroid(z_lst,A_st,C_a,t_sk,t_sp,h)


#moment of inertia calculation
def MMI(z_lst,y_lst,A_st,C_a,t_sk,t_sp,h,zc,yc):
    incl = np.sqrt((C_a - h) ** 2 + h ** 2)

    stff_z= sum([(i - yc) ** 2 * A_st for i in y_lst])
    triangle_z = ((2*incl)**3*t_sk/12)*((h)/incl)**2
    spar_z = (2*t_sp/12)*(2*h)**3
    semic_z = 0.5*np.pi*h**3*t_sk/8
    Izz= stff_z+triangle_z+spar_z+semic_z

    stff_y = sum([(i - zc) ** 2 * A_st for i in z_lst])
    triangle_y = 2 / 12 * t_sk * (incl) ** 3 * ((C_a - h) / incl) ** 2 + 2 * incl * t_sk * (-h - (C_a - h) / 2 - zc) ** 2
    spar_y = (2 * h / 12) * (t_sk) ** 3 + 2 * h * t_sp * (-h - zc) ** 2
    semic_y = h ** 3 * t_sk * (np.pi / 2 - 4 / np.pi) + np.pi * h * t_sk * ((-h + 2 * h / np.pi) - zc) ** 2

    Iyy= stff_y + spar_y + semic_y +triangle_y
    #Ixx = 1/12*l_a*C_a**3

    return Izz,Iyy

Izz,Iyy=MMI(z_lst,y_lst,A_st,C_a,t_sk,t_sp,h,zc,yc)

#          Integration tool--------------------checked
def integrate(a,b,n,fx): #enter funtion in the form lamda x:f(x)

    f =  fx
    h = (b-a)/n
    A = 0.5*h*(f(a)+f(b))
    for i in range(1,n):
        A += h*(f(a+i*h))
    return A

#        shear flow calculated using integration
#        a and b are the limits of integration
incl = np.sqrt((C_a - h) ** 2 + h ** 2)

def shearfl(Iyy,Izz,A_st,y_lst,z_lst,zc,h,t_sk,C_a,segment,b,Sy,Sz,n,incl):  #see diagram on paper

     qbl=[]
     qb1t = (-Sy / Izz) * (integrate(0, 0.5 * np.pi, 100, lambda x: h ** 2 * t_sk * np.sin(x)) + sum([A_st * i for i in (y_lst[0:3])]))+\
            (-Sz / Iyy) * (integrate(0, 0.5 * np.pi, 100, lambda x: (h-h*np.cos(x)-zc)*h*t_sk) + sum([A_st * (i-zc) for i in (z_lst[0:3])]))
     qb2t = (-Sy / Izz) * (integrate(0, h, 100, lambda x: t_sp * x)) + (-Sz / Iyy) * (integrate(0, h, 100, lambda x: t_sp * (-h-zc)))
     qb3t = (-Sy / Izz) * (integrate(0, incl, 100, lambda x: (h - (h / incl) * x)*t_sk) + sum( [A_st * i for i in (y_lst[3:9])])) + qb1t + qb2t +\
            (-Sz / Iyy) * (integrate(0, incl, 100, lambda x: (-h-((C_a-h)/incl)*x-zc)*t_sk) + sum( [A_st * (i-zc) for i in (z_lst[3:9])]))
     qb4t = (-Sy / Izz) * (integrate(0, incl, 100, lambda x: (- (h / incl) * x)*t_sk)+ sum([A_st * i for i in (y_lst[9:15])])) + qb3t +\
            (-Sz / Iyy) * (integrate(0, incl, 100, lambda x: (-C_a+((C_a-h)/incl)*x-zc)*t_sk)+ sum([A_st * (i-zc) for i in (z_lst[9:15])]))
     qb5t = (-Sy / Izz) * (integrate(0, -h, 100, lambda x: t_sp * x)) + (-Sz / Iyy) * (integrate(0, -h, 100, lambda x: t_sp * (-h-zc)))

     #qb6t = (-Sy / Izz) * (integrate(-0.5*np.pi, 0, 100, lambda x: h ** 2 * t_sk * np.sin(x)) + sum([A_st * i for i in (y_lst[15:17])])) + qb4t - qb5t

     h_seg=b/((len(segment)*n)-1)

     ds=0
     for segment in segment:
           for i in range(n):
                b=h_seg*(i+ds)
                # upper part semi circle
                if segment <=2:
                          qb=(-Sy/Izz)*(integrate(0,b,100, lambda x:h**2*t_sk*np.sin(x))+sum([A_st*i for i in (y_lst[0:segment+1])]))+\
                             (-Sz / Iyy) * (integrate(0,b, 100, lambda x: (h-h*np.cos(x)-zc)*h*t_sk) + sum([A_st * (i-zc) for i in (z_lst[0:segment+1])]))
                # upper spar part
                elif segment ==3:
                          qb = (-Sy / Izz) * (integrate(0, b, 100, lambda x: t_sp * x )) + (-Sz / Iyy) * (integrate(0, b, 100, lambda x: t_sp * (-h-zc)))
                # upper triangular part
                elif 4<=segment<=10:
                    if segment == 4:
                          qb = (-Sy / Izz) * (integrate(0, b, 100, lambda x: (h - (h / incl) * x)*t_sk)) + qb1t + qb2t +\
                               (-Sz / Iyy) * (integrate(0, b, 100, lambda x: ((-h-((C_a-h)/incl)*x)-zc)*t_sk))
                    else:
                          qb = (-Sy / Izz) * (integrate(0, b, 100, lambda x: (h - (h / incl) * x)*t_sk) + sum([A_st * i for i in (y_lst[3:segment-1])])) + qb1t + qb2t +\
                               (-Sz / Iyy) * (integrate(0, b, 100, lambda x: ((-h-((C_a-h)/incl)*x)-zc)*t_sk) + sum([A_st * (i-zc) for i in (z_lst[3:segment-1])]))
                # lower triangular part
                elif 11<=segment<=17:
                    if segment ==11:
                          qb = (-Sy / Izz) * (integrate(0, b, 100, lambda x: - (h / incl) * x*t_sk)) + qb3t +\
                               (-Sz / Iyy) * (integrate(0, b, 100, lambda x: ((-C_a+((C_a-h)/incl)*x)-zc)*t_sk))
                    else:
                          qb = (-Sy / Izz) * (integrate(0, b, 100, lambda x: - (h / incl) * x*t_sk) + sum([A_st * i for i in (y_lst[9:segment-2])])) + qb3t +\
                               (-Sz / Iyy) * (integrate(0, b, 100, lambda x: ((-C_a+((C_a-h)/incl)*x)-zc)*t_sk)+ sum([A_st * (i-zc) for i in (z_lst[9:segment-2])]))
                # lower spar part
                elif segment == 18:
                          qb   = (-Sy / Izz) * (integrate(0, -b, 100, lambda x: t_sp * x)) + (-Sz / Iyy) * (integrate(0, -b, 100, lambda x: t_sp * (-h-zc)))
                # lower part semi circle
                else:
                    b=-0.5*np.pi+b
                    if segment==19:
                          qb  = (-Sy / Izz) * (integrate(-0.5*np.pi,b, 100, lambda x: h ** 2 * t_sk * np.sin(x)) ) + qb4t -qb5t +\
                                (-Sz / Iyy) * (integrate(-0.5*np.pi,b, 100, lambda x: ((h-h*np.cos(x))-zc)*h*t_sk  ))
                    else:
                          qb = (-Sy / Izz) * (integrate(-0.5*np.pi,b, 100, lambda x: h ** 2 * t_sk * np.sin(x)) + sum([A_st * i for i in (y_lst[15:segment-4])])) + qb4t - qb5t +\
                               (-Sz / Iyy) * (integrate(-0.5*np.pi,b, 100, lambda x: ((h-h*np.cos(x))-zc)*h*t_sk)+sum([A_st * (i-zc) for i in (z_lst[15:segment-4])]))
                qbl.append(qb)
           ds+=n
           #print(qb5t)
     return qbl
n1=20

#Sy=int(input("Is Sy 1 or 0?"))
#Sz=int(input("Is Sz 1 or 0?"))
Sz=0
Sy=1

qb1 = shearfl(Iyy,Izz,A_st,y_lst,z_lst,zc,h,t_sk,C_a,[0,1,2],0.5*np.pi,Sy,Sz,n1,incl)
qb2 = shearfl(Iyy,Izz,A_st,y_lst,z_lst,zc,h,t_sk,C_a,[3],h,Sy,Sz,2*n1,incl)
qb3 = shearfl(Iyy,Izz,A_st,y_lst,z_lst,zc,h,t_sk,C_a,[4,5,6,7,8,9,10],incl,Sy,Sz,n1,incl)
qb4 = shearfl(Iyy,Izz,A_st,y_lst,z_lst,zc,h,t_sk,C_a,[11,12,13,14,15,16,17],incl,Sy,Sz,n1,incl)
qb5 = shearfl(Iyy,Izz,A_st,y_lst,z_lst,zc,h,t_sk,C_a,[18],h,Sy,Sz,2*n1,incl)
qb6 = shearfl(Iyy,Izz,A_st,y_lst,z_lst,zc,h,t_sk,C_a,[19,20,21],0.5*np.pi,Sy,Sz,n1,incl)


# points along line of integration (s)
s_co1 = np.linspace(0,np.pi*h/2,3*n1)
s_co2 = np.linspace(0,h,2*n1)
s_co3 = np.linspace(0,incl,7*n1)
s_co4 = s_co3
s_co5 = s_co2
s_co6 = s_co1

                                      #   z and y coordinate calculation
ds = 0
segment=[0,1,2]
zco1=[]
yco1=[]
h_seg=0.5*np.pi/((len(segment)*n1)-1)
for segment in segment:
    for i in range(n1):
        b = h_seg * (i + ds)
        zco1.append(-(h-h*np.cos(b)))
        yco1.append(h*np.sin(b))
    ds+=n1

zco6= zco1[::-1]
yco6= [-i for i in yco1[::-1]]

zco2= 2*n1*[-h]
yco2= s_co2
zco5= zco2
yco5= np.linspace(0,-h,2*n1)
zco3= np.linspace(-h,-C_a,7*n1)
yco3= np.linspace(h,0,7*n1)
zco4= zco3[::-1]
yco4= [-i for i in yco3[::-1]]


def second_integration(z_lst,q_lst):            #integration over number of points with known value,, approximated by trapezoidal rule, do for each wall!
                                                #q_lst is the shear flow of a number of points at a specific wall,
    step= z_lst[1] - z_lst[0]
    A = 0.5 * step * (q_lst[0] + q_lst[-1])
    for i in range(1,len(z_lst)):
        A += step*(q_lst[i])
    return A

#testing function
#A=second_integration([np.radians(i) for i in np.linspace(0,90,1000)],[np.sin(np.radians(i)) for i in np.linspace(0,90,1000)])

qblist1 =[[i/t_sk for i in qb1],[-i/t_sp for i in qb2],[-i/t_sp for i in qb5],[i/t_sk for i in qb6]]
zlist1  =[s_co1,s_co2,s_co5,s_co6]

qblist2= [[i/t_sp for i in qb2],[i/t_sp for i in qb5],[i/t_sk for i in qb3],[i/t_sk for i in qb4]]
zlist2= [s_co2,s_co5,s_co3,s_co4]

def second_int(zlist,qblist):
    A=0
    for i,j in zip(zlist,qblist):
        A+=second_integration(i,j)
    return A

a3=second_int(zlist1,qblist1)
b3=second_int(zlist2,qblist2)
a1= np.pi*h/t_sk + 2*h/t_sp
a2= -2*h/t_sp
b1= -2*h/t_sp
b2= 2*(incl/t_sk) + 2*h/t_sp

A = np.array([[a1,a2],[b1,b2]])
g = np.array([-a3,-b3])
solution = np.linalg.solve(A,g)

qs01 = solution[0]
qs02 = solution[1]
#print(qs01,qs02)
#print(np.allclose(np.dot(A,solution), g))
#qs01=1.283
#qs02=4.865

def Add_constant_shearflow(qb1,qb2,qb3,qb4,qb5,qb6,qs01,qs02):
     qb1 = [i + qs01 for i in qb1]
     qb6 = [i + qs01 for i in qb6]
     qb2 = [i - qs01 + qs02 for i in qb2]
     qb5 = [i - qs01 + qs02 for i in qb5]
     qb4 = [i + qs02 for i in qb4]
     qb3 = [i + qs02 for i in qb3]
     return qb1,qb2,qb3,qb4,qb5,qb6

qb1,qb2,qb3,qb4,qb5,qb6=Add_constant_shearflow(qb1,qb2,qb3,qb4,qb5,qb6,qs01,qs02)

def Cell_area(h,C_a):
    A1 = np.pi*h**2/2
    A2 = h*(C_a-h)
    return A1, A2   #m^2
A1,A2=Cell_area(h,C_a)

# -----------------------------------------------------  Moment around mid spar and shear center ---------------------------------------------------

M2 =  ((C_a-h) / incl)*h*second_integration(s_co3,qb3) + ((C_a-h) / incl)*h*second_integration(s_co4,qb4) + 2*h*second_integration(s_co1,qb1)

#if (Sy==1 and Sz==0):
 #print("Shear center location from the leading edge: ",-h-M2)

# ------------------------------------------------------------------------------------------------ ---------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------------------------------
#Calculate the shear stress due to shear

def Shear_stress_due_to_shear(qb1,qb2,qb3,qb4,qb5,qb6,t_sk,t_sp):        #compute shear stress in skins and spar by dividing by thickness, note spar shear flow not included in shear_flow_lst
    #Shear stress distribution per section due to Sy and Sz computed at the middle of each section
    tau_1 = [qb1[i]/t_sk for i in range(len(qb1))]
    tau_2 = [qb2[i]/t_sp for i in range(len(qb2))]
    tau_3 = [qb3[i]/t_sk for i in range(len(qb3))]
    tau_4 = [qb4[i]/t_sk for i in range(len(qb4))]
    tau_5 = [qb5[i]/t_sp for i in range(len(qb5))]
    tau_6 = [qb6[i]/t_sk for i in range(len(qb6))]                           
    return tau_1, tau_2, tau_3, tau_4, tau_5, tau_6

tau_1, tau_2, tau_3, tau_4, tau_5, tau_6 = Shear_stress_due_to_shear(qb1, qb2,qb3,qb4,qb5,qb6,t_sk,t_sp)


#Calculate the twist of the aileron after calculating the shear flow due to torsion and the twist rate at each spanwise location

T_lst = [-2,-1,-3]   #change to actual torque distribution along the span for a number of points based on Willems code   
G = 28e9  #shear modulus of the aileron
def Twist_of_aileron(T_lst,per_semicircle, per_triangle,t_sk,t_sp,h,l_a,A1,A2):    #Solves a set of 3 equations for unit torque applied, output q1, q2 and twsit_rate_times_G
    #For T_lst use the value of the torque at each location, for examply for adding a loop or using input out of a list of torques
    #For verification, one can change the perimiter of the circle and triangle by chaning the geometry and calculate the q1,q2 and twist rate and see wheter it makes sense or not
    #After computing twist rate, take distance from hinge line to shear center to compute the deflection of the hinge line
    k = 1/(2*A1) * (per_semicircle/t_sk + 2*h/t_sp) 
    l = 1/(2*A1)*-2*h/t_sp
    m = 1/(2*A2)*-2*h/t_sp
    n = 1/(2*A2) * (per_triangle/t_sk + 2*h/t_sp)
    B = np.array([[2*A1,2*A2,0],[k,l,-1],[m,n,-1]])
    twist_rate_lst = []
    q1_lst = []
    q2_lst = []
    x_theta_0 = l_a/2  #due to assumption around x, x_sc in middle [m]
    theta_0 = 0       #this is a boundary condition [rad]
    for i in range(len(T_lst)):
        w = np.array([T_lst[i],0,0])
        solution = np.linalg.solve(B,w)
        q1 = solution[0]
        q2 = solution[1]
        q1_lst.append(q1)
        q2_lst.append(q2)
        twist_rate_times_G = solution[2]
        twist_rate = twist_rate_times_G / G 
        twist_rate_lst.append(twist_rate)
    J = T_lst[-1] / twist_rate_times_G              #calculate the J for a combination of torque and twist rate
    dx = l_a/(len(T_lst)-1)         #step in x direction between the points where the torque is computed and thus where twist_rate is known
    n_steps = math.floor(x_theta_0/dx)      #number of full steps untill location of boundary condition reached, returns an integer
    twist_before_bc = sum([twist_rate_lst[j] for j in range(n_steps)]) * dx + theta_0       #twist of first section
    twist_lst = [twist_before_bc]
    twist_after_bc = 0
    for i in range(1,len(T_lst)):
        if i < n_steps:
            twist_before_bc = twist_before_bc - twist_rate_lst[i-1]*dx       #compute the twist of each section between two points (positive for positive twist rate)
            twist_lst.append(twist_before_bc)
            
        if i == n_steps:                                  #this is the section where the boundary condition is applied
            twist_lst.append(theta_0)                     #now the section of the boundary condition is reached, this entire section attains this value (neglecting the twist along the even smaller subsection if point of boundary condition falls in between two points)
        if i > n_steps:
            twist_after_bc = twist_after_bc + twist_rate_lst[i]*dx     #or -, plot if torque distribution is known. At the boundary condition, the sign of the twist should change 
            twist_lst.append(twist_after_bc)
    return q1_lst,q2_lst, J, twist_rate_lst, twist_lst         #J, twist rate and twist at every x location taken 
    
q1_lst, q2_lst, J, twist_rate_lst, twist_lst = Twist_of_aileron(T_lst,per_semicircle, per_triangle,t_sk,t_sp,h,l_a,A1,A2)



def Shear_stress_due_to_torsion(q1,q2,t_sk,t_sp):   
    tau_skin_cell_1_lst = [q1[i]/t_sk for i in range(len(q1))]
    tau_skin_cell_2_lst = [q2[i]/t_sk for i in range(len(q1))]
    tau_spar_lst = [(q2[i]-q1[i])/t_sp for i in range(len(q1))]
        
    return tau_skin_cell_1_lst, tau_skin_cell_2_lst, tau_spar_lst

tau_skin_cell_1_lst, tau_skin_cell_2_lst, tau_spar_lst = Shear_stress_due_to_torsion(q1_lst,q2_lst,t_sk,t_sp)

def Total_shear_stress(tau_1,tau_2,tau_3,tau_4,tau_5,tau_6, tau_skin_cell_1_lst,tau_skin_cell_2_lst,tau_spar_lst):
    """Calculates the total shear stress distribution tau_yz per section for every x location along the span and puts this in a list, so 6 lists for each spanwise location
       it returns a list of lists with every list a shear stress distribution at a specific x location.""" 
    total_shear_stress_distribution_at_every_x_loc = []
    for j in range(len(tau_skin_cell_1_lst)):
        tau_total_1_at_x_loc = [tau_1[i] + tau_skin_cell_1_lst[j] for i in range(len(tau_1))]
        tau_total_2_at_x_loc = [tau_2[i] + tau_spar_lst[j] for i in range(len(tau_2))]
        tau_total_3_at_x_loc = [tau_3[i] + tau_skin_cell_2_lst[j] for i in range(len(tau_3))]
        tau_total_4_at_x_loc = [tau_4[i] + tau_skin_cell_2_lst[j] for i in range(len(tau_4))]
        tau_total_5_at_x_loc = [tau_5[i] + tau_spar_lst[j] for i in range(len(tau_5))]
        tau_total_6_at_x_loc = [tau_6[i] + tau_skin_cell_1_lst[j] for i in range(len(tau_6))]
        total_shear_stress_distribution_at_every_x_loc.append(tau_total_1_at_x_loc + tau_total_2_at_x_loc + tau_total_3_at_x_loc + tau_total_4_at_x_loc + tau_total_5_at_x_loc + tau_total_6_at_x_loc)
    return total_shear_stress_distribution_at_every_x_loc  #example output[0] = shear stress distribution at first x location along the span, in the order of the sections 1,2,3,4,5,6

total_shear_stress_distribution_at_every_x_loc = Total_shear_stress(tau_1,tau_2,tau_3,tau_4,tau_5,tau_6, tau_skin_cell_1_lst,tau_skin_cell_2_lst,tau_spar_lst)


My = 1          #to be changed and computed by a function in Willems code, just a single value for a spefic x location
Mz = 1          #to be changed computed by a function in Willems code, just a single value for a spefic x location

#def Direct_stress_distribution(Mz,My,Iyy,Izz,zc,yco1,zco1,yco2,zco2,yco3,zco3,yco4,zco4,yco5,zco5,yco6,zco6):     #for a unit moment in x and y direction
#    """Computes the Direct stress distrubtion along the cross-section at each middle point of a section where the shear flow is calculated based on the Mx and My of a specific location along the span"""
#    sigma_xx_1 = [My*((zco1[i]+zco1[i+1])/2 - zc)/Iyy + Mz*(yco1[i]+yco1[i+1])/2 for i in range(len(zco1)-1)]            #zc is zlocation of centroid, compute y- and z-location of the middle of each section where the shear flow is computed
#    sigma_xx_2 = [My*((zco2[i]+zco2[i+1])/2 - zc)/Iyy + Mz*(yco2[i]+yco2[i+1])/2 for i in range(len(zco2)-1)]
#    sigma_xx_3 = [My*((zco3[i]+zco3[i+1])/2 - zc)/Iyy + Mz*(yco3[i]+yco3[i+1])/2 for i in range(len(zco3)-1)]
#    sigma_xx_4 = [My*((zco4[i]+zco4[i+1])/2 - zc)/Iyy + Mz*(yco4[i]+yco4[i+1])/2 for i in range(len(zco4)-1)]
#    sigma_xx_5 = [My*((zco5[i]+zco5[i+1])/2 - zc)/Iyy + Mz*(yco5[i]+yco5[i+1])/2 for i in range(len(zco5)-1)]
#    sigma_xx_6 = [My*((zco6[i]+zco6[i+1])/2 - zc)/Iyy + Mz*(yco6[i]+yco6[i+1])/2 for i in range(len(zco6)-1)]
#    sigma_xx_1.extend(sigma_xx_2 + sigma_xx_3 + sigma_xx_4 + sigma_xx_5 + sigma_xx_6)       #combine all section to get distribution in one list
#    return sigma_xx_1
#direct_stress_distribution = Direct_stress_distribution(Mz,My,Iyy,Izz,zc,yco1,zco1,yco2,zco2,yco3,zco3,yco4,zco4,yco5,zco5,yco6,zco6)
def Direct_stress_distribution(Mz,My,Iyy,Izz,zc,yco1,zco1,yco2,zco2,yco3,zco3,yco4,zco4,yco5,zco5,yco6,zco6):     #for a unit moment in x and y direction
    """Computes the Direct stress distrubtion along the cross-section at each point where the shear flow is calculated based on the Mx and My of a specific location along the span"""
    sigma_xx_1 = [My*(zco1[i] - zc)/Iyy + Mz*yco1[i] for i in range(len(zco1))]            #zc is zlocation of centroid, compute y- and z-location of the middle of each section where the shear flow is computed
    sigma_xx_2 = [My*(zco2[i] - zc)/Iyy + Mz*yco2[i] for i in range(len(zco2))]
    sigma_xx_3 = [My*(zco3[i] - zc)/Iyy + Mz*yco3[i] for i in range(len(zco3))]
    sigma_xx_4 = [My*(zco4[i]- zc)/Iyy + Mz*yco4[i] for i in range(len(zco4))]
    sigma_xx_5 = [My*(zco5[i] - zc)/Iyy + Mz*yco5[i] for i in range(len(zco5))]
    sigma_xx_6 = [My*(zco6[i] - zc)/Iyy + Mz*yco6[i] for i in range(len(zco6))]
    sigma_xx_1.extend(sigma_xx_2 + sigma_xx_3 + sigma_xx_4 + sigma_xx_5 + sigma_xx_6)       #combine all section to get distribution in one list
    return list(sigma_xx_1)
direct_stress_distribution = Direct_stress_distribution(Mz,My,Iyy,Izz,zc,yco1,zco1,yco2,zco2,yco3,zco3,yco4,zco4,yco5,zco5,yco6,zco6)


#Calculation of the Von-Mises stress distribution
#assume tau xy tau xz are neglegible compared to the shear acting in the yz plane, same for sigma_yy and sigma_zz --> calculation of von mises stress simplifies
tau_xy = 0
tau_xz = 0
sigma_yy = 0
sigma_zz = 0


n = len(T_lst)  #This should be changed to the variable that defines the number of points taken where the loads are discretised along the span including begin and end points, in my program it matches the size of the torque list but that also changes according to the number of points along x chosen
def Von_Mises_stress_distribution(direct_stress_distribution,shear_stress_distribution,n):     #use total shear stresses!
    """This function calculates the Von Mises stress distrubtion for every x-location taken along the span, based on the direct stress calculation and the total shear stress calculation"""
    sigma_vm_distribution_at_every_x_loc = []
    for j in range(n):
        sigma_vm = [np.sqrt(direct_stress_distribution[i]**2 + 3*shear_stress_distribution[n-1][i]**2) for i in range(len(direct_stress_distribution))]
        sigma_vm_distribution_at_every_x_loc.append(sigma_vm)
    return sigma_vm_distribution_at_every_x_loc  #example sigma_vm_distribution[0] is the distribution at the first x-location taken along the span

Von_Mises_stress_distribution_at_every_x_loc = Von_Mises_stress_distribution(direct_stress_distribution,total_shear_stress_distribution_at_every_x_loc,n)

print()

#-------------------------------------------------------------------------------------------------------------------------------------------------

#for i,j in zip([yco2,yco3,yco4,yco5,yco6],[zco2,zco3,zco4,zco5,zco6]):
#    yco1.extend(i)
#    zco1.extend(j)
#  
#
#def plot_shearflow(zco,yco,Von_Mises_stress_distribution_at_every_x_loc):
#    n = 2           #distribution calculated at first point
#    plt.scatter(zco1,yco1,c=Von_Mises_stress_distribution_at_every_x_loc[n],cmap='jet')
#    plt.xlabel('z[m]', fontsize=16)
#    plt.ylabel('y[m]', fontsize=16)
#    #plt.legend()
#    c=plt.colorbar()
#    c.set_label('Von Mises Stress [N/m^2]')
#    plt.title('Von Mises Stress Distribution')
#    #plt.savefig("flow.png")
#    plt.show()
#plot_shearflow(zco1,yco1,Von_Mises_stress_distribution_at_every_x_loc)





