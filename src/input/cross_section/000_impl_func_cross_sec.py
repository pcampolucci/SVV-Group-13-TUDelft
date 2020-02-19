# input

import numpy as np

C_a = 0.547              #chord length aileron [m]
l_a = 2.771              #span of the aileron [m]
x1 = 0.153               #x-location of hing 1 [m]
x1 = 1.281               #x-location of hing 2 [m]
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



def Stiffener_spacing(h,C_a):
    perimeter_semicircle = np.pi*h
    perimeter_triangle = 2*np.sqrt((C_a-h)*2+h*2)
    stiffener_spacing_semicircle = perimeter_semicircle/6
    stiffener_spacing_triangle = perimeter_triangle/14
    return stiffener_spacing_triangle, stiffener_spacing_semicircle, perimeter_semicircle, perimeter_triangle        #[m]


#calculation of boom locations
def Boom_locations(h,C_a,stiffener_spacing_semicircle):
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

delta_st_triangle, delta_st_semicircle,per_semicircle,per_triangle = Stiffener_spacing(h,C_a,n_st)          # [m]
z_lst, y_lst = Boom_locations(h,C_a,delta_st_semicircle)     # [m]

#stiffener area calculation
def Stiffener_A(w_st,h_st,t_st):
    A_st = w_st*t_st+(h_st-t_st)*t_st
    return A_st
A_st=Stiffener_A(w_st,h_st,t_st)

#centroid
def Centroid(z_lst,A_st,C_a,t_sk,h):
    yc=0
    stff = sum([i * A_st for i in z_lst])
    triangle=-2*(h+(C_a-h)/2)*np.sqrt((C_a-h)*2+h*2)*t_sk
    semic=  -h*(1-2/np.pi)*np.pi(h*2-(h-t_sk)*2)/2
    spar=-t_sk*h*2*h

    zc=(stff+semic+triangle+spar)/(17*A_st+np.sqrt((C_a-h)*2+h*2)*t_sk*2+np.pi(h*2-(h-t_sk)*2)/2+t_sk*2*h)
    return (yc , zc)

yc,zc=Centroid(z_lst,A_st,C_a,t_sk,h)
#print(yc,zc)

#moment of inertia calculation
def MMI(z_lst,y_lst,A_st,C_a,t_sk,h,zc,yc):
    inclined_section = np.sqrt((C_a - h) * 2 + h * 2)

    stff_z= sum([(i - yc) ** 2 * A_st for i in y_lst])
    triangle_z = ((2*inclined_section)*3*t_sk/12)((h)/inclined_section)*2 + 2 * inclined_section * t_sk * (h/2 ) * 2
    spar_z = (2*h/12)(t_sk)*3
    semic_z = np.pi*h**3*t_sk/2
    Izz= stff_z+triangle_z+spar_z+semic_z

    stff_y = sum([(i-zc)**2*A_st for i in z_lst])
    triangle_y = ((2*inclined_section) * 3 * t_sk / 12) * ((C_a-h) / inclined_section) * 2 + 2*inclined_section*t_sk*(-(h+(C_a-h)/2)-zc)**2
    spar_y = (t_sk/12)*(2*h)*3 + 2*h*t_sk(-h-zc)**2
    semic_y = np.pi*h*3*t_sk/2 + (-h(1-2/np.pi)-zc)*2*np.pi(h*2-(h-t_sk)*2)/2
    Iyy= stff_y+triangle_y+spar_y+semic_y

    Ixx = 1/12*l_a*C_a**3

    return Izz, Ixx, Iyy

Izz,Ixx,Iyy = MMI(z_lst,y_lst,A_st,C_a,t_sk,h,zc,yc)