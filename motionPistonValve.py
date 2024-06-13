import math
import numpy
from scipy.optimize import fsolve
import matplotlib.pyplot as plot
from functionsPistonValve import *
from findValveParameters import *

# All data are in meters 
# x = geometry coordinate
# z = position (variable with time)
# t = time
# _p = piston
# _v = valve
# _op1 = opening 1
# _op2 = opening 2

P_filename = 'p_alpha_spilling.csv'

# Geometry (taken from TDC)
x_p_0 = 0.09351
x_v_0 = 0.05901
x_op1 = 0.085
x_op2 = 0.243
H_p = 0.162
H_v = 0.23
L_p = 0.092
L_v = 0.0335
L_vb = 0.123
L_op1 = 0.02
L_op2 = 0.02
D_p = 0.1025*2
D_v = 0.04*2
D_vi = 0.0285*2
Alpha_upper_exp, Alpha_lower_exp, P_upper_exp, P_lower_exp = importPressure(P_filename)

# Operating conditions
f = 50 # frequency (1/s)
delta = math.pi*55/180
t_0 = 0
z_p_0 = x_p_0+(H_p-L_p)/2       # piston at mid
z_v_t = x_v_0 + (H_v-2*L_v-L_vb)
z_v_0 = 0.0954
delta = findValveParameters(f,t_0,H_v,L_v,L_vb,z_v_t,z_v_0)
delta_deg = math.degrees(delta)
# z_v_0 = position(delta/(2*math.pi*f),f,delta,t_0,z_p_0,x_v_0+(H_v-L_vb-2*L_v)/2,H_p,H_v,L_p,L_v,L_vb)[1][0]

# Find valve parameters
'''
T_input = numpy.array([0.0, 90.0, 180.0, 270.0])*math.pi/180/f
Z_v_0_input = numpy.array([95.4, 90.5, 63.5, 65.8])/1000
res0 = findValveParameters([x_v_0,H_v,delta],Z_v_0_input,T_input,f,t_0,z_p_0,H_p,L_p,L_v,L_vb)
valve = lambda var: findValveParameters(var,Z_v_0_input,T_input,f,t_0,z_p_0,H_p,L_p,L_v,L_vb)[0][1:4]
valve_param = fsolve(valve,[x_v_0,H_v,delta])
x_v_0, H_v, delta = valve_param
z_v_0 = findValveParameters(valve_param,Z_v_0_input,T_input,f,t_0,z_p_0,H_p,L_p,L_v,L_vb)[1]
#print(valve_param)
#print(findValveParameters(valve_param,Z_v_0_input,T_input,f,t_0,z_p_0,H_p,L_p,L_v,L_vb))
'''
# Calculations
N = 200
T = numpy.linspace(0,1/f,N)
Alpha = T*360*f
Z_p = numpy.zeros((N,2))
Z_v = numpy.zeros((N,4))
OD_1 = numpy.zeros((N,1))
OD_2 = numpy.zeros((N,1))
for i,t in enumerate(T):
    Z_p[i,:], Z_v[i,:], = position(t,f,delta,t_0,z_p_0,z_v_0,H_p,H_v,L_p,L_v,L_vb)
    OD_1[i], OD_2[i] = openings(Z_p[i,:],Z_v[i,:],x_op1,x_op2,L_op1,L_op2,x_p_0,H_p)

# Plot motion of piston and valve and opening degree
fig1, ax = plot.subplots(1,3,figsize=(10,4),layout='constrained')
ax[0].plot(Alpha,Z_p[:,0],color='blue')
ax[0].plot(Alpha,Z_p[:,1],color='blue')
ax[0].fill_between(Alpha,Z_p[:,0],Z_p[:,1],color='blue',alpha=0.5)
ax[0].plot([Alpha[0],Alpha[-1]],[x_p_0, x_p_0],color='black')
ax[0].plot([Alpha[0],Alpha[-1]],[x_p_0+H_p, x_p_0+H_p],color='black')
ax[0].set_xlabel('Crank angle (°)')
ax[0].set_ylabel('Piston position (m)')
ax[1].plot(Alpha,Z_v[:,0],color='red')
ax[1].plot(Alpha,Z_v[:,1],color='red')
ax[1].plot(Alpha,Z_v[:,2],color='red')
ax[1].plot(Alpha,Z_v[:,3],color='red')
ax[1].fill_between(Alpha,Z_v[:,0],Z_v[:,1],color='red',alpha=0.5)
ax[1].fill_between(Alpha,Z_v[:,2],Z_v[:,3],color='red',alpha=0.5)
ax[1].plot([Alpha[0],Alpha[-1]],[x_v_0, x_v_0],color='black')
ax[1].plot([Alpha[0],Alpha[-1]],[x_v_0+H_v, x_v_0+H_v],color='black')
ax[1].set_xlabel('Crank angle (°)')
ax[1].set_ylabel('Valve position (m)')
ax[2].plot(Alpha,OD_1,color='green',marker='o',markersize=2)
ax[2].plot(Alpha,OD_2,color='purple',marker='o',markersize=2)
ax[2].plot([Alpha[0],Alpha[-1]],[0,0],color='black')
ax[2].plot([Alpha[0],Alpha[-1]],[1,1],color='black')
ax[2].set_xlabel('Crank angle (°)')
ax[2].set_ylabel('Opening degree (-)')
ax[2].legend(('op_1','op_2'))
fig1.suptitle('Motion of piston and valve and opening degree')
plot.tight_layout()
# Plot instantaneous position piston and valve 
i = math.ceil(225/360*N)
fig2 = plot.figure()
plot.plot([0,D_p,D_p,0,0],[x_p_0,x_p_0,x_p_0+H_p,x_p_0+H_p,x_p_0],color='black')
plot.plot([D_p,D_p+D_v,D_p+D_v,D_p,D_p],[x_v_0,x_v_0,x_v_0+H_v,x_v_0+H_v,x_v_0],color='black')
plot.plot([D_p*(1-1/30),D_p*(1+1/30),D_p*(1+1/30),D_p*(1-1/30),D_p*(1-1/30)],\
          [x_op1,x_op1,x_op1+L_op1,x_op1+L_op1,x_op1],color='black')
plot.fill_between([D_p*(1-1/30),D_p*(1+1/30)],[x_op1,x_op1],[x_op1+L_op1,x_op1+L_op1],color='black')
plot.plot([D_p*(1-1/30),D_p*(1+1/30),D_p*(1+1/30),D_p*(1-1/30),D_p*(1-1/30)],\
          [x_op2,x_op2,x_op2+L_op2,x_op2+L_op2,x_op2],color='black')
plot.fill_between([D_p*(1-1/30),D_p*(1+1/30)],[x_op2,x_op2],[x_op2+L_op2,x_op2+L_op2],color='black')
plot.plot([0,D_p,D_p,0,0],[Z_p[i,0],Z_p[i,0],Z_p[i,1],Z_p[i,1],Z_p[i,0]],color='blue')
plot.fill_between([0,D_p],[Z_p[i,0],Z_p[i,0]],[Z_p[i,1],Z_p[i,1]],color='blue',alpha=0.5)
plot.plot([D_p,D_p+D_v,D_p+D_v,D_p+D_vi,D_p+D_vi,D_p+D_v,D_p+D_v,D_p,D_p,D_p+D_v-D_vi,D_p+D_v-D_vi,D_p,D_p],\
          [Z_v[i,0],Z_v[i,0],Z_v[i,1],Z_v[i,1],Z_v[i,2],Z_v[i,2],Z_v[i,3],Z_v[i,3],Z_v[i,2],Z_v[i,2],Z_v[i,1],Z_v[i,1],Z_v[i,0]],color='red')
plot.fill_betweenx([Z_v[i,0],Z_v[i,1],Z_v[i,1],Z_v[i,2],Z_v[i,2],Z_v[i,3]],\
                   [D_p,D_p,D_p+D_v-D_vi,D_p+D_v-D_vi,D_p,D_p],\
                    [D_p+D_v,D_p+D_v,D_p+D_vi,D_p+D_vi,D_p+D_v,D_p+D_v],color='red',alpha=0.5)
#plot.title('Position of piston and valve at t = %.4f ms' % T[i])
plot.title('Position of piston and valve at alpha = %.4f °' % Alpha[i])
ax2 = plot.gca()
ax2.set_xticks([])
ax2.set_yticks([])
print([Z_p[i,0],Z_v[i,0]])

# Plot experimental pressure as a function of crank angle
fig3 = plot.figure()
plot.plot(Alpha_upper_exp,P_upper_exp,color='red')
plot.plot(Alpha_lower_exp,P_lower_exp,color='blue')
plot.xlabel('Crank angle (°)')
plot.ylabel('Pressure (bar)')
plot.legend(('upper chamber','lower chamber'))
plot.show()
print('OK')