import math
import numpy
from scipy.optimize import fsolve
import matplotlib.pyplot as plot
from functionsPistonValve import *
from findValveParameters import *
from fluidProperties import *

# All data are in meters 
# x = geometry coordinate
# z = position (variable with time)
# t = time
# _p = piston
# _v = valve
# _op1 = opening 1
# _op2 = opening 2

P_filename = 'p_alpha_spilling.csv'
calcValveDelay = True
calcTimestep = True
calcAverageVelocity = True

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
D_op1 = 0.055
L_op2 = 0.02
D_op2 = 0.055
D_p = 0.1025*2
D_v = 0.04*2
D_vi = 0.0285*2
Alpha_upper_exp, Alpha_lower_exp, P_upper_exp, P_lower_exp = importPressure(P_filename)

# Operating conditions
fluidName = 'Water'
EoS = 'HEOS'
P_in = 2e5          # intake pressure (Pa)
T_in = 120.5+273.15 # intake temperature (K)
P_dis = 6e5         # discharge pressure (Pa)
m = 250/3600       # mass flow rate (kg/s)
rotSpeed = 1500     # rotational speed (rpm)
maxTranslation = 0.1*L_op1
doubleEffect = True
delta_deg = 55
delta = math.radians(delta_deg)
t_0 = 0
z_v_0 = 0.0954  # initial position valve (to determine valve delay)
alpha_deg = 0    # position to show

# Determine valve delay
z_p_0 = x_p_0 + (H_p-L_p)/2       # piston at mid
if calcValveDelay:
    delta = findValveParameters(H_v,L_v,L_vb,x_v_0,z_v_0)
    delta_deg = math.degrees(delta)
    print('Calculated valve delay: ',delta_deg)

# Determine piston-valve motion
f = rotSpeed/60
N = 200
T = numpy.linspace(0,1/f,N)
Alpha = T*360*f
Z_p = numpy.zeros((N,2))
Z_v = numpy.zeros((N,4))
OD_1 = numpy.zeros((N,))
OD_2 = numpy.zeros((N,))
OD_cyl_1 = numpy.zeros((N,))
OD_cyl_2 = numpy.zeros((N,))
for i,t in enumerate(T):
    Z_p[i,:], Z_v[i,:], = position(t,f,delta,t_0,z_p_0,x_v_0,H_p,H_v,L_p,L_v,L_vb)
    OD_1[i], OD_2[i] = openings(Z_p[i,:],Z_v[i,:],x_op1,x_op2,L_op1,L_op2,x_p_0,H_p)[0:2]
    if calcAverageVelocity:
        OD_cyl_1[i], OD_cyl_2[i] = openings(Z_p[i,:],Z_v[i,:],x_op1,x_op2,L_op1,L_op2,x_p_0,H_p)[2:4]


# Determine max timestep for max position change
if calcTimestep:
    v_p_max, v_v_max = velocity(0,f,delta,t_0,H_p,H_v,L_p,L_v,L_vb)[4:6]
    dt_max = min(maxTranslation/v_p_max,maxTranslation/v_v_max)
    print('Max timestep (s): ', dt_max)

# Determine average velocities
if calcAverageVelocity:
    M_cycle = m/f/2
    A_eff_op1 = OD_cyl_1*(2*D_op1)*L_op1
    A_eff_op2 = OD_cyl_2*(2*D_op2)*L_op2
    rho_in, rho_dis, T_cyl_dis = fluidProperties(fluidName,EoS,P_in,T_in,P_dis)
    V_avg_in_op1 = M_cycle/(rho_in*sum(A_eff_op1)*(1/f)/N)
    V_avg_dis_op1 = M_cycle/(rho_dis*sum(A_eff_op1)*(1/f)/N)
    V_avg_in_op2 = M_cycle/(rho_in*sum(A_eff_op1)*(1/f)/N)
    V_avg_dis_op2 = M_cycle/(rho_dis*sum(A_eff_op1)*(1/f)/N)
    A_cyl_op1 = D_p*(Z_p[:,0]-x_p_0)
    A_cyl_op2 = D_p*(x_p_0+H_p-Z_p[:,1])
    V_cyl_in_op1 = V_avg_in_op1*A_eff_op1/A_cyl_op1
    V_cyl_in_op2 = V_avg_in_op2*A_eff_op2/A_cyl_op2
    V_cyl_dis_op1 = V_avg_dis_op1*A_eff_op1/A_cyl_op1
    V_cyl_dis_op2 = V_avg_dis_op2*A_eff_op2/A_cyl_op2
    fig0, ax0 = plot.subplots(1,2,figsize=(8,4),layout='constrained')
    ax0[0].plot(Alpha,V_cyl_in_op1,color='blue')
    ax0[0].plot(Alpha,V_cyl_dis_op1,color='red')
    ax0[0].set_ylim(0, max(V_cyl_in_op1[i],V_cyl_dis_op1[i],V_cyl_in_op2[i],V_cyl_dis_op2[i])*2)
    ax0[0].set_title('Upper chamber')
    ax0[0].set_xlabel('Crank angle (°)')
    ax0[0].set_ylabel('Velocity (m/s)')
    ax0[0].legend(('V_in','V_dis'))
    ax0[1].plot(Alpha,V_cyl_in_op2,color='blue')
    ax0[1].plot(Alpha,V_cyl_dis_op2,color='red')
    ax0[1].set_ylim(0, max(V_cyl_in_op2[i],V_cyl_dis_op2[i],V_cyl_in_op1[i],V_cyl_dis_op1[i])*2)
    ax0[1].set_title('Lower chamber')
    ax0[1].set_xlabel('Crank angle (°)')
    ax0[1].set_ylabel('Velocity (m/s)')
    ax0[1].legend(('V_in','V_dis'))
    i = math.ceil(alpha_deg/360*N)
    print('V_cyl_in_op1 = ', V_cyl_in_op1[i], 'at alpha_deg = ', alpha_deg)
    print('V_cyl_dis_op1 = ', V_cyl_dis_op1[i], 'at alpha_deg = ', alpha_deg)
    print('V_cyl_in_op2 = ', V_cyl_in_op2[i], 'at alpha_deg = ', alpha_deg)
    print('V_cyl_dis_op2 = ', V_cyl_dis_op2[i], 'at alpha_deg = ', alpha_deg)
    


# Plot motion of piston and valve and opening degree
fig1, ax1 = plot.subplots(1,3,figsize=(10,4),layout='constrained')
ax1[0].plot(Alpha,Z_p[:,0],color='blue')
ax1[0].plot(Alpha,Z_p[:,1],color='blue')
ax1[0].fill_between(Alpha,Z_p[:,0],Z_p[:,1],color='blue',alpha=0.5)
ax1[0].plot([Alpha[0],Alpha[-1]],[x_p_0, x_p_0],color='black')
ax1[0].plot([Alpha[0],Alpha[-1]],[x_p_0+H_p, x_p_0+H_p],color='black')
ax1[0].set_xlabel('Crank angle (°)')
ax1[0].set_ylabel('Piston position (m)')
ax1[1].plot(Alpha,Z_v[:,0],color='red')
ax1[1].plot(Alpha,Z_v[:,1],color='red')
ax1[1].plot(Alpha,Z_v[:,2],color='red')
ax1[1].plot(Alpha,Z_v[:,3],color='red')
ax1[1].fill_between(Alpha,Z_v[:,0],Z_v[:,1],color='red',alpha=0.5)
ax1[1].fill_between(Alpha,Z_v[:,2],Z_v[:,3],color='red',alpha=0.5)
ax1[1].plot([Alpha[0],Alpha[-1]],[x_v_0, x_v_0],color='black')
ax1[1].plot([Alpha[0],Alpha[-1]],[x_v_0+H_v, x_v_0+H_v],color='black')
ax1[1].set_xlabel('Crank angle (°)')
ax1[1].set_ylabel('Valve position (m)')
ax1[2].plot(Alpha,OD_1,color='green',marker='o',markersize=2)
ax1[2].plot(Alpha,OD_2,color='purple',marker='o',markersize=2)
ax1[2].plot([Alpha[0],Alpha[-1]],[0,0],color='black')
ax1[2].plot([Alpha[0],Alpha[-1]],[1,1],color='black')
ax1[2].set_xlabel('Crank angle (°)')
ax1[2].set_ylabel('Opening degree (-)')
ax1[2].legend(('op_1','op_2'))
fig1.suptitle('Motion of piston and valve and opening degree')
plot.tight_layout()

# Plot motion of piston and valve and opening degree together
fig2 = plot.figure()
plot.plot(Alpha,Z_p[:,0],color='blue')
plot.plot(Alpha,Z_p[:,1],color='blue')
plot.fill_between(Alpha,Z_p[:,0],Z_p[:,1],color='blue',alpha=0.5)
#plot.plot([Alpha[0],Alpha[-1]],[x_p_0, x_p_0],color='black')
#plot.plot([Alpha[0],Alpha[-1]],[x_p_0+H_p, x_p_0+H_p],color='black')
plot.plot(Alpha,Z_v[:,0],color='red')
plot.plot(Alpha,Z_v[:,1],color='red')
plot.plot(Alpha,Z_v[:,2],color='red')
plot.plot(Alpha,Z_v[:,3],color='red')
plot.fill_between(Alpha,Z_v[:,0],Z_v[:,1],color='red',alpha=0.5)
plot.fill_between(Alpha,Z_v[:,2],Z_v[:,3],color='red',alpha=0.5)
#plot.plot([Alpha[0],Alpha[-1]],[x_v_0, x_v_0],color='black')
#plot.plot([Alpha[0],Alpha[-1]],[x_v_0+H_v, x_v_0+H_v],color='black')
plot.plot([Alpha[0],Alpha[-1]],[x_op1, x_op1],color='black')
plot.plot([Alpha[0],Alpha[-1]],[x_op1+L_op1, x_op1+L_op1],color='black')
plot.fill_between([Alpha[0],Alpha[-1]],[x_op1, x_op1],[x_op1+L_op1, x_op1+L_op1],color='black',alpha=0.5)
plot.plot([Alpha[0],Alpha[-1]],[x_op2, x_op2],color='black')
plot.plot([Alpha[0],Alpha[-1]],[x_op2+L_op2, x_op2+L_op2],color='black')
plot.fill_between([Alpha[0],Alpha[-1]],[x_op2, x_op2],[x_op2+L_op2, x_op2+L_op2],color='black',alpha=0.5)
plot.xlabel('Crank angle (°)')
plot.ylabel('Axial coordinate (m)')
plot.title('Opening degree of cylinder - valve orifices')

# Plot instantaneous position piston and valve
i = math.ceil(alpha_deg/360*N)
fig3 = plot.figure()
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
plot.title('Position of piston and valve at alpha = %.2f °' % Alpha[i])
ax3 = plot.gca()
ax3.set_xticks([])
ax3.set_yticks([])
print([Z_p[i,0],Z_v[i,0]])

# Plot experimental pressure as a function of crank angle
fig4 = plot.figure()
plot.plot(Alpha_upper_exp,P_upper_exp,color='red')
plot.plot(Alpha_lower_exp,P_lower_exp,color='blue')
plot.xlabel('Crank angle (°)')
plot.ylabel('Pressure (bar)')
plot.legend(('upper chamber','lower chamber'))

# Plot opening degree and pressure in the cylinder as a function of the crank angle
fig5, ax5 = plot.subplots(1,2,figsize=(8,4),layout='constrained')
ax5[0].plot(Alpha,OD_1,color='green',marker='o',markersize=2)
ax5[0].set_xlabel('Crank angle (°)')
ax5[0].set_ylabel('Opening degree (-)',color='green')
ax5[0].tick_params(axis='y',labelcolor='green')
ax5_01 = ax5[0].twinx()
ax5_01.plot(Alpha_lower_exp,P_lower_exp,color='blue')
ax5_01.set_ylabel('Pressure (bar)',color='blue')
ax5_01.tick_params(axis='y',labelcolor='blue')
ax5[1].plot(Alpha,OD_2,color='purple',marker='o',markersize=2)
ax5[1].set_xlabel('Crank angle (°)')
ax5[1].set_ylabel('Opening degree (-)',color='purple')
ax5[1].tick_params(axis='y',labelcolor='purple')
ax5_11 = ax5[1].twinx()
ax5_11.plot(Alpha_upper_exp,P_upper_exp,color='red')
ax5_11.set_ylabel('Pressure (bar)',color='red')
ax5_11.tick_params(axis='y',labelcolor='red')
fig5.suptitle('Opening degree and pressure in the cylinder')
plot.tight_layout()

fig6 = plot.figure()
plot.plot(Alpha,OD_1*100,color='green')
plot.plot(Alpha,OD_2*100,color='purple')
plot.plot([Alpha[0],Alpha[-1]],[0,0],color='black')
plot.plot([Alpha[0],Alpha[-1]],[100,100],color='black')
plot.xlabel('Crank angle (°)')
plot.ylabel('Opening degree (%)')
plot.legend(('lower chamber','upper chamber'),loc=(0.6, 0.8))
plot.title('Opening degree of cylinder - valve orifices')

plot.show()
print('OK')

