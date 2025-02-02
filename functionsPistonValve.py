def position(t,f,delta,t_0,z_p_0,x_v_0,H_p,H_v,L_p,L_v,L_vb):
    import math
    z_p = (H_p-L_p)/2*(math.sin(2*math.pi*f*(t-t_0))) + z_p_0
    z_v = (H_v-2*L_v-L_vb)/2*(math.sin(2*math.pi*f*(t-t_0)-delta)) + x_v_0 + (H_v-2*L_v-L_vb)/2
    Z_p = [z_p, z_p + L_p]
    Z_v = [z_v, z_v + L_v, z_v + L_v + L_vb, z_v + 2*L_v + L_vb]
    return Z_p, Z_v

def position_orig(t,f,delta,t_0,z_p_0,x_v_0,H_p,H_v,L_p,L_v,L_vb):
    import math
    z_p = (H_p-L_p)/2*(math.sin(2*math.pi*f*(t-t_0))) + z_p_0
    z_v = (H_v-2*L_v-L_vb)/2*(math.sin(2*math.pi*(f+delta)*(t-t_0)-delta)) + x_v_0 + (H_v-2*L_v-L_vb)/2
    Z_p = [z_p, z_p + L_p]
    Z_v = [z_v, z_v + L_v, z_v + L_v + L_vb, z_v + 2*L_v + L_vb]
    return Z_p, Z_v

def velocity(t,f,delta,t_0,H_p,H_v,L_p,L_v,L_vb):
    import math
    v_p = (H_p-L_p)/2*2*math.pi*f*(math.cos(2*math.pi*f*(t-t_0)))
    v_v = (H_v-2*L_v-L_vb)/2*2*math.pi*(f+delta)*(math.cos(2*math.pi*(f+delta)*(t-t_0)-delta))
    t_max_v_p = t_0
    t_max_v_v = delta/(2*math.pi*(f+delta))+t_0
    v_p_max = (H_p-L_p)/2*2*math.pi*f
    v_v_max = (H_v-2*L_v-L_vb)/2*2*math.pi*(f+delta)
    return v_p, v_v, t_max_v_p, t_max_v_v, v_p_max, v_v_max

def openings(Z_p,Z_v,x_op1,x_op2,L_op1,L_op2,x_p_0,H_p):
    OD_1 = [0, 0]
    OD_2 = [0, 0]
    # Piston opening1
    if Z_p[0]<=max(x_p_0,x_op1) and Z_p[1]>=x_op1+L_op1:
        OD_1[0] = 0
    elif Z_p[1]<=x_op1 or Z_p[0]>=x_op1+L_op1:
        OD_1[0] = 1
    else:
        OD_1[0] = 0.5
    # Valve opening1
    if Z_v[0]<=x_op1 and Z_v[1]>=x_op1+L_op1:
        OD_1[1] = 0
    elif Z_v[1]<=x_op1 or Z_v[0]>=x_op1+L_op1:
        OD_1[1] = 1
    else:
        OD_1[1] = 0.5
    # Combined opening1
    od_1 = OD_1[0]*OD_1[1]
    if OD_1[1]<=0:
        od_cyl_1 = 0
    else: 
        od_cyl_1 = OD_1[0]
    if od_1>0 and od_1<1:
        ub = [x_op1+L_op1]
        if Z_p[0]>=max(x_op1,x_p_0):
            ub.append(Z_p[0])
        if Z_v[0]>=x_op1:
            ub.append(Z_v[0])
        lb = [x_op1] #[max(x_op1,x_p_0)]
        if Z_p[1]<x_op1+L_op1:
            lb.append(Z_p[1])
        if Z_v[1]<x_op1+L_op1:
            lb.append(Z_v[1])
        od_1 = max(0,(min(ub)-max(lb))/L_op1)
    if od_cyl_1>0 and od_cyl_1<1:
        ub_cyl = [x_op1+L_op1]
        if Z_p[0]>=x_op1:
            ub_cyl.append(Z_p[0])
        lb_cyl = [x_op1]
        if Z_p[1]<=min(x_p_0+H_p,x_op1+L_op1):
            lb_cyl.append(Z_p[1])
        od_cyl_1 = max(0,(min(ub_cyl[0:2])-max(lb_cyl[0:2]))/L_op1)

    # Piston opening2
    if Z_p[0]<=x_op2 and Z_p[1]>=min(x_p_0+H_p,x_op2+L_op2):
        OD_2[0] = 0
    elif Z_p[1]<=x_op2 or Z_p[0]>=x_op2+L_op2:
        OD_2[0] = 1
    else:
        OD_2[0] = 0.5
    # Valve opening2      
    if Z_v[2]<=x_op2 and Z_v[3]>=x_op2+L_op2:
        OD_2[1] = 0
    elif Z_v[3]<=x_op2 or Z_v[2]>=x_op2+L_op2:
        OD_2[1] = 1
    else:
        OD_2[1] = 0.5
    # Combined opening2
    od_2 = OD_2[0]*OD_2[1]
    if OD_2[1]<=0:
        od_cyl_2 = 0
    else:
        od_cyl_2 = OD_2[0]
    if od_2>0 and od_2<1:
        ub = [x_op2+L_op2] #[min(x_p_0+H_p,x_op2+L_op2)]
        if Z_p[0]>=x_op2:
            ub.append(Z_p[0])
        if Z_v[2]>=x_op2:
            ub.append(Z_v[2])
        lb = [x_op2]
        if Z_p[1]<=min(x_p_0+H_p,x_op2+L_op2):
            lb.append(Z_p[1])
        if Z_v[3]<=x_op2+L_op2:
            lb.append(Z_v[3])
        od_2 = max(0,(min(ub)-max(lb))/L_op2)
    if od_cyl_2>0 and od_cyl_2<1:
        ub_cyl = [x_op2+L_op2] #[min(x_p_0+H_p,x_op2+L_op2)]
        if Z_p[0]>=x_op2:
            ub_cyl.append(Z_p[0])
        lb_cyl = [x_op2]
        if Z_p[1]<=min(x_p_0+H_p,x_op2+L_op2):
            lb_cyl.append(Z_p[1])
        od_cyl_2 = max(0,(min(ub_cyl[0:2])-max(lb_cyl[0:2]))/L_op2)
    return od_1, od_2, od_cyl_1, od_cyl_2

def importPressure(filename):
    import pandas as pd
    import numpy as np
    pressureData = pd.read_csv(filename)
    Alpha_upper = pressureData.iloc[:71,0].to_numpy()
    Alpha_lower = pressureData.iloc[72:,0].to_numpy()
    P_upper = pressureData.iloc[:71,1].to_numpy()
    P_lower = pressureData.iloc[72:,1].to_numpy()
    return Alpha_upper, Alpha_lower, P_upper, P_lower