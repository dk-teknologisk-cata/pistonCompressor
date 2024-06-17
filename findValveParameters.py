'''
def findValveParameters(var,Z_v_0_input,T_input,f,t_0,z_p_0,H_p,L_p,L_v,L_vb):
    import numpy
    import math
    from functionsPistonValve import position
    x_v_0, H_v, delta = var
    z_v_0 = position(delta/(2*math.pi*f),f,delta,t_0,z_p_0,x_v_0+(H_v-L_vb-2*L_v)/2,H_p,H_v,L_p,L_v,L_vb)[1][0]
    Z_v_0 = numpy.zeros(4)
    for i,t in enumerate(T_input):
        Z_v_0[i] = position(t,f,delta,t_0,z_p_0,z_v_0,H_p,H_v,L_p,L_v,L_vb)[1][0]
    res = Z_v_0-Z_v_0_input
    return res, z_v_0
'''
def findValveParameters(H_v,L_v,L_vb,x_v_0,z_v_0):
    import math
    delta = -math.asin((x_v_0+(H_v-2*L_v-L_vb)/2-z_v_0)/((H_v-2*L_v-L_vb)/2))
    return delta 