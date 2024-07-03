def findValveParameters(H_v,L_v,L_vb,x_v_0,z_v_0):
    import math
    delta = -math.asin((x_v_0+(H_v-2*L_v-L_vb)/2-z_v_0)/((H_v-2*L_v-L_vb)/2))
    return delta 