def fluidProperties(fluidName,EoS,P_inlet,T_inlet,P_outlet):
    import CoolProp
    import sys
    fluid = CoolProp.AbstractState(EoS,fluidName)
    eta_c = 0.8
    '''
    P_inlet = 2e5
    T_inlet = 120.5+273.15#119.62+273.15
    P_outlet = 6e5
    '''
    P_cyl_intake = 1.8e5
    P_cyl_discharge = 6.2e5

    # Determine temperature conditions from isentropic efficiency
    fluid.update(CoolProp.PQ_INPUTS,P_inlet,1)
    SH_inlet = T_inlet-fluid.T()
    if SH_inlet < 0:
        print('Specified temperature gives two-phase conditions at inlet (SH = ',SH_inlet,'). Exiting calculation.')
        sys.exit()
    fluid.update(CoolProp.PT_INPUTS,P_inlet,T_inlet)
    rho_inlet = fluid.rhomass()
    s_inlet = fluid.smass()
    h_inlet = fluid.hmass()
    fluid.update(CoolProp.PSmass_INPUTS,P_outlet,s_inlet)
    T_is_outlet = fluid.T()
    h_is_outlet = fluid.hmass()
    h_outlet = h_inlet + (h_is_outlet-h_inlet)/eta_c
    fluid.update(CoolProp.HmassP_INPUTS,h_outlet,P_outlet)
    rho_outlet = fluid.rhomass()
    T_outlet = fluid.T()
    s_outlet = fluid.smass()
    fluid.update(CoolProp.PSmass_INPUTS,P_cyl_intake,s_inlet)
    T_cyl_intake = fluid.T()
    fluid.update(CoolProp.PSmass_INPUTS,P_cyl_discharge,s_outlet)
    T_cyl_discharge = fluid.T()
    return rho_inlet, rho_outlet
    '''
    print('Inlet: P = ',P_inlet,' (Pa); ','T = ',T_inlet, '(K)')
    print('Cylinder intake: P = ',P_cyl_intake,' (Pa); ','T = ',T_cyl_intake, '(K)')
    print('Outlet: P = ',P_outlet,' (Pa); ','T = ',T_outlet, '(K)')
    print('Cylinder discharge: P = ',P_cyl_discharge,' (Pa); ','T = ',T_cyl_discharge, '(K)')
    fluid.update(CoolProp.PQ_INPUTS,P,0)
    T = fluid.T() 
    print('Saturation temperature (K): ',T)
    '''
