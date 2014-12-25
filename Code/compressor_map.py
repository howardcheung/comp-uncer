from __future__ import division #Make integer 3/2 give 1.5 in python 2.x
from CoolProp.CoolProp import Props
import pandas as PD
from math import isnan

#Credit: This file is based on ACHP's compressor.py

class CompressorClass():
    """
    Compressor Model based on 10-coefficient Model from `ANSI/AHRI standard 540 <http://www.ahrinet.org/App_Content/ahri/files/standards%20pdfs/ANSI%20standards%20pdfs/ANSI-ARI-540-2004%20latest.pdf>`_

    
    Required Parameters:
        
    ===========   ==========  ========================================================================
    Variable      Units       Description
    ===========   ==========  ========================================================================
    M             varied      A numpy-like list of compressor map coefficients for mass flow
    P             varied      A numpy-like list of compressor map coefficients for electrical power
    Ref           N/A         A string representing the refrigerant
    Tin_r         K           Refrigerant inlet temperature
    pin_r         kPa         Refrigerant suction pressure (absolute)
    pout_r        kPa         Refrigerant discharge pressure (absolute)
    fp            --          Fraction of electrical power lost as heat to ambient
    Vdot_ratio    --          Displacement Scale factor
    ===========   ==========  ========================================================================
    
    All variables are of double-type unless otherwise specified
        
    """
    def __init__(self,**kwargs):
        #Load up the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
    def Update(self,**kwargs):
        #Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
    def get_set_coeffs(self, **kwargs):
        filename = kwargs['filename']
        df=PD.read_csv(filename,header=0)
        m_coeffs = df['Mass Flow'].values
        pwr_coeffs = df['Power'].values
        kwargs={
              'M': m_coeffs,
              'P': pwr_coeffs
              }
        self.__dict__.update(kwargs)
        
    def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items with indices:
                [0] Description of value
                
                [1] Units of value
                
                [2] The value itself
        """

        if hasattr(self,'fp'):
            fp = self.fp
        else:
            fp = 999

        if hasattr(self,'Vdot_ratio'):
            Vdot_ratio = self.Vdot_ratio
        else:
            Vdot_ratio = 1  
        
        return [
            ('M1','-',self.M[0]),
            ('M2','-',self.M[1]),
            ('M3','-',self.M[2]),
            ('M4','-',self.M[3]),
            ('M5','-',self.M[4]),
            ('M6','-',self.M[5]),
            ('M7','-',self.M[6]),
            ('M8','-',self.M[7]),
            ('M9','-',self.M[8]),
            ('M10','-',self.M[9]),
            ('P1','-',self.P[0]),
            ('P2','-',self.P[1]),
            ('P3','-',self.P[2]),
            ('P4','-',self.P[3]),
            ('P5','-',self.P[4]),
            ('P6','-',self.P[5]),
            ('P7','-',self.P[6]),
            ('P8','-',self.P[7]),
            ('P9','-',self.P[8]),
            ('P10','-',self.P[9]),
            ('Heat Loss Fraction','-',fp),
            ('Displacement scale factor','-',Vdot_ratio),
            ('Power','W',self.W),
            ('Mass flow rate','kg/s',self.mdot_r),
            ('Inlet Temperature','K',self.Tin_r),
            ('Inlet Pressure','kPa',self.pin_r),
            ('Outlet Temperature','K',self.Tout_r),
            ('Outlet Pressure','kPa',self.pout_r),
            ('Inlet Enthalpy','J/kg',self.hin_r),
            ('Outlet Enthalpy','J/kg',self.hout_r),
            ('Overall isentropic efficiency','-',self.eta_oi),
            ('Pumped flow rate','m^3/s',self.Vdot_pumped)
         ]
        
    def Calculate(self):
        #Local copies of coefficients
        P=self.P
        M=self.M
        
        #Calculate suction superheat and dew temperatures
        self.Tsat_s_K=Props('T','P',self.pin_r,'Q',1.0,self.Ref)
        self.Tsat_d_K=Props('T','P',self.pout_r,'Q',1.0,self.Ref)
        self.DT_sh_K=self.Tin_r-self.Tsat_s_K
        
        #Convert saturation temperatures in K to F
        Tsat_s = self.Tsat_s_K * 9./5. - 459.67
        Tsat_d = self.Tsat_d_K * 9./5. - 459.67
    
        #Apply the 10 coefficient ARI map to saturation temps in F
        power_map = P[0] + P[1] * Tsat_s + P[2] * Tsat_d + P[3] * Tsat_s**2 + P[4] * Tsat_s * Tsat_d + P[5] * Tsat_d**2 + P[6] * Tsat_s**3 + P[7] * Tsat_d * Tsat_s**2 + P[8] * Tsat_d**2*Tsat_s + P[9] * Tsat_d**3
        mdot_map = M[0] + M[1] * Tsat_s + M[2] * Tsat_d + M[3] * Tsat_s**2 + M[4] * Tsat_s * Tsat_d + M[5] * Tsat_d**2 + M[6] * Tsat_s**3 + M[7] * Tsat_d * Tsat_s**2 + M[8] * Tsat_d**2*Tsat_s + M[9] * Tsat_d**3

        self.mdot_map_lbm_h = mdot_map
    
        # Convert mass flow rate to kg/s from lbm/h
        mdot_map *= 0.000125998

        self.power_map = power_map
        self.mdot_map = mdot_map
    
        # Add more mass flow rate to scale
        if hasattr(self,'Vdot_ratio'):
            #this tuning assumes linearity between mass flowrate and displacement
            mdot_map*=self.Vdot_ratio
            power_map*=self.Vdot_ratio

        P1 = self.pin_r
        P2 = self.pout_r
        T1_actual = self.Tsat_s_K + self.DT_sh_K
        try:
            #print "density call in compressor with", self.Tsat_s_K + 20.0/9.0*5.0,P1
            v_map = 1 / Props('D', 'T', self.Tsat_s_K + 20.0/9.0*5.0, 'P', P1, self.Ref)
        except: 
            print "P1, self.Ref,self.Tsat_s_K",P1, self.Ref,self.Tsat_s_K
            raise
        v_actual = 1 / Props('D', 'T', self.Tsat_s_K + self.DT_sh_K, 'P', P1, self.Ref)
        F = 0.75
        mdot = (1 + F * (v_map / v_actual - 1)) * mdot_map
    
        T1_map = self.Tsat_s_K + 20 * 5 / 9
        s1_map = Props('S', 'T', T1_map, 'P', P1, self.Ref)
        h1_map = Props('H', 'T', T1_map, 'P', P1, self.Ref)
        h2s_map = Props('H', 'S', s1_map, 'P', P2, self.Ref)
    
        s1_actual = Props('S', 'T', T1_actual, 'P', P1, self.Ref)
        h1_actual = Props('H', 'T', T1_actual, 'P', P1, self.Ref)
        h2s_actual = Props('H', 'S', s1_actual, 'P', P2, self.Ref)
    
        #Shaft power based on 20F superheat calculation from fit overall isentropic efficiency
        power = power_map * (mdot / mdot_map) * (h2s_actual - h1_actual) / (h2s_map - h1_map)
    
        h2 = power/1000 * (1 - self.fp) / mdot + h1_actual#
        self.eta_oi=mdot*(h2s_actual-h1_actual)/(power/1000)
        self.Tout_r = Props('T', 'H', h2, 'P',P2 , self.Ref)  # T_hp(self.Ref, h2, P2, T1_map + 40) #Plus 20 for guess value for discharge temp
        if isnan(self.Tout_r ):
            print "Problem solving for exit temperature Tout_r in Compressor.py"
        self.sout_r = Props('S','T',self.Tout_r,'P',P2,self.Ref) * 1000
        self.sin_r = Props('S','T',self.Tin_r,'P',P1,self.Ref) * 1000
        self.hout_r = h2 * 1000
        self.hin_r = h1_actual * 1000
        self.mdot_r=mdot
        if hasattr(self,'Pwr_tune'):
            self.W = power*self.Pwr_tune #self.Pwr_tune is applied after displacement factor
        else:
            self.W = power #only corrected by Vdot_ratio, if set
        self.CycleEnergyIn=power*(1-self.fp) #energy that goes into the cycle, excluding heat losses from the compressor itself
        self.Vdot_pumped=mdot/Props('D','T',self.Tin_r,'P',P1,self.Ref)
        
if __name__=='__main__':
    def F2K(T_F):
        """Convert temperature in Fahrenheit to Kelvin, code from ACHP"""
        return 5./9.*(T_F+459.67)
    
    Comp=CompressorClass()
    kwds = {
            'filename':'../Data/mfg_data_sheets/compressor/H23A463DBL.csv',
            }
    Comp.get_set_coeffs(**kwds)
    ref = 'R22'
    kwds = {
            'Ref':ref,
            'Tin_r':F2K(-20)+20*5./9.,
            'pin_r':Props('P','T',F2K(-20),'Q',1.0,ref),
            'pout_r':Props('P','T',F2K(80),'Q',1.0,ref),
            'fp':0.15#Fractionofelectricalpowerlostasheattoambient
            #'Vdot_ratio':1.0,#DisplacementScalefactor
            #'Pwr_tune':1.0#powertuningfactor,appliedafterVdot_ratio
            }
    Comp.Update(**kwds)
    Comp.Calculate()
    print "Map data >> flow in lbm/hr: ", Comp.mdot_map_lbm_h, "power in W", Comp.power_map
    print "Corrected data >> flow in lbm/hr: ", Comp.mdot_r/0.000125998, "power in W", Comp.W, "discharge temp [C]", Comp.Tout_r-273.15
