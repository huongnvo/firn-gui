'''
Firn Air Model
'''

import Tkinter
import tkFileDialog
from Tkinter import Menu
from Tkinter import *
import os
import csv
import json
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.sparse.linalg as splin
from scipy.sparse import spdiags
from scipy.sparse.linalg import lsqr
from scipy.integrate import cumtrapz
import math
import ModelParameters.Gasses as MPG
import ModelParameters.Sites as MPS
import ModelParameters.Plotting as plots
import ModelParameters.Diffusivity as MPD
import ModelParameters.density as MPRHO
import time
import airModel

class firnAirModel:
    def __init__(self, cc, dataDir, dataSave, text): 
        # load globals
        self.g        = 9.81                # m/s^2
        self.rho_i    = 917.0               # kg/m^3
        self.R        = 8.314               # J/mol/K
        self.M_air    = 28.97e-3            # kg/mol
        self.rho_bco  = 815.                # kg/m^3
        self.p_0      = 1.01325e5           # Standard Amtmospheric Pressure, Pa
        self.T_0      = 273.15              # Standard Temp, K
        self.sPerYear = 365.25 * 24 * 3600  # seconds per year
    
        self.cc       = cc                  # .json file
        self.dataDir  = dataDir             # data directory 
        self.dataSave = dataSave            # results directory
        self.text     = text
        
        reload(MPG)
        reload(MPS)
        reload(MPD)
        reload(plots)
        reload(MPRHO)
        
        tic = time.time()
            
        self.d = self.firnair()
        
        nodes = self.d['nodes']
        d15N2 = self.d['d15N2']
        d15   = d15N2[:,-1] - 1 #
        slope = (d15[100] - d15[50]) / (nodes[100] - nodes[50])
    
        # print plots
        #fig1 = plt.figure(1)
        #plt.clf()
        #plt.plot(nodes, d15N2[:,-1])
        #plt.show()
                   
        # print time     
        elapsed = time.time() - tic
        elapsed_min = elapsed / 60.
        mins = np.floor(elapsed_min)
        secs = (elapsed_min-mins) * 60
        string = mins, 'min', secs, 'sec elapsed'
        self.write(string)
    
    
    def w(self, por, zz, Accu, rho_interface, T, p_a):
        '''
        Function for downward advection of air and total air content

        Arguments:
        por -- 
        z_edges -- 
        Accu -- 
        rho_interface -- 
        T -- 
        p_a -- 
        z_nodes -- 
        dz -- 

        Returns:
        w_ad --
        bubble_pres -- 
        '''
        z_edges = zz['z_edges']
        z_nodes = zz['z_nodes']
        dz = zz['dz']
        
        por_tot_interface = np.interp(z_edges, z_nodes, por['por_tot'])
        por_cl_interface  = np.interp(z_edges, z_nodes, por['por_cl'])
        por_op_interface  = np.interp(z_edges, z_nodes, por['por_op'])
        teller_co = np.argmax(por_cl_interface)
        w_ice = Accu * self.rho_i / rho_interface 
        
        if self.cc['ad_method'] == 'ice_vel':
            w_ad    = w_ice
            trapped = 0.0
            bubble_pres = np.zeros_like(z_edges)
        
        elif self.cc['ad_method'] == 'Christo':
        ### Christo's Method from his thesis (chapter 5). This (maybe) could be vectorized to speed it up.
        
            bubble_pres = np.zeros_like(z_edges)
            dscl        = np.append(0, np.diff(por['por_cl']) / dz)    
            C           = np.exp(self.M_air * self.g * z_edges / (self.R * T))
            strain      = np.gradient(np.log(w_ice), dz)
            s           = por_op_interface + por_cl_interface
            
            for teller1 in range (0, teller_co + 1): 
                integral  = np.zeros(teller1 + 1)
                integral2 = np.zeros(teller1 + 1)
                
                for teller2 in range(0, teller1 + 1):
                    integral[teller2] = dscl[teller2] * C[teller2] * (s[teller2] / s[teller1]) / (1 + np.trapz(strain[teller2:teller1 + 1], dx = dz)) 
                    if dscl[teller2] == 0:
                        dscl[teller2] = 1e-14
                    integral2[teller2] = dscl[teller2]
                    
                bubble_pres[teller1] = (dz * np.sum(integral)) / (dz * np.sum(integral2))
            
            bubble_pres[teller_co + 1:] = bubble_pres[teller_co] * (s[teller_co] / s[teller_co + 1:]) / (w_ice[teller_co + 1:] / w_ice[teller_co])
            
            bubble_pres[0] = 1
            #print 'bubble pressure = %s' % bubble_pres
            
            flux     = w_ice[teller_co] * bubble_pres[teller_co] * por['por_cl'][teller_co]
            velocity = np.minimum(w_ice ,((flux + (1e-10) - w_ice * bubble_pres * por_cl_interface) / ((por_op_interface + 1e-10) *C)))        
            w_ad     = velocity

        return w_ad, bubble_pres
    
    def A(self, P): 
        '''
        Power-law scheme (from Pakantar, eq. 5.34)

        Arguments:
        P -- 

        Returns:
        A -- 
        '''
        A = np.maximum((1 - 0.1 * np.abs(P)) ** 5, np.zeros(np.size(P)))
        return A    
        
    def F_upwind(self, F): 
        '''
        Upwinding scheme

        Arguments:
        F --

        Returns:
        F_upwind -- 
        '''
        F_upwind = np.maximum(F, 0)
        return F_upwind
    
    def solver(self, a, b): #routine to solve Ax=b
        '''
        Routine to solve Ax = B

        Arguments:
        a -- 
        b -- 

        Returns:
        phi_t --
        '''
        a_U = a['a_U']
        a_D = a['a_D']
        a_P = a['a_P']

        nz    = np.size(b)
        Diags = (np.append([a_U, -a_P], [a_D], axis=0))
        cols  = np.array([1, 0, -1])      
        big_A = spdiags(Diags, cols, nz, nz, format = 'csc')
        big_A = big_A.T        
        rhs   = -b    
        phi_t = splin.spsolve(big_A, rhs)    
        return phi_t
    
    def FirnAir_SS(self, gaschoice): 
        '''
        Steady state analysis

        Arguments: 
        gaschoice

        Returns:
        phi --
        a_P_out --
        bubble_pres --
        z_nodes -- 
        ''' 
        loadgas     = True 
        depth       = self.cc['depth'] # m
        ad_method   = self.cc['ad_method'] # advection method 
        sitechoice  = self.cc['sitechoice'] # Set up parameters for different sites.
        z_res       = self.cc['z_res'] #resolution of grid, m
        stpsperyear = self.cc['steps'] #If this is for transient this number must (for now) be the same time steps as the input density/depth files. Make sure that it is a float.
        
        p_a, T, Accu_0, czd, z_co, LIZ, rho0, hemisphere = MPS.sites(sitechoice) 
        D_x, M, deltaM, conc1, d_0, omega = MPG.gasses(gaschoice, sitechoice, T, p_a, self.dataDir, hemisphere, loadgas)
        zz, nz_P, nz_fv = self.space(depth, z_res) #call on space function to set up spatial grid
        rhoHL = MPRHO.rhoHLAnalytic(self.R, T, self.rho_i, rho0, self.rho_bco, zz['z_nodes'], Accu_0) 
        rho_co, por, bcoRho, LIDRho = self.porosity(rhoHL, T)
        diffu, d_eddy = self.diffusivity(rho_co, por, z_co, czd, LIZ, d_0, D_x, T, p_a, zz['z_nodes'], rhoHL) #get diffusivity profiles
        
        Accu_m = Accu_0 # Accumulation in m/year
        Accu_0 = Accu_0 / self.sPerYear # accumulation in m/second
        
        time_yr   = conc1[:,0] # Atmospheric measurements times
        time_yr_s = time_yr * self.sPerYear  
        gas_org   = conc1[:,1] # Atmospheric measurement concentrations   
              
        nt, model_time, model_time_years, dt = self.time(time_yr, stpsperyear)

        rho_prof  = rhoHL # check this line

        dcon  = 1.0 
        diffu = diffu * dcon   
        gas   = np.interp(model_time, time_yr_s, gas_org) #interpolate atmospheric gas history to model time.
        bc_u, bc_d, bc_u_0 = self.boundaries(gas_org) #set boundary and initial conditions: bc_u is atmospheric condition, bc_d is zero gradient.
        
        phi_0 = np.zeros(nz_P) #phi is mixing ratio of gas.
        phi_0[:] = bc_u_0
        
        ConcPath  = os.path.join(self.dataSave, 'concentration.csv') #Save location.
        DepthPath = os.path.join(self.dataSave, 'depth.csv') #Save location.
        
        a, a_P_0, b_0, dZ_d, S_P, dZ, S_C_0, w_edges, bubble_pres = self.calculate(por, zz, diffu, d_eddy, deltaM, T, omega, None, phi_0, rho_prof, Accu_0, p_a, dt)
    
        s = (nz_P, nt)
        phi = np.zeros(s)
        a_P_out = np.zeros(s)
        
        phi_t = phi_0
        
        self.writeToFile(ConcPath, "w", 0, phi_t) 
        self.writeToFile(DepthPath, "w", 0, zz['z_nodes'])
                        
        for i_time in range(0, nt): #factor things
            
            bc_u_0  = gas[i_time]
            bc_type = 1.
            bc_u    = np.array([bc_u_0, bc_type])
            
            b = b_0 + a_P_0 * phi_t       
            
            #Up boundary
            a['a_P'][0] = 1 
            a['a_U'][0] = 0
            a['a_D'][0] = 0
            b[0]   = bc_u[0]
            
            #Down boundary
            a['a_P'][-1] = 1 
            a['a_D'][-1] = 0
            a['a_U'][-1] = 1
            b[-1]   = -dZ_d[-1] * bc_d[0] # probably does not matter as long as it is zero flux.
            
            phi_t = self.solver(a, b)
            
            phi[:,i_time] = phi_t
            a_P_out[:,i_time] = a['a_P']
                  
            a['a_P'] = a['a_U'] + a['a_D'] + a_P_0 - S_P * dZ
            
            S_C = S_C_0 * phi_t
            b_0 = S_C * dZ

            self.writeToFile(ConcPath, "a", model_time_years[i_time], phi_t)
            
        return phi, a_P_out, bubble_pres, zz['z_nodes']
         
    def FirnAir_TR(self, gaschoice, jj):
        '''
        Transient analysis

        Arguments:
        gaschoice --
        jj --

        Returns:
        phi --
        diffu_hold -- 
        rho_hold -- 
        z_nodes --
        '''
        ConcPath  = os.path.join(self.dataSave, 'conc_out_%s.csv') #Save location.
        RhoPath   = os.path.join(self.dataSave, 'rho_out.csv') #Save location.
        DiffuPath = os.path.join(self.dataSave, 'diffu_out.csv') #Save location.
        ZPath     = os.path.join(self.dataSave, 'Znodes_out.csv') #Save location.

        depth       = self.cc["depth"] # m
        ad_method   = self.cc["ad_method"] #advection method
        p_a         = self.cc['pressure'] #pressure at site. Could be vectorized to change through time.
        czd         = self.cc['ConZoneDepth'] #could also be a vector
        z_co_vec    = f_BCO[2,:] #2 for mart, 4 for 815
        LIZ_vec     = f_LID[2,:]
        rho0        = self.cc['rho0']
        hemisphere  = 'SCENARIO'
        sitechoice  = self.cc['sitechoice']
        z_res       = self.cc['z_resolution']
        stpsperyear = self.cc['StepsPerYear'] #If this is for transient this number must (for now) be the same time steps as the input density/depth files. Make sure that it is a float.

        zz, nz_P, nz_fv = self.space(depth, z_res)
        #D_x, M, deltaM, conc1, d_0, omega = MPG.gasses(gaschoice, sitechoice, T, p_a, self.dataDir, hemisphere, loadgas) #loadgas is False for userdata and true for not
        #gas = np.interp(model_time, time_yr_s, gas_org) #interpolate atmospheric gas history to model time.
        #bc_u, bc_d, bc_u_0 = self.boundaries(gas_org) #set boundary and initial conditions: bc_u is atmospheric condition, bc_d is zero gradient.        
  
        if (self.cc['UserData']):
            loadgas   = False
            f_depth   = np.loadtxt(os.path.join(self.dataDir,'depth.csv'), delimiter = ',', skiprows = 0)
            f_density = np.loadtxt(os.path.join(self.dataDir,'density.csv'), delimiter = ',', skiprows = 0)
            f_dcon    = np.loadtxt(os.path.join(self.dataDir,'Dcon.csv'), delimiter = ',', skiprows = 0)        
            f_temp    = np.loadtxt(os.path.join(self.dataDir,'temp.csv'), delimiter = ',', skiprows = 0)        
            f_clim    = np.loadtxt(os.path.join(self.dataDir,'Clim.csv'), delimiter = ',', skiprows = 0)
            f_BCO     = np.loadtxt(os.path.join(self.dataDir,'BCO.csv'), delimiter = ',', skiprows = 0)
            f_LID     = np.loadtxt(os.path.join(self.dataDir,'LID.csv'), delimiter = ',', skiprows = 0)
            f_gas     = np.loadtxt(os.path.join(self.dataDir,'GasHistory.csv'), delimiter = ',', skiprows = 0)
            
        
            self.write("Gasses loaded")
            
            Accu_vec = f_clim[:,1]
            T_vec    = f_clim[:,2]
            T_DX     = T_vec[0]
            
            time_yr    = f_clim[:,0]
            time_yr_s  = time_yr * self.sPerYear
            gas_org    = f_gas[jj + 1,:] #this needs to match up with the order of gasses specified in the config.
            
            stpsperyear = 1 / (time_yr[2] - time_yr[1])
            nt, model_time, model_time_years, dt  = self.time(time_yr, stpsperyear)
    
            if nt > len(time_yr): #this is because occasionally np.arrange will include end points.
                model_time = np.arange(time_yr[0] * self.sPerYear, time_yr[-1] * self.sPerYear, dt)
                model_time_years = model_time / self.sPerYear
                nt, model_time, model_time_years, dt  = np.size(model_time) #number of time steps        
                            
            Accu_vec = np.interp(model_time,time_yr,Accu_vec)
            T_vec = np.interp(model_time,time_yr,T_vec)
            D_x, M, deltaM, conc1, d_0, omega = MPG.gasses(gaschoice, sitechoice, T, p_a, self.dataDir, hemisphere, loadgas) #loadgas is False for userdata and true for not
            gas = np.interp(model_time, time_yr_s, gas_org) #interpolate atmospheric gas history to model time.
            bc_u, bc_d, bc_u_0 = self.boundaries(gas_org) #set boundary and initial conditions: bc_u is atmospheric condition, bc_d is zero gradient. 
            phi_0 = np.zeros(nz_P) #phi is mixing ratio of gas.
            phi_0[:] = bc_u_0
            rho_co, por, bcoRho, LIDRho = self.porosity(rho_prof, T) #get porosity
        
        else:    
            p_a, T, Accu_0, czd, z_co, LIZ, rho0, hemisphere = MPS.sites(sitechoice)   
            Accu_m = Accu_0 #Accumulation in m/year 
            Accu_0 = Accu_0 / self.sPerYear #accumulation in m/second
            
            loadgas = True
                    
            time_yr   = conc1[:,0] # Atmospheric measurements times   
            time_yr_s = time_yr * self.sPerYear  
            gas_org   = conc1[:,1] # Atmospheric measurement concentrations
            
            nt, model_time, model_time_years, dt = self.time(time_yr, stpsperyear)
            
            rhoHL = MPRHO.rhoHLAnalytic(self.R, T, self.rho_i, rho0, self.rho_bco, z_nodes, Accu_m) # Get density profile from H&L analytic
            rho_co, por, bcoRho, LIDRho = self.porosity(rhoHL, T)  
            rho_prof = rhoHL           

            D_x, M, deltaM, conc1, d_0, omega = MPG.gasses(gaschoice, sitechoice, T, p_a, self.dataDir, hemisphere, loadgas) #loadgas is False for userdata and true for not
            gas = np.interp(model_time, time_yr_s, gas_org) #interpolate atmospheric gas history to model time.
            bc_u, bc_d, bc_u_0 = self.boundaries(gas_org) #set boundary and initial conditions: bc_u is atmospheric condition, bc_d is zero gradient.   

        phi_0 = np.zeros(nz_P) #phi is mixing ratio of gas.
        phi_0[:] = bc_u_0
        rho_co, por, bcoRho, LIDRho = self.porosity(rho_prof, T) #get porosity
        
        # time = 0                 
        s          = (nt,nz_P)
        phi        = np.zeros(s)
        diffu_hold = np.zeros(s)
        rho_hold   = np.zeros(s)
        phi_t      = phi_0

        self.writeToFile(ConcPath, "w", model_time_years[i_time], phi_t)
        self.writeToFile(RhoPath, "w", model_time_years[i_time], rho_prof)
        self.writeToFile(DiffuPath, "w", model_time_years[i_time], diffu)
        self.writeToFile(DiffuPath, "w", model_time_years[i_time], z_nodes)  

        for i_time in range(1, nt): #6/18/14: need to fix nt so that it is consistent with output from firnmodel.py
    
            if self.cc['UserData']: #the 1: is because the first value in the csv rows is the time.
                rho_prof = np.interp(z_nodes, f_depth[i_time, 1:], f_density[i_time, 1:]) #density profile, interpolating onto consistent grid (z_nodes); 6/18/14: need to make sure that z_nodes go deep enough to track.       
                dconint  = np.interp(z_nodes, f_depth[i_time, 1:], f_dcon[i_time, 1:]) #this is the diffusivity constant
                Tprof    = np.interp(z_nodes, f_depth[i_time, 1:], f_temp[i_time, 1:])
                T        = Tprof[0]
                Accu_m   = Accu_vec[i_time] #Accumulation in m/year 
                Accu_0   = Accu_m / self.sPerYear #accumulation in m/second
                              
            z_co = min(z_nodes[rho_prof >= (bcoRho)]) #close-off depth; bcoRho is close off density
            LIZ  = min(z_nodes[rho_prof >= (LIDRho)]) #lock in depth; LIDRho is lock-in density
        
            a, a_P_0, b_0, dZ_d, S_P, dZ, S_C_0, w_edges, bubble_pres  = self.calculate(por, zz, diffu, d_eddy, deltaM, T, omega, Tprof, phi_t, rho_prof, Accu_0, p_a)
            
            bc_u_0  = gas[i_time] 
            bc_type = 1.
            bc_u    = np.array([bc_u_0, bc_type])
            
            b = b_0 + a_P_0 * phi_t       
            
            #Upper boundary
            a['a_P'][0] = 1 
            a['a_U'][0] = 0
            a['a_D'][0] = 0
            b[0]   = bc_u[0]
            
            #Down boundary
            a['a_P'][-1] = 1 
            a['a_D'][-1] = 0
            a['a_U'][-1] = 1
            b[-1]   = dZ_u[-1] * bc_d[0]
    
            phi_t = self.solver(a, b)
            
            phi[i_time,:] = phi_t
            diffu_hold[i_time,:] = diffu
            rho_hold[i_time,:] = rho_prof
            

            self.writeToFile(ConcPath, "a", model_time_years[i_time], phi_t)
            self.writeToFile(RhoPath, "a", model_time_years[i_time], rho_prof)
            self.writeToFile(DiffuPath, "a", model_time_years[i_time], diffu)
            
        return phi, diffu_hold, rho_hold, zz['z_nodes']
            
    def writeToFile(self, file, option, x, y):
        '''
        Write results to .csv files

        Arguments:
        file -- file name that the user wishes to save to
        option -- 'w' or 'a'
        x -- first column (time)
        y -- second column (data)
        '''
        with open(file, option) as f:
            writer = csv.writer(f)
            writer.writerow(np.append(x, y))

    def space(self, depth, z_res):
        '''

        Arguments:
        depth --
        z_res --

        Returns:
        dz --
        z_edges --
        z_nodes --
        nodes --
        nz_P --
        nz_fv --
        '''
        nz_fv = np.around(depth / z_res)
        nz_P  = nz_fv + 2
        
        dz      = depth / nz_fv
        z_edges = dz * np.arange(0, nz_fv + 1)
            
        z_nodes = np.concatenate(([z_edges[0]], z_edges[0:-1] + np.diff(z_edges) / 2, [z_edges[-1]]))
        #nodes   = np.size(z_nodes)
        
        zz = {}
        zz['dz'] = dz
        zz['z_edges'] = z_edges
        zz['z_nodes'] = z_nodes
        return zz, nz_P, nz_fv
        
    def porosity(self, rho_prof, T):
        '''

        Arguments:
        rho_prof --
        T --

        Returns:
        rho_co --
        por -- 
        bcoRho --
        LIDRho --
        '''   
        bcoRho = 1 / (1 / (self.rho_i) + T * 6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
        LIDRho = bcoRho - 14.0 #LIZ depth (Blunier and Schwander, 2000)
        
        ## Porosity, from Goujon et al., 2003, equations 9 and 10
        por_tot = 1 - rho_prof / self.rho_i # Total porosity
        rho_co  = bcoRho #use Martinerie close-off criteria
        #rho_co = 0.815 # User chosen close off-density (put in site-specific in sites?)
        por_co = 1 - rho_co / self.rho_i # Porosity at close-off
        alpha  = 0.37 # constant determined in Goujon
        por_cl = alpha * por_tot * (por_tot / por_co) ** (-7.6)
        ind    = por_cl > por_tot
        por_cl[ind] = por_tot[ind]
        
        por_op = por_tot - por_cl # Open Porosity
        por_op[por_op <= 0] = 1.e-10
        
        por = {}
        por['por_co'] = por_co
        por['por_tot'] = por_tot
        por['por_cl'] = por_cl
        por['por_op'] = por_op

        return rho_co, por, bcoRho, LIDRho
          
    def diffusivity(self, rho_co, por, z_co, czd, LIZ, d_0, D_x, T, p_a, z_nodes, rhoHL): #rhoprof is density profile
        ''' 

        Arguments:
        rho_co --
        por -- 
        z_co --
        czd --
        LIZ --
        d_0 --
        D_x --
        T --
        p_a --
        z_nodes --
        rhoHL --

        Returns:
        diffu --
        d_eddy -- 
        '''        
        #if rhoprof is None:
        #    rhoprof = rhoHL
        
        ## Constants
        d_eddy_sc = d_0 #Eddy diffusivity in the convective zone
        h = z_nodes
        dind = np.min(np.where(z_nodes>LIZ))
         
        ## Use Severinghaus relationship from Cuffey and Paterson
        d_0_sev = d_0 * 1.7
        #d_0_sev=d_0
        diffu_full_sev = D_x * d_0_sev * ((self.p_0 / p_a) * (T / self.T_0) ** 1.85 * (2.00 * (1 - (rhoHL / self.rho_i)) - 0.167))
        diffu_full_sev = diffu_full_sev - diffu_full_sev[dind]
        diffu_full_sev[diffu_full_sev <= 0] = 1e-15
        
        ## Use Schwander 1988, Eq. 2 Diffusivity (does not work very well) use 4e2
        ## for d_0
        k_sch = self.p_0 / p_a * (T / 253.16) ** 1.85 # Constant given in Schwander
        #diffu_full_sch =3.72*0.5*k_sch*(23.7*por_tot-2.84)*31.5 # Schwander' diffusivity relationship (for CO2). 31.5 is unit conversion. Added extra 3.72* 9/12/13
        diffu_full_sch = k_sch * (23.7 * por['por_tot'] - 2.84) / (1000 ** 2) # Schwander' diffusivity relationship (for CO2). 1/1000**2 is unit conversion. Added extra 3.72* 9/12/13
        diffu_full_sch = diffu_full_sch - diffu_full_sch[dind]
        diffu_full_sch[diffu_full_sch < 0] = 1.e-15
        
        ## Use Freitag, 2002, Eq 15 Diffusivity use 9e2 for d_0
        d_0_fre = d_0 * 4.9
        diffu_full_fre = D_x * d_0_fre * por['por_op'] ** 2.1
        diffu_full_fre = diffu_full_fre - diffu_full_fre[dind]
        diffu_full_fre[diffu_full_fre <= 0] = 1e-15
        
        sitechoice = self.cc['sitechoice']
        
        ## Use Christo's diffusivity data from NEEM-EU
        if sitechoice == 'NEEM':
            diffu_data = np.loadtxt(os.path.join(self.dataDir,'c_diffu_NEEM.txt'))
            h = diffu_data[:,0]
            diffu_full_data = D_x * d_0 * diffu_data[:,1]
    
        elif sitechoice == 'WAIS':
            diffu_data = np.loadtxt(os.path.join(self.dataDir,'c_diffu_WAIS.txt'))
            h = diffu_data[:,0]
            diffu_full_data = D_x * d_0 * diffu_data[:,1]
            
        elif sitechoice == 'SCENARIO':
            diffu_data = np.loadtxt(os.path.join(self.dataDir,'c_diffu_NEEM.txt'))
            h = diffu_data[:,0]
            diffu_full_data = D_x * d_0 * diffu_data[:,1]        

        diffu_full_data = np.interp(z_nodes, h, diffu_full_data)
        
        ## Add in high diffusivity in convective zone and low diffusivity below LIZ
        if self.cc["diffu"] == "data":
            diffu_full = diffu_full_data #change this line to change your choice of diffusivity
        
        if self.cc["diffu"] == "Sev":
            diffu_full = diffu_full_sev #change this line to change your choice of diffusivity
            
        if self.cc["diffu"] == "Fre":
            diffu_full = diffu_full_fre #change this line to change your choice of diffusivity
        
        if self.cc["diffu"] == "Sch":
            diffu_full = diffu_full_sch #change this line to change your choice of diffusivity
        
        #Add eddy diffusivity terms: convective zone and non-diffusive zone
        d_eddy = np.zeros(np.size(diffu_full))
        ind = np.nonzero(z_nodes < czd)
        d_eddy_surf = 2.426405E-5 #Kawamura, 2006
        H_scale = czd
        d_eddy_up = d_eddy_surf * np.exp(-1 * z_nodes / H_scale)

        ind  = np.flatnonzero(z_nodes > LIZ)
        ind2 = np.flatnonzero(z_nodes < z_co)
        ind3 = np.intersect1d(ind, ind2)

        d_eddy[ind3] = diffu_full[ind]
        d_eddy       = d_eddy + d_eddy_up

        diffu_full[ind] = 1e-15 #set molecular diffusivity equal to zero after LIZ - eddy diffusivity term drives diffusion below
        
        diffu = diffu_full
             
        return diffu, d_eddy # diffu_full_fre, diffu_full_sch, diffu_full_sev, diffu_full_data  
                
    def boundaries(self, gas_org):
        '''

        Arguments:
        gas_org --

        Returns:
        bc_u --
        bc_d --
        bc_u_0 --
        '''
        bc_u_0  = gas_org[0] #this is concentration on upper boundary (i.e. atmosphere) It is not actually necessary here, because it gets looped over.
        bc_type = 1
        bc_u    = np.concatenate(([bc_u_0], [bc_type]))
    
        bc_d_0  = 0 
        bc_type = 2
        bc_d    = np.concatenate(([bc_d_0], [bc_type]))
         
        return bc_u, bc_d, bc_u_0

    def calculate(self, por, zz, diffu, d_eddy, deltaM, T, omega, Tprof, phi_t, rho_prof, Accu_0, p_a, dt):    
        '''
        Arguments:
        por --
        diffu -- 
        d_eddy --
        deltaM --
        T --
        omega --
        Tprof --
        phi_t --
        rho_prof --
        Accu_0 --
        p_a --
        dt --

        Returns:
        a -- 
        a_P_0 --
        b_0 --
        dZ_d --
        S_P --
        dZ --
        S_C_0 --
        w_edges --
        bubble_pres --
        '''
        mode = self.cc['runtype']
        dz = zz['dz']
        por_op = por['por_op']

        diffu_P  = diffu * por_op
        d_eddy_P = d_eddy * por_op

        dZ   = np.concatenate(([1], np.diff(zz['z_edges']), [1]))
        
        dZ_u = np.diff(zz['z_nodes'])
        dZ_u = np.append(dZ_u[0], dZ_u)
    
        dZ_d = np.diff(zz['z_nodes'])
        dZ_d = np.append(dZ_d, dZ_d[-1])
    
        f_u = np.append(0, (1 -(zz['z_nodes'][1:] - zz['z_edges']) / dZ_u[1:]))
        f_d = np.append(1 - (zz['z_edges'] - zz['z_nodes'][0:-1]) / dZ_d[0:-1], 0)
        
        diffu_U = np.append(diffu_P[0], diffu_P[0:-1])
        diffu_D = np.append(diffu_P[1:], diffu_P[-1])
    
        diffu_u =  1/ ((1 - f_u) / diffu_P + f_u / diffu_U)
        diffu_d =  1/ ((1 - f_d) / diffu_P + f_d / diffu_D)
    
        d_eddy_U = np.append(d_eddy_P[0], d_eddy_P[0:-1])
        d_eddy_D = np.append(d_eddy_P[1:], d_eddy_P[-1])
    
        d_eddy_u =  1/ ((1 - f_u) / d_eddy_P + f_u / d_eddy_U)
        d_eddy_d =  1/ ((1 - f_d) / d_eddy_P + f_d / d_eddy_D)
        
        if self.cc['gravity']=="off" and self.cc['thermal']=="off":
            print 'gravity and thermal are off'
            S_C_0 = 0.0
    
        elif self.cc['gravity'] =='on' and self.cc['thermal']=='off':
            S_C_0 = (-diffu_d + diffu_u) * (deltaM * self.g / (self.R * T)) / dz #S_C is independent source term in Patankar
    
        elif self.cc['gravity'] == 'on' and self.cc['thermal'] == 'on' and mode == 'steady':
            print 'thermal on'
            dTdz = np.zeros(np.size(diffu_d))
            dTdz[0:100] = -0.0 #K/m. Negative gradient here means that it is colder deep (z is positive down)  
            S_C_0 = (diffu_d - diffu_u) * ((-deltaM * self.g / (self.R * T)) + (omega * dTdz)) / dz #S_C is independent source term in Patankar

        elif self.cc['gravity'] == 'on' and self.cc['thermal'] == 'on' and mode == 'transient':
            print "thermal on"
            if self.cc["UserData"]:
                dTdz = np.gradient(Tprof) / dz
            else:
                dTdz = np.ones(np.size(diffu_d))
                dTdz = self.cc["Tgrad"]
                
            S_C_0 = (diffu_d - diffu_u) * ((-deltaM * self.g / (self.R * T)) + (omega * dTdz)) / dz 
        
        S_C = S_C_0 * phi_t #this line might be the troublesome one! Should it be phi_0 instead?
        
        S_P = 0.0
        
        b_0 = S_C * dZ
        
        rho_interface = np.interp(zz['z_edges'], zz['z_nodes'], rho_prof)
        
        w_edges, bubble_pres = self.w(por, zz, Accu_0, rho_interface, T, p_a)

        w_u = np.append(w_edges[0], w_edges)
        w_d = np.append(w_edges, w_edges[-1])
        
        D_u = ((diffu_u + d_eddy_u) / dZ_u) 
        D_d = ((diffu_d + d_eddy_d) / dZ_d)
    
        
        F_u =  w_u * por_op
        F_d =  w_d * por_op
        
        P_u = F_u / D_u
        P_d = F_d / D_d
        
        a_U = D_u * self.A(P_u) + self.F_upwind( F_u)
        a_D = D_d * self.A(P_d) + self.F_upwind(-F_d)
    
        a_P_0 = por_op * dZ / dt
                                
        a_P =  a_U + a_D + a_P_0 - S_P * dZ
        a = {}
        a['a_U'] = a_U
        a['a_P'] = a_P
        a['a_D'] = a_D
           
        return a, a_P_0, b_0, dZ_d, S_P, dZ, S_C_0, w_edges, bubble_pres

    def time(self, time_yr, stpsperyear):
        '''

        Arguments: 
        time_yr --
        stpsperyear --

        Returns:
        nt --
        model_time --
        model_time_years --
        dt --
        '''
        yrs         = time_yr[-1] - time_yr[0]
        time_total  = yrs * self.sPerYear #total model run time in seconds
        t_steps     = yrs * stpsperyear
        dt          = time_total / t_steps #time step size.
        
        model_time  = np.arange(time_yr[0] * self.sPerYear, time_yr[-1] * self.sPerYear + dt, dt) #set model time steps        
        model_time_years = model_time / self.sPerYear
        nt = np.size(model_time) #number of time steps

        return nt, model_time, model_time_years, dt

    def firnair(self): #FIX GASCHOICE OPTIONS
        '''

        Returns:
        d --
        '''
        gaschoice_all = ["d15N2","CO2"]
        #gaschoice     = self.cc['gaschoice']
        nogas = len(gaschoice_all) #FIX -- right now, it's the length of the string
        d = {}
        for jj in xrange(nogas):
            print jj         
            gaschoice = gaschoice_all[jj]
            print jj, gaschoice      
            runtype = self.cc["runtype"]         
            if runtype == 'transient':
                phi, diffu_hold, rho_hold, z_nodes = self.FirnAir_TR(gaschoice, jj)
                #print 'maximum = %s' % np.max(phi)
                
            elif runtype == 'steady':
                phi, a_P_out, bubble_pres, z_nodes = self.FirnAir_SS(gaschoice)
                #print 'maximum = %s' % np.max(phi)
            
            print 'maximum = %s' % np.max(phi)
            d[gaschoice] = phi
            print '%s done' % gaschoice
            d['nodes'] = z_nodes           
        return d
    
    def write(self, txt): 
            self.text.insert(Tkinter.END, str(txt))

#if __name__ == "__main__":
    #firnAir = firnAirModel(sys.argv[1:], sys.argv[2:], sys.argv[3:]) #need to be changed (sys.argv[1:])
