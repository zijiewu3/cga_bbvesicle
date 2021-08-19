import crease_ga as cga
import numpy as np
import sys

class scatterer_generator:
    '''
    shape specific descriptors (shape_params):
    ------------------------------------------
    k:
        number of layers. Total number of input
        parameters need to be k*3+3.

    Input parameters to be predicted:
    ------------------------------------------
    for each k_i, excluding the center layer:
    rho_i:
        electron density of layer i
    sigma_i:
        sd of gaussian distribution ("width" of 
        layer).
    epsilon_i:
        mean of gaussian distribution ("position"
         of layer).

    for center layer:
    R_0:
        "radius" of the vesicle as defined in B&B
    sigma_center:
        "width" of center layer.

    bg:
        -log10(bg)
    '''
    def __init__(self,k = 3,
                 minvalu = (-10,30,-200,-10,30,30,110,30,0.1),
                 maxvalu = (10,200,-30,10,200,200,700,200,4)):
        self.k = k
        self._numvars = (self.k-1)*3+3
        self.minvalu = minvalu
        self.maxvalu = maxvalu

    @property
    def numvars(self):
        return self._numvers
    
    def converttoIQ(self,qrange,param):
        k = self.k
        if (len(param)-3)/3 != k-1:
            print("Total number of parameters is incorrect")
            sys.exit(1)
        rhos = np.zeros(k)
        sigmas = np.zeros(k)
        epsilons = np.zeros(k)

        for i in range(k-1):
            rhos[i] = param[i*3]
            sigmas[i] = params[i*3+1]
            epsilons[i] = params[i*3+2]


        bg = param[-1]
        sigma_center = param[-2]
        R0 = param[-3]
        
        rhos[-1] = -1
        sigmas[-1] = sigma_center
        epsilons[-1] = 0

        IQid = np.zeros((len(qrange)))
        for qi,q in enumerate(qrange):
            for k1 in range(k):
                for k2 in range(k1,k):
                    itemsp = (R0+epsilon[k1])*(R0+epsilon[k2])
                    item1 = rho[k1] * rho[k2] * sigma[k1] * sigma[k2]
                    item2 = np.exp(-q**2*(sigma[k1]**2+sigma[k2]**2)/2)
                    item3 = np.cos(q*(epsilon[k1]-epsilon[k2]))
                    IQid[qi] += q**(-2)*itemsp*item1*item2*item3
        IQid /= IQid[0]
        for i in IQid:
            if i < 0:
                #if a negative value shows up, effectively throw this
                #individual away
                IQid = 100*np.ones(len(qrange))
        IQid += 10**(-bg)

        return IQid
        #Plug into the equation
