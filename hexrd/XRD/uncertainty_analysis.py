# ============================================================
# Copyright (c) 2012, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. 
# Written by Joel Bernier <bernier2@llnl.gov> and others. 
# LLNL-CODE-529294. 
# All rights reserved.
# 
# This file is part of HEXRD. For details, see https://github.com/joelvbernier/hexrd.
# 
# Please also see the file LICENSE.
# 
# This program is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free Software
# Foundation) version 2.1 dated February 1999.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of the 
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this program (see file LICENSE); if not, write to
# the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA or visit <http://www.gnu.org/licenses/>.
# ============================================================
"""
follows Nonlinear parameter estimation, Bard.
"""
import numpy as num
from scipy.linalg import inv,eig
import scipy.optimize
import numpy.random
import scipy.stats

def findChisqr(conf_level, df, x0 = 16.):
    """
    uses scipy.optimize.fsolve to find the chisquared value for a particular level of confidence.
    """
    def _chisqrdiff(chisq, df, target):
        return scipy.stats.chisqprob(chisq,df) - target
    target_chisqr = 1 - conf_level
    chisqr = scipy.optimize.fsolve(_chisqrdiff, x0 = x0, args = (df, target_chisqr))
    return chisqr

def computeAllparameterUncertainties(params,cov_x,confidence_level):
    u_is = []
    for i in range(len(params)):
        u_i = parameterUncertainty(params, i, cov_x, confidence_level)
        u_is.append(u_i)
    return num.array(u_is)
def parameterUncertainty(params, index, cov_x, confidence_level):
    """
    computes uncertainty in parameters using covariance matrix. See Bard p 188.
    example Bard eq (7-10-6)
    \begin{equation}
    Pr[(\hat{\theta} - \theta^*)^T V^{-1} (\hat{\theta} - \theta^*) < 15.987] = .9
    \end{equation}
    $\theta^*, \hat{\theta}$ tested objective minimizer, true objective minimizer
    
    Here 0.9 is the confidence level, 15.987 is the chisquared value corresponding to this confidence level.
    the quantity
    \[
    (\hat{\theta} - \theta^*)^T V^{-1} (\hat{\theta} - \theta^*)
    \]
    is distributed as $\chi^2$ with $l$ degress of freedom, where $l$ is the number of parameters $\theta$.
    Technically only valid for normal, unbiased distribution.
        
    Assumes the inverse covariance matrix is a good estimate for the Hessian of the objective function, good for normally distributions, where the Hessian H is
    H_{ij}(\theta) = \frac{\partial^2 \Phi}{\partial \theta^2}
    """
    def _bilinear_distance(step, direction, A, target):        
        x = step*direction 
        Ax = num.dot(A,x)
        return num.dot(Ax,x) - target
    df = len(params)
    direction = num.zeros(df)
    direction[index]=1
    cov_x_inv = inv(cov_x)
    chi_sqr = findChisqr(confidence_level, df)    
    sol = scipy.optimize.fsolve(_bilinear_distance, x0 = 0.001, args = (direction,cov_x_inv,chi_sqr), warning = False)
    return sol

def propagateUncertainty(model,uncertainties, step=1e-8,*args):
    #print 'propUncertainty',model,uncertainties,args
    parameters = args[0]
    #print 'args',args
    assert(len(parameters)==len(uncertainties))
    #extra_args = ()
    #if len(args)==2:
    extra_args = args[1:len(args)]
    #elif len(args)>2:
    #    raise Exception, 'propagateUncertainty has extra arguments for function'
    sum_squares = 0.
    for i in range(len(parameters)):
        dfdx = directionalDerivative(model,i,step,parameters, *extra_args)
        ux = uncertainties[i]
        sum_squares += (dfdx*ux)**2
    uf = num.sqrt(sum_squares)
    return uf
        
def directionalDerivative(model,index,step = 1e-8,*args):
    #print 'directDerivative',model,args
    parameters = args[0]
    #extra_args = ()
    extra_args = args[1:len(args)]
    #if len(args)==2:
    #    extra_args = args[1]
    #elif len(args)>2:
    #    raise Exception, 'directionalDerivative has extra arguments for function'
        
    direction = num.zeros(len(parameters))
    direction[index]=1

    f0 = model(parameters, *extra_args)
    f1 = model(parameters+direction*step, *extra_args)
    #print f0,f1
    df_dp = (f1-f0)/step
    return df_dp


def getLsqUncertainty(lsq_func,conf_level ,*lsq_func_args, **lsq_func_kwargs):
    out = lsq_func(*lsq_func_args, **lsq_func_kwargs)
    cov = out[1]
    params = out[0]
    uparams = computeAllparameterUncertainties(params,cov,conf_level)
    return uparams, params, cov


def scan_initial_condition_3space(estimator, kwname, seedpt, delta, nsteps, filename, *estimator_args, **estimator_kwargs):
    """will need external wrapper"""
    
    axis1 = num.linspace(-delta[0], delta[0], nsteps)
    axis2 = num.linspace(-delta[1], delta[1], nsteps)
    axis3 = num.linspace(-delta[2], delta[2], nsteps)
    all_params = []
    all_params_dict = {}
    for p1 in axis1:
        print 'p1,',p1
        for p2 in axis2:
            for p3 in axis3:
                testpt = seedpt + num.array([p1,p2,p3])
                out = wrapf(estimator, testpt, kwname, *estimator_args, **estimator_kwargs)
                params = out[0]
                all_params.append(params)
                all_params_dict[(p1,p2,p3)] = params
    all_params_ = num.array(all_params)
    meansolution = num.mean(all_params_, 0)
    cov_ = num.cov(all_params_, rowvar = 0)
    from File_funcs import write_lines
    from Vector_funcs import Mag
    f = open(filename,'w')
    for key in all_params_dict.keys():
        p1,p2,p3 = key
        param = Mag(all_params_dict[key] - meansolution)
        data = []
        data.extend(key)
        data.extend(param)
        write_lines(f,data)
    f.close()
    print 'wrote', f.name
    pass
def scan_initial_condition_12space(estimator, kwname, seedpt, delta, nsteps, filename, *estimator_args, **estimator_kwargs):
    """will need external wrapper

    scans the estimator function given a variety of initial conditions
    scalar output for each initial condition is the two norm difference between that particular result and the
    
    import uncertainty_analysis
    import grain_wrap as gw
    #have a bunch of grains
    for nf in rubyGrains.keys():
        rgrain = rubyGrains[nf]
        rubyGrain_wrap = gw.grainStrainAnalysis(rgrain)
        fname = "fitRotationStar_Angle_weights_initial_cond_plot" +str(nf)
        uncertainty_analysis.scan_initial_condition_3space(rubyGrain_wrap.fitRotationStar_angle_weights, "Rparams", num.array([0,0,0], dtype = float), delta = num.array([1.2, 1.2, 1.2]), nsteps = nsteps, filename = fname , Uparams = 'default', report = False)



    """
    
    axis1 = num.linspace(-delta[0], delta[0], nsteps)
    axis2 = num.linspace(-delta[1], delta[1], nsteps)
    axis3 = num.linspace(-delta[2], delta[2], nsteps)
    axis4 = num.linspace(-delta[3], delta[3], nsteps)
    axis5 = num.linspace(-delta[4], delta[4], nsteps)
    axis6 = num.linspace(-delta[5], delta[5], nsteps)
    axis7 = num.linspace(-delta[6], delta[6], nsteps)
    axis8 = num.linspace(-delta[7], delta[7], nsteps)
    axis9 = num.linspace(-delta[8], delta[8], nsteps)
    axis10 = num.linspace(-delta[9], delta[9], nsteps)
    axis11 = num.linspace(-delta[10], delta[10], nsteps)
    axis12 = num.linspace(-delta[11], delta[11], nsteps)

        
    
    all_params = []
    all_params_dict = {}
    for p1 in axis1:
        print 'p1,',p1
        for p2 in axis2:
            for p3 in axis3:
                for p4 in axis4:
                    for p5 in axis5:
                        for p6 in axis6:
                            for p7 in axis7:
                                for p8 in axis8:
                                    for p9 in axis9:
                                        for p10 in axis10:
                                            for p11 in axis11:
                                                for p12 in axis12:
                                                    testpt = seedpt + num.array([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12])
                                                    out = wrapf(estimator, testpt, kwname, *estimator_args, **estimator_kwargs)
                                                    params = out[0]
                                                    all_params.append(params)
                                                    all_params_dict[(p1,p2,p3)] = params
    
    all_params_ = num.array(all_params)
    meansolution = num.mean(all_params_, 0)
    cov_ = num.cov(all_params_, rowvar = 0)
    from File_funcs import write_lines
    from Vector_funcs import Mag
    f = open(filename,'w')
    for key in all_params_dict.keys():
        p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12 = key
        param = Mag(all_params_dict[key] - meansolution)
        data = []
        data.extend(key)
        data.extend(param)
        write_lines(f,data)
    f.close()
    print 'wrote', f.name
   
def scan_initial_condition_6space(estimator, kwname, seedpt, delta, nsteps, filename, *estimator_args, **estimator_kwargs):
    """will need external wrapper

    scans the estimator function given a variety of initial conditions
    scalar output for each initial condition is the two norm difference between that particular result and the
    
    import uncertainty_analysis
    import grain_wrap as gw
    #have a bunch of grains
    for nf in rubyGrains.keys():
        rgrain = rubyGrains[nf]
        rubyGrain_wrap = gw.grainStrainAnalysis(rgrain)
        fname = "fitRotationStar_Angle_weights_initial_cond_plot" +str(nf)
        uncertainty_analysis.scan_initial_condition_3space(rubyGrain_wrap.fitRotationStar_angle_weights, "Rparams", num.array([0,0,0], dtype = float), delta = num.array([1.2, 1.2, 1.2]), nsteps = nsteps, filename = fname , Uparams = 'default', report = False)



    """
    
    axis1 = num.linspace(-delta[0], delta[0], nsteps)
    axis2 = num.linspace(-delta[1], delta[1], nsteps)
    axis3 = num.linspace(-delta[2], delta[2], nsteps)
    axis4 = num.linspace(-delta[3], delta[3], nsteps)
    axis5 = num.linspace(-delta[4], delta[4], nsteps)
    axis6 = num.linspace(-delta[5], delta[5], nsteps)
    
    all_params = []
    all_params_dict = {}
    for p1 in axis1:
        print 'p1,',p1
        for p2 in axis2:
            for p3 in axis3:
                for p4 in axis4:
                    for p5 in axis5:
                        for p6 in axis6:
                            testpt = seedpt + num.array([p1,p2,p3,p4,p5,p6])
                            out = wrapf(estimator, testpt, kwname, *estimator_args, **estimator_kwargs)
                            params = out[0]
                            all_params.append(params)
                            all_params_dict[(p1,p2,p3)] = params
    
    all_params_ = num.array(all_params)
    meansolution = num.mean(all_params_, 0)
    cov_ = num.cov(all_params_, rowvar = 0)
    from File_funcs import write_lines
    from Vector_funcs import Mag
    f = open(filename,'w')
    for key in all_params_dict.keys():
        p1,p2,p3,p4,p5,p6 = key
        param = Mag(all_params_dict[key] - meansolution)
        data = []
        data.extend(key)
        data.extend(param)
        write_lines(f,data)
    f.close()
    print 'wrote', f.name
      
def wrapf(func, arg1, kwname, *args, **kwargs):
    kwargs[kwname]=arg1
    print 'kwargs',kwargs
    return func(*args,**kwargs)
def wrapf_slot1(arg1, func, kwname, *args, **kwargs):
    kwargs[kwname]=arg1
    #print 'kwargs',kwargs
    return func(*args,**kwargs)
    
    
                
    

""" rest is misc for now """
#class model_func:
#    def __init__(self):
#        self.params = None
#        self.
def scalar_model_function_data(model,parameters, x_array):
    """
    model - scalar valued function, with element of x_array as independant input variables.
    """
    y_array = num.zeros(len(x_array))
    for i in range(len(x_array)):
        x = x_array[i]
        y = model(parameters,x)
        y_array[i]=y
    return y_array
def scalar_model_function_leastsq(parameters,model,x,y):
    difference = scalar_model_function_data(model,parameters,x)-y
    return difference

def scalar_model_function_data_noise(model,parameters, x_array, noisefactor = .1):
    import random
    y_array = num.zeros(len(x_array))
    for i in range(len(x_array)):
        x = x_array[i]
        y = model(parameters,x)
        y_array[i]=y+numpy.random.normal(0)*noisefactor*y
    return y_array


def Objective_Sum_Squares_Scalar(model_func,parameters,x,y):
    """
    x: independent variables into model function (list of numpy array)
    y: measured data (numpy array)
    """
    difference = scalar_model_function_data(model_func,parameters,x)-y
    residual = difference*difference
    return num.sum(residual)

def Hessian_Func(obj_func,model_func,params, x, y, thetasp = 10**-8):
    l = len(params)
    jacobian = num.zeros([l,l])
    for i in range(l):
        for j in range(l):
            dparami = num.zeros(l); dparamj = num.zeros(l)
            dparami[i] = 1        ; dparamj[j] = 1
            dPhi_dj = (obj_func(model_func,params+dparamj*thetasp,x,y) - obj_func(model_func,params,x,y))/thetasp
            d2Phi_didj = ((obj_func(model_func,params+dparamj*thetasp+dparami*thetasp,x,y) - obj_func(model_func,params + dparami*thetasp,x,y))/thetasp - dPhi_dj)/thetasp
            jacobian[i,j]=d2Phi_didj
    return jacobian

def d2Phi_dtheta_dw_func(obj_func,model_func,params,x,y,thetasp = 10**-8, wsp = 10**-8):
    columns = len(y)
    rows = len(params)
    temp = num.zeros([rows,columns])
    for i in range(rows):
        for j in range(columns):
            dw = num.zeros(columns)
            dtheta = num.zeros(rows)
            dw[j]=1
            dtheta[i]=1
            dPhi_dw = (obj_func(model_func,params,x,y+dw*wsp) - obj_func(model_func,params,x,y))/wsp
            d2Phi_dtheta_dw = ((obj_func(model_func, params+dtheta*thetasp,x,y+dw*wsp) - obj_func(model_func, params+dtheta*thetasp,x,y))/wsp - dPhi_dw)/thetasp
            temp[i,j]=d2Phi_dtheta_dw
    return temp
            

def Covariance_Matrix(obj_func, model_func, params, x, y, thetasp = 10**-8, wsp = 10**-8):
    H = Hessian_Func(obj_func, model_func, params,x,y,thetasp = thetasp)
    d2func = d2Phi_dtheta_dw_func(obj_func, model_func,params,x,y, thetasp = thetasp, wsp = wsp)
    Hinv = inv(H)
    V = num.dot(Hinv, d2func)
    V = num.dot(V,d2func.T)
    V = num.dot(V,Hinv)
    return V
"""
if __name__ == '__main__':
    'example'
    from uncertainty_analysis import *
    import numpy as num
    from scipy.linalg import eig
    import random as rand
    
    def example_model_function(parameters,x):
        p1,p2,p3 = parameters
        x1,x2,x3 = x
        return x1**p1+x2**p2+x3**p3
    
    def construct_position_data(numentries):
        temp = []
        for i in range(numentries):
            temp.append(num.array([rand.random(),rand.random(),rand.random()])*10)
        return temp
    
    position_data = construct_position_data(10)
    
    params = num.array([2,2,2])
    
    model_data = scalar_model_function_data(example_model_function,params, position_data)
    simulated_data = scalar_model_function_data_noise(example_model_function,params,position_data, noisefactor = .01)
    simulated_data2 = scalar_model_function_data_noise(example_model_function,params,position_data, noisefactor = .01)
    simulated_data3 = scalar_model_function_data_noise(example_model_function,params,position_data, noisefactor = .01)
    
    print Objective_Sum_Squares_Scalar(example_model_function,params,position_data,simulated_data)
    print Objective_Sum_Squares_Scalar(example_model_function,params,position_data,model_data)
    
    
    x = position_data[:]
    y = simulated_data
    y = model_data
    H = Hessian_Func(Objective_Sum_Squares_Scalar, example_model_function,params,x,y)
    d2func = d2Phi_dtheta_dw_func(Objective_Sum_Squares_Scalar, example_model_function,params,x,y)
    v_theta = Covariance_Matrix(Objective_Sum_Squares_Scalar, example_model_function, params, x, y)
    
    #find_best_params
    from scipy.optimize import leastsq
    out = leastsq(scalar_model_function_leastsq,x0 = (2,2,2), args = (example_model_function,x,y),full_output = 1)
    pfinal = out[0]
    cov = out[1]
    v_theta = Covariance_Matrix(Objective_Sum_Squares_Scalar, example_model_function, pfinal, x, y,10**-6,10**-6)
    #note: v_theta ~= cov
    print v_theta, cov
    #H = Hessian_Func(Objective_Function, a_multiple_data_scalar_model_function,params,x,model_data)
    #
    H = Hessian_Func(Objective_Sum_Squares_Scalar, example_model_function,pfinal,x,y)
    
    conf_level = .9 #90 percent
    
    sol0 = parameterUncertainty(params, 0, v_theta, conf_level)
    sol1 = parameterUncertainty(params, 1, v_theta, conf_level)
    sol2 = parameterUncertainty(params, 2, v_theta, conf_level)
    print sol0,'\n',sol1,'\n',sol2
    
    print 'uncertainties',computeAllparameterUncertainties(params, v_theta, conf_level)
    
    sol0 = parameterUncertainty(pfinal, 0, v_theta, conf_level)
    sol1 = parameterUncertainty(pfinal, 1, v_theta, conf_level)
    sol2 = parameterUncertainty(pfinal, 2, v_theta, conf_level)
    
    Objective_Sum_Squares_Scalar(example_model_function, pfinal,x,y)
    pedge0 = pfinal + [sol0,0,0]
    pedge1 = pfinal + [0,sol1,0]
    pedge2 = pfinal + [0,0,sol2]
    print Objective_Sum_Squares_Scalar(example_model_function, pedge0,x,y)
    print Objective_Sum_Squares_Scalar(example_model_function, pedge1,x,y)
    print Objective_Sum_Squares_Scalar(example_model_function, pedge2,x,y)
    
    
"""
