#!/usr/bin/env python

"""
A class to implement Simultaneous Perturbation Stochastic Approximation.
"""
# from lib2to3.pytree import BasePattern
from xmlrpc.client import TRANSPORT_ERROR
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
from sympy import residue
from ESP_Basic_Model import Water_curve
from ESP_Class_2022 import *
# from ESP_Basic_Model import *
# from pandas.core import api
# from sympy import E

# from sklearn.model_selection import train_test_split
# import multiprocessing as mp
import time
# import datetime


num_cores = 40      # this computer has 40 core cpu
each_core = 20

base_pump = 'GC6100'
pump_name = 'P100'
dpi = 300

QBEM_default = {'TE2700': 4500, 'DN1750': 3000, 'GC6100': 8543.69, 'P100': 11000, 'Flex31': 5000}    # ori
QBEM_default = {'TE2700': 4398, 'DN1750': 4630, 'GC6100': 8668, 'P100': 11498, 'Flex31': 6193}    
QBEM_base = 6000

# factor = 1.2 # factor used in gas bubble size
# DB_GVF_EFF = 1.02    ## factor used in gas bubble size
# DB_NS_EFF = 2    ## factor used in gas bubble size
# alphaG_crit_coef = 0.5
# alphaG_crit_critical = 0.5
# CD_Gas_EFF = 0.4    # Original drag coefficiet for gas velocity in modified Sun for INT flow
# CD_Liquid_EFF = 0.4    # Original drag coefficiet for gas velocity in modified Sun for INT flow
# CD_INT_EFF = 2
# CD_GV_EFF = 2 
train_type = 'ALL_geometry'
curve_type = 'SGL'
Coef_list = {
                "FTI": 3, "FTD": 3, "FTI_coef": 0.1*1.00169619e-02, "FTD_coef": 0.1*7.27961770e-01, "F_leakage":0.25, "SGMU_coef":175*3.35147357e+01, 
                "SGM_coef":0.1 , "SGM_coef_2":0.25, "SGMU_coef_2": 0.01, "QBEM_VISL_coef":0.01*5.76383343e+01,
                "emulsion_E": 3.0, "emulsion_SN": 0.01, "emulsion_WE":0.1, "emulsion_RE": 0.1, "emulsion_ST1": 2.5, "emulsion_ST2":0.2,
                "alphaG_crit_coef":4,"alphaG_crit_critical":0.5, "alphaG_crit_critical2":0.25,  "factor": 1.4, "DB_GVF_EFF": 1.02, "DB_NS_EFF": 2, "CD_Gas_EFF": 0.62,    
                "CD_Liquid_EFF": 1, "CD_INT_EFF": 1.5, "CD_GV_EFF": 2, "SGL":'zhu_2018', "GL":'None', "flg": 'F', "train_type": 'ALL_geometry'
            }

Coef_list = {
                "FTI": 3, "FTD": 3, "FTI_coef": 0.01, "FTD_coef": 0.01, "F_leakage":0.25, "SGMU_coef":175, 
                "SGM_coef":0.1 , "SGM_coef_2":0.25, "SGMU_coef_2": 0.5, "QBEM_VISL_coef":0.01,
                "emulsion_E": 3.0, "emulsion_SN": 0.01, "emulsion_WE":0.1, "emulsion_RE": 0.1, "emulsion_ST1": 2.5, "emulsion_ST2":0.2,
                "alphaG_crit_coef":4,"alphaG_crit_critical":0.5, "alphaG_crit_critical2":0.25,  "factor": 1.4, "DB_GVF_EFF": 1.02, "DB_NS_EFF": 2, "CD_Gas_EFF": 0.62,    
                "CD_Liquid_EFF": 1, "CD_INT_EFF": 1.5, "CD_GV_EFF": 2, "SGL":'zhu_2018', "GL":'None', "flg": 'F', "train_type": 'ALL_geometry'
            }


ESP_default = {'TE2700':
            {
                "R1": 0.017496,     "R2": 0.056054,     "TB": 0.00272,      "TV": 0.00448,      "RD1": 0.056,
                "RD2": 0.017496,    "YI1": 0.012194,    "YI2": 0.007835,    "VOI": 0.000016119, "VOD": 0.000011153,
                "ASF": 0.00176464,  "ASB": 0.00157452,  "AB": 0.001319,     "AV": 0.001516,     "ADF": 0.001482,
                "ADB": 0.000935,    "LI": 0.076,        "LD": 0.08708,      "RLK": 0.056209,    "LG": 0.00806,
                "SL": 0.00005,      "EA": 0.000254,     "ZI": 5,            "ZD": 9,            "B1": 19.5,
                "B2": 24.7,         "NS": 1600,         "DENL": 1000,       "DENG": 11.2,       "DENW": 1000,       "VISL": 0.001,
                "VISG": 0.000018,   "VISW": 0.001,      "ST": 0.073,        "N": 3500,          "SGM": 0.3,
                "QL": 2700,         "QG": 500,           "GVF": 10,          "WC": 0.,           "SN": 1,
                "P": 150,           "T":100,            "SGL":'zhu_2018', "GL":'None',        "flgz": 'F',
                "train_type": 'ALL_geometry'
            },
    'GC6100':
        {
                "R1": 0.027746,     "R2": 0.050013,     "TB": 0.0019862,    "TV": 0.002894,     "RD1": 0.0547,
                "RD2": 0.017517,    "YI1": 0.017399,    "YI2": 0.013716,    "VOI": 1.512E-5,    "VOD": 1.9818E-5,
                "ASF": 8.9654E-4,   "ASB": 9.4143E-4,   "AB": 1.0333E-3,    "AV": 1.769E-3,     "ADF": 2.0486E-3,
                "ADB": 1.0301E-3,   "LI": 0.0529,       "LD": 0.0839,       "RLK": 4.35E-2,   "LG": 0.0015475,
                "SL": 3.81E-4,      "EA": 0.0003,       "ZI": 7,            "ZD": 8,            "B1": 33.375,
                "B2": 41.387,       "NS": 3220,         "DENL": 1000,       "DENG": 11.2,       "DENW": 1000,       "VISL": 0.001,
                "VISG": 0.000018,   "VISW": 0.001,      "ST": 0.073,        "N": 3600,          "SGM": 0.3,
                "QL": 6100,         "QG": 50,           "GVF": 10,          "WC": 0.,           "SN": 1,
                "P": 250,           "T":100,            "SGL":'zhu_2018', "GL":'None',        "flgz": 'F'
        },
        # original R1 = 0.014351   TB = 0.0025896   "RLK": 6.1237E-2, 
        # new "RLK": 4.35E-2


        'Flex31':
            {
                "R1": 0.018,        "R2": 0.039403,     "TB": 0.0018875,    "TV": 0.0030065,    "RD1": 0.042062,
                "RD2": 0.018841,    "YI1": 0.015341,    "YI2": 0.01046,     "VOI": 0.000010506, "VOD": 0.000010108,
                "ASF": 6.888e-4,    "ASB": 6.888e-4,    "AB": 6.888e-4,     "AV": 8.626e-4,     "ADF": 8.626e-4, 
                "ADB": 8.626e-4,    "LI": 0.04252,      "LD": 0.060315,     "RLK": 0.033179,    "LG": 0.005, 
                "SL": 0.000254,     "EA": 0.0003,       "ZI": 6,            "ZD": 8,            "B1": 21.99,        
                "B2": 57.03,        "NS": 2975,         "DENL": 1000,       "DENG": 11.2,       "DENW": 1000,       "VISL": 0.001,      
                "VISG": 0.000018,   "VISW": 0.001,      "ST": 0.073,        "N": 3600,          "SGM": 0.3,
                "QL": 4000,         "QG": 50,           "GVF": 10,          "WC": 0.,           "SN": 1,
                "P": 160,           "T":100,            "SGL":'zhu_2018', "GL":'None',        "flgz": 'F'
            },
        # ori   "ASF": 6.888e-4,           "ASB": 6.888e-4,           "AB": 0.0020663,    "AV": 0.0025879,    "ADF": 0, 

        'DN1750':
        {
                "R1": 1.9875E-2,    "R2": 3.5599E-2,    "TB": 1.7E-3,       "TV": 3.12E-3,      "RD1": 0.04,
                "RD2": 0.01674,     "YI1": 1.3536E-2,   "YI2": 7.13E-3,     "VOI": 6.283E-6,    "VOD": 7.063E-6,
                "ASF": 6.8159E-04,  "ASB": 6.549E-04,   "AB": 6.9356E-04,   "AV": 7.1277E-04,   "ADF": 1.0605E-03,
                "ADB": 6.6436E-04,  "LI": 0.039,        "LD": 5.185E-02,    "RLK": 0.04,        "LG": 0.01,
                "SL": 0.00005,      "EA": 0.000254,     "ZI": 6,            "ZD": 8,            "B1": 20.3,
                "B2": 36.2,         "NS": 2815,         "DENL": 1000,       "DENG": 11.2,       "DENW": 1000,"VISL": 0.001,
                "VISG": 0.000018,   "VISW": 0.001,      "ST": 0.073,        "N": 3500,          "SGM": 0.3,
                "QL": 1750,         "QG": 50,           "GVF": 10,          "WC": 0.,           "SN": 3,
                "P": 150,           "T":100,            "SGL":'zhu_2018', "GL":'None',        "flgz": 'F'
        },

    'P100':
        {
                "R1": 0.023793966,  "R2": 0.05097,      "TB": 0.0023,       "TV": 0.00284,      "RD1": 0.052424,
                "RD2": 0.025349,    "YI1": 0.02315,     "YI2": 0.01644,     "VOI": 2.9E-5,      "VOD": 2.61E-5,
                "ASF": 0.00083,     "ASB": 0.001277143, "AB": 0.00186,      "AV": 0.00224,      "ADF": 0.002506,
                "ADB": 0.00142,     "LI": 0.04336,      "LD": 0.0810175,    "RLK": 0.045465,    "LG": 0.009605,
                "SL": 0.0001,       "EA": 0.000254,     "ZI": 7,            "ZD": 7,            "B1": 35.315,
                "B2": 38.17,        "NS": 3448,         "DENL": 1000,       "DENG": 11.2,       "DENW": 1000,       "VISL": 0.001,
                "VISG": 0.000018,   "VISW": 0.001,      "ST": 0.073,        "N": 3600,          "SGM": 0.3,
                "QL": 9000,         "QG": 50,           "GVF": 10,          "WC": 0.,           "SN": 1,
                "P": 150,           "T":100,            "SGL":'zhu_2018', "GL":'None',        "flgz": 'F'
        }
        
    }

class SimpleSPSA ( object ):
    """Simultaneous Perturbation Stochastic Approximation. 
    """
    # These constants are used throughout
    alpha = 0.602
    gamma = 0.101
    
    

    def __init__ ( self, loss_function, a_par = 1e-6, noise_var=0.01, args=(), \
            min_vals=None, max_vals=None, param_tolerance=None, \
            function_tolerance=None, max_iter=5000000,residual=1e-8 ):
        """The constructor requires a loss function and any required extra 
        arguments. Optionally, boundaries as well as tolerance thresholds can
        be specified.
        
        :param loss_function: The loss (or cost) function that will be minimised.
            Note that this function will have to return a scalar value, not a 
            vector.
        :param a_par: This is the ``a`` parameter, which controls the scaling of
            the gradient. It's value will have to be guesstimated heuristically.
        :param noise_var: The noise variance is used to scale the approximation
            to the gradient. It needs to be >0.
        :param args: Any additional arguments to ``loss_function``.
        :param min_vals: A vector with minimum bounds for parameters
        :param max_vals: A vector with maximum bounds for parameters
        :param param_tolerance: A vector stating the maximum parameter change
            per iteration.
        :param function_tolerance: A scalar stating the maximum change in 
            ``loss_function`` per iteration.
        :return: None
        """
        self.args = args
        self.loss = loss_function
        self.min_vals = min_vals
        self.max_vals = max_vals
        self.param_tolerance = param_tolerance
        self.function_tolerance = function_tolerance
        self.c_par = noise_var
        self.max_iter = max_iter
        self.big_a_par = self.max_iter/10.
        self.a_par = a_par
        self.residual=residual
        
    def calc_loss ( self, theta ):
        """Evalute the cost/loss function with a value of theta"""
        retval = self.loss ( theta, *(self.args ) )
        return retval

    def minimise ( self, theta_0, ens_size=2, report=500 ):
        """The main minimisation loop. Requires a starting value, and optionally
        a number of ensemble realisations to estimate the gradient. It appears
        that you only need two of these, but the more the merrier, I guess.
        
        :param theta_0: The starting value for the minimiser
        :param ens_size: Number of relaisations to approximate the gradient.
        :return: A tuple containing the parameters that optimise the function,
            the function value, and the number of iterations used.
        """
        n_iter = 0
        num_p = theta_0.shape[0]
        print ("Starting theta=", theta_0)
        theta = theta_0
        # theta1 = theta
        # theta2 = theta
        # theta3 = theta
        j_old = self.calc_loss ( theta )
        # j1 = j_old
        # j2 = j_old
        # j3 = j_old
        # i1 = 0
        # i2 = 0
        # i3 = 0
        j_ori = self.calc_loss ( theta_0 )
        # Calculate the initial cost function
        theta_saved = theta_0*100
        j_list = []
        j_list.append(j_old)
        while  (np.linalg.norm(theta_saved-theta)/np.linalg.norm(theta_saved) >\
                self.residual) and (n_iter < self.max_iter):
            # The optimisation carried out until the solution has converged, or
            # the maximum number of itertions has been reached.
            theta_saved = theta # Store theta at the start of the iteration
                                # as we may well be restoring it later on.
            # Calculate the ak and ck scalars. Note that these require
            # a degree of tweaking
            ak = self.a_par/( n_iter + 1 + self.big_a_par)**self.alpha
            ck = self.c_par/( n_iter + 1 )**self.gamma  
            ghat = 0.  # Initialise gradient estimate
            for j in np.arange ( ens_size ):
                # This loop produces ``ens_size`` realisations of the gradient
                # which will be averaged. Each has a cost of two function runs.
                # Bernoulli distribution with p=0.5
                delta = (np.random.randint(0, 2, num_p) * 2 - 1)
                # Stochastic perturbation, innit
                theta_plus = theta + ck*delta
                theta_plus = np.minimum ( theta_plus, self.max_vals )
                theta_minus = theta - ck*delta
                theta_minus = np.maximum ( theta_minus, self.min_vals )
                # Funcion values associated with ``theta_plus`` and 
                #``theta_minus``
                j_plus = self.calc_loss ( theta_plus )
                j_minus = self.calc_loss ( theta_minus )
                # Estimate the gradient
                ghat = ghat + ( j_plus - j_minus)/(2.*ck*delta)
            # Average gradient...
            ghat = ghat/float(ens_size)
            # The new parameter is the old parameter plus a scaled displacement
            # along the gradient.
            not_all_pass = True
            this_ak = ( theta*0 + 1 )*ak
            theta_new = theta
            while not_all_pass:
                out_of_bounds = np.where ( np.logical_or ( \
                    theta_new - this_ak*ghat > self.max_vals, 
                    theta_new - this_ak*ghat < self.min_vals ) )[0]
                theta_new = theta - this_ak*ghat
                if len ( out_of_bounds ) == 0:
                    theta = theta - this_ak*ghat
                    not_all_pass = False
                else:
                    this_ak[out_of_bounds] = this_ak[out_of_bounds]/2.
            
            # The new value of the gradient.
            j_new = self.calc_loss ( theta )
            j_sort = j_list.copy()
            j_sort.sort()
            # if j_new<j_sort[0]:
            #     theta3=theta2
            #     theta2=theta1
            #     theta1=theta
            #     j3 = j2
            #     j2 = j1
            #     j1 = j_new
            #     i3 = i2
            #     i2 = i1 
            #     i1 = n_iter
            j_list.append(j_new)
            
            # check fitting residual and adjust a
            # if n_iter % 1 == 0:
            if j_new>j_ori:
                self.a_par/=5
                ak = self.a_par/( n_iter + 1 + self.big_a_par)**self.alpha
                # theta = theta_0
                print('fitting worse a/=5, a={data}'.format(data=self.a_par))
            #     else:
            #         self.a_par*=1.2
            #         ak = self.a_par/( n_iter + 1 + self.big_a_par)**self.alpha
            #         j_ori = j_new
            #         # theta_0 = theta
            #         print('fitting better a*=1.2, a={data}'.format(data=self.a_par))

            # Be chatty to the user, tell him/her how it's going...
            if n_iter % report == 0:
                print ("\tIter %05d" % n_iter, np.sqrt(j_new * self.c_par**2), ak, ck)
                print ("\tTrained_parameter" , theta)
            # Functional tolerance: you can specify to ignore new theta values
            # that result in large shifts in the function value. Not a great
            # way to keep the results sane, though, as ak and ck decrease
            # slowly.
            if self.function_tolerance is not None:    
                if np.abs ( j_new - j_old ) > self.function_tolerance:
                    print ("\t No function tolerance!", np.abs ( j_new - j_old ))
                    theta = theta_saved
                    continue
                else:
                    j_old = j_new
            # You can also specify the maximum amount you want your parameters
            # to change in one iteration.
            if self.param_tolerance is not None:
                theta_dif = np.abs ( theta - theta_saved ) 
                if not np.all ( theta_dif < self.param_tolerance ):
                    print ("\t No param tolerance!", theta_dif < \
                        self.param_tolerance)
                    theta = theta_saved
                    continue
            # Ignore results that are outside the boundaries
            if (self.min_vals is not None) and (self.max_vals is not None):      
                i_max = np.where ( theta >= self.max_vals )[0]
                i_min = np.where ( theta <= self.min_vals )[0]
                if len( i_max ) > 0:
                    theta[i_max] = self.max_vals[i_max]*0.9
                if len ( i_min ) > 0:
                    theta[i_min] = self.min_vals[i_min]*1.1
            if report == 1:
                plt.plot ( theta, '-r' )
                plt.title ( "Iter %08d, J=%10.4G" % ( n_iter, j_new ))
                plt.grid ( True )
                plt.savefig ("SPSA_%08d.png" % n_iter, dpi=72 )
                plt.close()
            n_iter += 1
        # return ( theta, j_new, n_iter, j_list, theta1, theta2, theta3, j1, j2, j3, i1, i2, i3)
        return ( theta, j_new, n_iter, j_list)

def connect_db(db):
    """
    :param db: the data base name in text
    :return: a connection with database and a cursor
    """
    conn = sqlite3.connect(db)
    c = conn.cursor()
    return conn, c

def disconnect_db(conn):
    """
    :param conn: a sqlite database connection
    :return: None
    """
    conn.commit()   #apply changes to the database
    conn.close()

def  SPSA_match(Train_parameter, Input, Target, noise_var, a_par, min_vals, max_vals, max_iter, report, residual=1e-8):
    # Pump curve before training
    global pump_name
    # HP = ESP_fit_test(Input, Train_parameter)

    opti = SimpleSPSA ( ESP_fit_loss, a_par = a_par, args=(Input, Target, noise_var), \
        noise_var=noise_var, min_vals=min_vals, max_vals = max_vals,max_iter=max_iter, residual=residual )
    # ( Train_parameter, j_opt, niter, J_list, Train_parameter1, Train_parameter2, Train_parameter3, j1, j2, j3, i1, i2, i3 ) = opti.minimise (Train_parameter, report=report)
    ( Train_parameter, j_opt, niter, J_list) = opti.minimise (Train_parameter, report=report)
    print (Train_parameter, j_opt, niter)
    # print (Train_parameter1, j1, i1)
    # print (Train_parameter2, j2, i2)
    # print (Train_parameter3, j3, i3)

    fig, ax2 = plt.subplots(dpi=dpi, figsize = (3.33,2.5), nrows=1, ncols=1)

    for i in range(len(J_list)):
        J_list[i]=J_list[i]
    ax2.plot(J_list, linewidth=0.75)

    ax2.set_xlabel('Iteration', fontsize=8)
    ax2.set_ylabel('Loss', fontsize=8)

    if pump_name == 'Flex31':
        title='Training history MTESP'
    else:
        title='Training history '+pump_name
    ax2.set_title(title, fontsize=8)
    ax2.xaxis.set_tick_params(labelsize=8)
    ax2.yaxis.set_tick_params(labelsize=8)

    fig.tight_layout()
    # fig.savefig('/Users/haiwenzhu/Desktop/work/006 code/001 ESP/Python/ESP model py 2021/SPSA/'+title)
    fig.savefig('SPSA/'+str(title)+'.jpg')
    
    # return Train_parameter, Train_parameter1, Train_parameter2, Train_parameter3
    return Train_parameter

def ESP_parameter(Train_parameter):
    
    '''From base pump'''
    QBEM=QBEM_default[base_pump]
    ESP_input = ESP_default[base_pump]
    Coef_list_new = Coef_list.copy()
    


    if train_type == 'coefficient' or train_type == 'multi_CPU':
        '''ESP model coefficient'''
        ESP_output = ESP_default.copy()
        QBEM_out = QBEM_default.copy()
        Coef_list_new = Coef_list.copy()

        '''viscosity'''
        # # Coef_list_new['FTI']=Coef_list['FTI']*Train_parameter[0]
        # # Coef_list_new['FTD']=Coef_list['FTD']*Train_parameter[1]
        # Coef_list_new['FTI_coef']=Coef_list['FTI_coef']*Train_parameter[0]
        # Coef_list_new['FTD_coef']=Coef_list['FTD_coef']*Train_parameter[1]
        # Coef_list_new['SGMU_coef']=Coef_list['SGMU_coef']*Train_parameter[2]
        # Coef_list_new['QBEM_VISL_coef']=Coef_list['QBEM_VISL_coef']*Train_parameter[3]
        # Coef_list_new['FTI']=Coef_list['FTI']*Train_parameter[4]
        # Coef_list_new['FTD']=Coef_list['FTD']*Train_parameter[5]

        # QBEM = QBEM_default[pump_name]
        # QBEM = QBEM*Train_parameter[6]
        # QBEM_out=QBEM_default.copy()
        # QBEM_out [pump_name] = QBEM

        # Coef_list_new['SGMU_coef_2']=Coef_list['SGMU_coef_2']*Train_parameter[7]
        # Coef_list_new['SGM_coef']=Coef_list['SGM_coef']*Train_parameter[8]
        # Coef_list_new['SGM_coef_2']=Coef_list['SGM_coef_2']*Train_parameter[9]
        # # Coef_list_new['SGMU_coef_2']=1

        '''gas liquid'''
        QBEM = QBEM_default[pump_name]
        QBEM = QBEM*Train_parameter[0]
        QBEM_out=QBEM_default.copy()
        QBEM_out [pump_name] = QBEM
        Coef_list_new['FTI']=Coef_list['FTI']*Train_parameter[1]
        Coef_list_new['FTD']=Coef_list['FTD']*Train_parameter[2]
        Coef_list_new['SGMU_coef']=Coef_list['SGMU_coef']*Train_parameter[3]
        Coef_list_new['SGMU_coef_2']=Coef_list['SGMU_coef_2']*Train_parameter[4]

        # Coef_list_new['emulsion_E']=Coef_list['emulsion_E']*Train_parameter[4]
        # Coef_list_new['emulsion_SN']=Coef_list['emulsion_SN']*Train_parameter[5]
        # Coef_list_new['emulsion_WE']=Coef_list['emulsion_WE']*Train_parameter[6]
        # Coef_list_new['emulsion_RE']=Coef_list['emulsion_RE']*Train_parameter[7]
        # Coef_list_new['emulsion_ST1']=Coef_list['emulsion_ST1']*Train_parameter[8]
        # Coef_list_new['emulsion_ST2']=Coef_list['emulsion_ST2']*Train_parameter[9]
        Coef_list_new['alphaG_crit_coef']=Coef_list['alphaG_crit_coef']*Train_parameter[10]
        Coef_list_new['alphaG_crit_critical']=Coef_list['alphaG_crit_critical']*Train_parameter[11]
        Coef_list_new['alphaG_crit_critical2']=Coef_list['alphaG_crit_critical2']*Train_parameter[12]
        Coef_list_new['factor']=Coef_list['factor']*Train_parameter[13]
        Coef_list_new['DB_GVF_EFF']=Coef_list['DB_GVF_EFF']*Train_parameter[14]
        Coef_list_new['DB_NS_EFF']=Coef_list['DB_NS_EFF']*Train_parameter[15]
        Coef_list_new['CD_Gas_EFF']=Coef_list['CD_Gas_EFF']*Train_parameter[16]
        Coef_list_new['CD_Liquid_EFF']=Coef_list['CD_Liquid_EFF']*Train_parameter[17]
        Coef_list_new['CD_INT_EFF']=Coef_list['CD_INT_EFF']*Train_parameter[18]
        Coef_list_new['CD_GV_EFF']=Coef_list['CD_GV_EFF']*Train_parameter[19]
    elif train_type == 'water_QBEM':
        
        QBEM = QBEM_default[pump_name]
        QBEM = QBEM*Train_parameter[0]
        QBEM_out=QBEM_default.copy()
        QBEM_out [pump_name] = QBEM

        ESP_output = ESP_default.copy()
        Coef_list_new = Coef_list.copy()
        Coef_list_new['SGM_coef_2']=Coef_list['SGM_coef_2']*Train_parameter[1]

    elif train_type == '4_geometry':
        '''Train ESP 4 key coefficient for geometry'''
        # trained parameter
        coeff_D=Train_parameter[0]
        coeff_L= Train_parameter[1]
        coeff_B = Train_parameter[2]
        QBEM = QBEM*Train_parameter[3]
        coeff_ZI = 1
        coeff_ZD = 1

        ESP_input = ESP_default[base_pump]
        ESP = ESP_input.copy()
        ESP['R1'] = ESP_input['R1']*coeff_D
        ESP['R2'] = ESP_input['R2']*coeff_D
        ESP['TB'] = ESP_input['TB']*coeff_D
        ESP['TV'] = ESP_input['TV']*coeff_D
        ESP['RD1'] = ESP_input['RD1']*coeff_D
        ESP['RD2'] = ESP_input['RD2']*coeff_D
        ESP['B2'] = ESP_input['B2']*coeff_B
        ESP['ZI'] = ESP_input['ZI']*coeff_ZI
        ESP['ZD'] = ESP_input['ZD']*coeff_ZD
        # ESP['ZI'] = ESP_input['ZI']*1
        # ESP['ZD'] = ESP_input['ZD']*1
        ESP['YI1'] = ESP_input['YI1']*coeff_D
        ESP['YI2'] = ESP_input['YI2']*coeff_D
        ESP['LI'] = ESP_input['LI']*coeff_D
        ESP['LD'] = ESP_input['LD']*coeff_D
        ESP['VOI'] = ESP_input['VOI']*coeff_L*coeff_D
        ESP['VOD'] = ESP_input['VOD']*coeff_L*coeff_D
        ESP['ASF'] = ESP_input['ASF']*coeff_D
        ESP['ASB'] = ESP_input['ASB']*coeff_D
        ESP['AB'] = ESP_input['AB']*coeff_D
        ESP['AV'] = ESP_input['AV']*coeff_D
        ESP['ADF'] = ESP_input['ADF']*coeff_D
        ESP['ADB'] = ESP_input['ADB']*coeff_D
        ESP['RLK'] = ESP_input['RLK']*coeff_D
        ESP['LG'] = ESP_input['LG']*coeff_D
        ESP['SL'] = ESP_input['SL']*coeff_D
        ESP['B1'] = ESP_input['B1']*coeff_B
        ESP['NS'] = ESP_input['NS']
        ESP_output = ESP_default.copy()
        ESP_output [pump_name] = ESP
        QBEM_out=QBEM_default.copy()
        QBEM_out [pump_name] = QBEM


    elif train_type == 'ALL_geometry':
        '''Train ESP all geometry'''
        ESP = ESP_input.copy()
        QBEM = QBEM*Train_parameter[0]
        ESP['R1'] = ESP_input['R1']*Train_parameter[1]
        ESP['R2'] = ESP_input['R2']*Train_parameter[2]
        ESP['TB'] = ESP_input['TB']*Train_parameter[3]
        ESP['TV'] = ESP_input['TV']*Train_parameter[4]
        ESP['RD1'] = ESP_input['RD1']*Train_parameter[5]
        ESP['RD2'] = ESP_input['RD2']*Train_parameter[6]
        ESP['B2'] = ESP_input['B2']*Train_parameter[7]
        ESP['ZI'] = ESP_input['ZI']*Train_parameter[8]
        ESP['ZD'] = ESP_input['ZD']*Train_parameter[9]
        ESP['YI1'] = ESP_input['YI1']*Train_parameter[10]
        ESP['YI2'] = ESP_input['YI2']*Train_parameter[11]
        ESP['LI'] = ESP_input['LI']*Train_parameter[12]
        ESP['LD'] = ESP_input['LD']*Train_parameter[13]
        ESP['VOI'] = ESP_input['VOI']*Train_parameter[14]
        ESP['VOD'] = ESP_input['VOD']*Train_parameter[15]
        ESP['ASF'] = ESP_input['ASF']*Train_parameter[16]
        ESP['ASB'] = ESP_input['ASB']*Train_parameter[17]
        ESP['AB'] = ESP_input['AB']*Train_parameter[18]
        ESP['AV'] = ESP_input['AV']*Train_parameter[19]
        ESP['ADF'] = ESP_input['ADF']*Train_parameter[20]
        ESP['ADB'] = ESP_input['ADB']*Train_parameter[21]
        ESP['RLK'] = ESP_input['RLK']*Train_parameter[22]
        ESP['LG'] = ESP_input['LG']*Train_parameter[23]
        ESP['SL'] = ESP_input['SL']*Train_parameter[24]
        ESP['B1'] = ESP_input['B1']*Train_parameter[25]
        ESP['NS'] = ESP_input['NS']
        ESP_output = ESP_default.copy()
        ESP_output [pump_name] = ESP
        QBEM_out=QBEM_default.copy()
        QBEM_out [pump_name] = QBEM


    return ESP_output, QBEM_out, Coef_list_new

def ESP_fit_test(Input, Train_parameter):

    

    ESP_ALL, QBEM_ALL, Coef_list_new = ESP_parameter(Train_parameter=Train_parameter)

    # HP = []
    loss = 0
    for Pump in Input.Pump.unique():
        input_i = Input[Input['Pump']==Pump]
        _, HP_i = ESP_validation(ESP_ALL[Pump],QBEM_ALL[Pump],input_i,Pump,False,Coef_list_new).ESP_curve(QL=input_i["QL_bpd"], 
            QG=input_i['QG_bpd'],
            VISO_in = input_i["TargetVISL_cp"]/1000, 
            DENO_in = input_i["DENL_kgm3"], 
            DENG_std=1.225,
            WC=input_i["TargetWC_%"], 
            RPM= input_i["RPM"],
            VISG_std=0.000018,
            O_W_ST=0.035,
            G_L_ST=0.073,
            P=input_i["Ptank_psi"], 
            T=input_i["Tin_F"], 
            curve_type = curve_type)
        # HP_i = HP_i.tolist()
        # HP += HP_i
        loss += np.sum((input_i['DP_psi']-HP_i)**2/noise_var**2)

    return loss

def ESP_fit_loss(Train_parameter, Input, Target, noise_var):

    ESP_ALL, QBEM_ALL, Coef_list_new = ESP_parameter(Train_parameter=Train_parameter)

    # e1_list = []
    # e2_list = []
    # e3_list = []
    # for Pump in Input.Pump.unique():
    #     input_i = Input[Input['Pump']==Pump]
        
    #     e1, e2, e3, e4, e5, e6 = ESP_validation(ESP_ALL[Pump],QBEM_ALL[Pump],input_i,Pump,False,Coef_list_new).error_analysis()
    #     e1_list.append(e4/noise_var**2)
    #     e2_list.append(e5/noise_var**2)
    #     e3_list.append(e6/noise_var**2)

    # loss = sum(e1_list)
    # loss = sum(e2_list)
    # loss = sum(e3_list)

    loss = 0
    for Pump in Input.Pump.unique():
        input_i = Input[Input['Pump']==Pump]
        _, HP_i = ESP_validation(ESP_ALL[Pump],QBEM_ALL[Pump],input_i,Pump,False,Coef_list_new).ESP_curve(QL=input_i["QL_bpd"], 
            QG=input_i['QG_bpd'],
            VISO_in = input_i["TargetVISL_cp"]/1000, 
            DENO_in = input_i["DENL_kgm3"], 
            DENG_std=1.225,
            WC=input_i["TargetWC_%"], 
            RPM= input_i["RPM"],
            VISG_std=0.000018,
            O_W_ST=0.035,
            G_L_ST=0.073,
            P=input_i["Ptank_psi"], 
            T=input_i["Tin_F"], 
            curve_type = curve_type)
        loss += np.sum((input_i['DP_psi']-HP_i)**2/noise_var**2)
    # loss = np.sum(((ESP_fit_test(Input, Train_parameter)-Target))**2/noise_var**2)
    return loss

def test():
    
    '''water QBEM'''
    conn, c = connect_db('ESP.db')
    Exp_data = pd.read_sql_query("SELECT * FROM Catalog_All;", conn)
    Exp_data = Exp_data[Exp_data.QL_bpd != 0]
    Exp_data = Exp_data.reset_index(drop=True)
    disconnect_db(conn)
    Input = Exp_data
    Target = Exp_data.DP_psi


    train_type = 'water_QBEM'
    noise_var=0.1
    a_par=1e-4
    max_iter=5
    report=10

    Train_parameter = np.ones(2)
    min_vals = [0.991, 0.1]
    max_vals = [1.001, 10]
    qbem = []

    _, Train_parameter1, _, _ = SPSA_match(Train_parameter, Input, Input.DP_psi, noise_var, a_par, min_vals, max_vals, max_iter, report)


    for pump in Input.Pump.unique():
        df = Input[Input.Pump == pump]
        base_pump = pump
        pump_name = pump
        _, Train_parameter1, _, _ = SPSA_match(Train_parameter, df, df.DP_psi, noise_var, a_par, min_vals, max_vals, max_iter, report)
        ESP, QBEM, Coef = ESP_parameter(Train_parameter1)
        fig1, ax = plt.subplots(dpi = dpi, figsize = (4,3))
        ESP_validation(ESP_GEO=ESP[pump],QBEM=QBEM[pump],Exp_data=df,pump_name=pump,bx=ax).water_validation()
        print(pump, QBEM_default[pump]*Train_parameter1[0], Coef_list['SGM_coef']*Train_parameter1[1])
        qbem.append([pump, QBEM_default[pump]*Train_parameter1[0], Coef_list['SGM_coef']*Train_parameter1[1]])
    
    print(qbem)
    plt.show()

if __name__ == "__main__":

    Coef_list = Coef_list

    
    noise_var=0.1
    a_par=1e-6
    max_iter=200
    report=100

    Train_parameter = np.ones(26)*(1)   # all geometry
    min_vals=np.ones(26)*(0.8)
    max_vals = np.ones(26)*(1.2)

    train_type = 'coefficient'
    noise_var=0.1
    a_par=1e-7
    max_iter=500
    report=10

    Train_parameter = np.ones(10)*(1)

    # FTI_coef, FTD_coef, SGMU_coef, QBEM_VISL_coef, FTI, FTD, QBEM, SGMU_coef_2, SGM_coef, SGM_coef——2
    min_vals = [0.01, 0.01, 0.001, 0.01, 0.1, 0.1, 0.5, 0.001, 0.001, 0.001]
    max_vals = [10,   10,   1000,  1000, 10,   10, 2,   1000,  1000,  1000]


    # Train_parameter = np.ones(3)*(1)   # 4 geometry
    # min_vals=np.ones(3)*(0.2)
    # max_vals = np.ones(3)*(5)


    # Train_parameter = np.ones(19)*(1)   # coefficient
    # min_vals=np.ones(19)*(0.2)
    # max_vals = np.ones(19)*(5)


    '''GC6100 gas liquid'''
    # match
    base_pump = 'GC6100'
    pump_name = 'GC6100'
    output_file = 'SPSA/gas liquid/'
    TargetRPM = 3000
    TargetP_psi = 150
    curve_type = 'GL'

    conn, c = connect_db('ESP.db')
    Exp_data = pd.read_sql_query("SELECT * FROM All_pump "
                              + "ORDER BY TargetQG_bpd, QL_bpd"
                              + ";", conn)
    # Exp_data = Exp_data[(Exp_data.Pump == pump_name) & (Exp_data.Test == 'Mapping') & 
    #             (Exp_data.TargetRPM == TargetRPM) & (Exp_data.TargetP_psi == TargetP_psi)]
    Exp_data = Exp_data[(Exp_data.Pump == pump_name) & (Exp_data.Test == 'Mapping') & 
                 (Exp_data.TargetP_psi == TargetP_psi) & (Exp_data.TargetQG_bpd < 6)]
    Exp_data = Exp_data[(Exp_data.DP_psi > 0)]
    Exp_data=Exp_data.reset_index(drop=True)
    disconnect_db(conn)


    Input = Exp_data
    Target = Exp_data.DP_psi
    print('Train data: ', Input.shape[0])
    
    train_type = 'coefficient'
    noise_var=0.1
    a_par=1e-9
    max_iter=200
    report=10

    Train_parameter = np.ones(20)*(1)
    min_vals = np.ones(20)*(0.5)
    max_vals = np.ones(20)*(10)
    min_vals [0] = 0.9
    max_vals [0] = 1.2

    Train_parameter1= SPSA_match(Train_parameter, Input, Target, noise_var, a_par, min_vals, max_vals, max_iter, report)
    # Train_parameter1 = Train_parameter

    error = []
    ESP_ALL_new,QBEM_ALL_new,Coef_list_new = ESP_parameter(Train_parameter1)
    
    for Pump in Exp_data['Pump'].unique():
        df_1 = Exp_data[Exp_data['Pump']==Pump]
        for RPM in df_1['RPM'].unique():
            df_2 = df_1[df_1['RPM']==RPM]
            fig1, (ax, ax2) = plt.subplots(dpi = dpi, figsize = (8,3), nrows=1, ncols=2)
            fig2, (bx, bx2) = plt.subplots(dpi = dpi, figsize = (8,3), nrows=1, ncols=2)
            # ori pump
            
            ESP_validation(ESP_default[Pump],QBEM_default[Pump],df_2,Pump,ax,Coef=Coef_list).GL_validation()
            ax.set_title('Before train (GL): '+Pump+' at '+str(int(df_2['RPM'].mean()))+ ' RPM', fontsize=8)
            e1_ori, e2_ori, e3_ori, e4_ori, e5_ori, e6_ori = ESP_validation(ESP_default[Pump],QBEM_default[Pump],df_2,Pump,bx,Coef=Coef_list).error_analysis()
            bx.set_title(r'Before train error analysis: e1: %.d, e2: %.d, e3: %.d' % (e1_ori, e2_ori, e3_ori), fontsize=8)
            
            # new pump
            ESP_validation(ESP_ALL_new[Pump],QBEM_ALL_new[Pump],df_2,Pump,ax2,Coef=Coef_list_new).GL_validation()
            ax2.set_title('After train (Oil): '+Pump+' at '+str(int(df_2['RPM'].mean()))+ ' RPM', fontsize=8)
            e1, e2, e3, e4, e5, e6 = ESP_validation(ESP_ALL_new[Pump],QBEM_ALL_new[Pump],df_2,Pump,bx2,Coef=Coef_list_new).error_analysis()
            bx2.set_title(r'After train error analysis: e1: %.d, e2: %.d, e3: %.d' % (e1, e2, e3), fontsize=8)
            error.append([Pump, RPM, 'new', e1, e2, e3, e4, e5, e6, 'ori', e1_ori, e2_ori, e3_ori, e4_ori, e5_ori, e6_ori ])
            fig1.tight_layout()
            fig2.tight_layout()
            fig1.savefig(output_file+'Gas liquid performance for '+str(Pump)+' at ' +str(RPM) + ' RPM.jpg')
            fig2.savefig(output_file+'Gas liquid error for '+str(Pump)+' at ' +str(RPM) + ' RPM.jpg')
    
    error = pd.DataFrame(error)
    error.to_excel(output_file+'error.xls')

    Train_parameter1 = pd.DataFrame(Train_parameter1)
    Train_parameter1.to_excel(output_file+'Train_parameter_coefficient.xls')

    
    
    print(error)
    

    plt.show()

    # match
    ''' viscosity '''
    # base_pump = 'DN1750'
    # pump_name = 'DN1750'
    # # TargetRPM = 3500
    # # QBEM_default[pump_name] = 11000
    # output_file = 'SPSA/oil curve/'

    # train_type = 'coefficient'
    # noise_var=0.1
    # a_par=1e-7
    # max_iter=2000
    # report=1000
    # residuel=1e-8

    # Train_parameter = np.ones(10)*(1)

    # # FTI_coef, FTD_coef, SGMU_coef, QBEM_VISL_coef, FTI, FTD, QBEM, SGMU_coef_2, SGM_coef, SGM_coef——2
    # min_vals = [0.01, 0.01, 0.001, 0.01, 0.1, 0.1, 0.2, 0.001, 0.001, 0.001]
    # max_vals = [10,   10,   1000,  1000, 10,   10, 5,   1000,  1000,  1000]

    # conn, c = connect_db('ESP.db')
    # Exp_data = pd.read_sql_query("SELECT * "
    #                           + "FROM df_Viscosity "
    #                         #   + "WHERE [TargetWC_%] = 0 "
    #                           + "ORDER BY Pump, RPM, TargetVISL_cp, QL_bpd"
    #                           + ";", conn)
    # disconnect_db(conn)
    
    # Exp_data.drop(Exp_data[Exp_data['Case']=='DN1750_Solano'].index, inplace=True)
    # Exp_data.drop(Exp_data[Exp_data['Case']=='DN1750_Banjar'].index, inplace=True)
    # Exp_data.drop(Exp_data[Exp_data['Case']=='DN1750_Solano_Ave'].index, inplace=True)
    # # Exp_data.drop(Exp_data[Exp_data['Case']=='P100'].index, inplace=True)
    # # Exp_data.drop(Exp_data[Exp_data['Case']=='DN1750_CFD'].index, inplace=True)
    # Exp_data.reset_index(drop=True, inplace=True)
    # # Exp_data = Exp_data[(Exp_data.Pump.str.startswith(pump_name)) & (Exp_data.RPM == TargetRPM)]
    # # Exp_data = Exp_data[Exp_data['Case']=='P100']
    # # Exp_data = Exp_data[(Exp_data.Pump.str.startswith(pump_name))]
    # Exp_data = Exp_data[Exp_data.QL_bpd != 0]
    # Exp_data = Exp_data[Exp_data.TargetVISL_cp <  700]
    # Exp_data = Exp_data.reset_index(drop=True)

    # # Test_data = Exp_data[(Exp_data.Pump.str.startswith(pump_name)) & (Exp_data.RPM == TargetRPM)]
    
    # # Exp_data = Exp_data[(Exp_data.Pump.str.startswith(pump_name)) & (Exp_data.RPM == TargetRPM) & ((Exp_data.TargetVISL_cp == 1000) | (Exp_data.TargetVISL_cp == 50))]
    # Exp_data = Exp_data.reset_index(drop=True)

    # Input = Exp_data
    # Target = Exp_data.DP_psi
    # print('Train data: ', Input.shape[0])


    # # fig1, (ax1, ax2) = plt.subplots(dpi = dpi, figsize = (6.66,2.5), nrows=1, ncols=2)
    # # fig2, (bx1, bx2) = plt.subplots(dpi = dpi, figsize = (6.66,2.5), nrows=1, ncols=2)
    # # ESP_ALL,QBEM_ALL,Coef_list = ESP_parameter(Train_parameter)
    # # update_global(ESP_ALL=ESP_ALL,QBEM_ALL=QBEM_ALL, Coef=Coef_list)
    # # ax1 = ESP_validation(pump_name, Test_data, ax1).error_analysis()
    # # bx1 = ESP_validation(pump_name, Test_data, bx1).Oil_validation()

    # time0 = time.time()
    # # time0 = datetime.datetime.now()
    # Train_parameter1= SPSA_match(Train_parameter, Input, Target, noise_var, 
    #                 a_par, min_vals, max_vals, max_iter, report, residuel)

    # # ESP_ALL,QBEM_ALL,Coef_list = ESP_parameter(Train_parameter1)
    # # update_global(ESP_ALL=ESP_ALL,QBEM_ALL=QBEM_ALL, Coef=Coef_list)
    # # ax2 = ESP_validation(pump_name, Test_data, ax2).error_analysis()
    # # bx2 = ESP_validation(pump_name, Test_data, bx2).Oil_validation()

    # time1 = time.time()

    # error = []
    # ESP_ALL_new,QBEM_ALL_new,Coef_list_new = ESP_parameter(Train_parameter1)
    # for Pump in Exp_data['Pump'].unique():
    #     df_1 = Exp_data[Exp_data['Pump']==Pump]
    #     for RPM in df_1['RPM'].unique():
    #         df_2 = df_1[df_1['RPM']==RPM]
    #         fig1, (ax, ax2) = plt.subplots(dpi = dpi, figsize = (8,3), nrows=1, ncols=2)
    #         fig2, (bx, bx2) = plt.subplots(dpi = dpi, figsize = (8,3), nrows=1, ncols=2)
    #         # ori pump
    #         ESP_validation(ESP_default[Pump],QBEM_default[Pump],df_2,Pump,ax,Coef=Coef_list).Oil_validation()
    #         ax.set_title('Before train (Oil): '+Pump+' at '+str(int(df_2['RPM'].mean()))+ ' RPM', fontsize=8)
    #         e1_ori, e2_ori, e3_ori, e4_ori, e5_ori, e6_ori = ESP_validation(ESP_default[Pump],QBEM_default[Pump],df_2,Pump,bx,Coef=Coef_list).error_analysis()
    #         bx.set_title(r'Before train error analysis: e1: %.d, e2: %.d, e3: %.d' % (e1_ori, e2_ori, e3_ori), fontsize=8)
    #         # new pump
    #         ESP_validation(ESP_ALL_new[Pump],QBEM_ALL_new[Pump],df_2,Pump,ax2,Coef=Coef_list_new).Oil_validation()
    #         ax2.set_title('After train (Oil): '+Pump+' at '+str(int(df_2['RPM'].mean()))+ ' RPM', fontsize=8)
    #         e1, e2, e3, e4, e5, e6 = ESP_validation(ESP_ALL_new[Pump],QBEM_ALL_new[Pump],df_2,Pump,bx2,Coef=Coef_list_new).error_analysis()
    #         bx2.set_title(r'After train error analysis: e1: %.d, e2: %.d, e3: %.d' % (e1, e2, e3), fontsize=8)
    #         error.append([Pump, RPM, 'new', e1, e2, e3, e4, e5, e6, 'ori', e1_ori, e2_ori, e3_ori, e4_ori, e5_ori, e6_ori ])
    #         fig1.tight_layout()
    #         fig2.tight_layout()
    #         fig1.savefig('SPSA/oil curve/Oil performance for '+str(Pump)+' at ' +str(RPM) + ' RPM')
    #         fig2.savefig('SPSA/oil curve/Oil error for '+str(Pump)+' at ' +str(RPM) + ' RPM')

    # error = pd.DataFrame(error)
    # Train_parameter1 = pd.DataFrame(Train_parameter1)
    # print(error)
    
    # error.to_excel(output_file+'error.xls')
    # Train_parameter1.to_excel(output_file+'Train_parameter_coefficient.xls')


    # print ('simulation time: ', (time1-time0))

    plt.show()
