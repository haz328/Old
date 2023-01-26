# -*- coding: utf-8 -*-

"""
Typical inputs: Q in bpd, N in rpm, angle in degree, all others in SI units
"""
#SN: stage number
#   AB		=	impeller blade surface area on both side (m2)				
#   AV		=	diffuser vane surface area on both side (m2)			
#   ASF	    =	shroud front surface area (m2)				
#   ASB	    =	shroud back surface area (m2)				
#   ADF	    =	diffuser front surface area (m2)					
#   ADB	    =	diffuser back surface area (m2)	
#   RD1     =   Diffuser inlet radius, input for sgl_calculate_jc by Jiecheng Zhang
#   RD2     =   Diffuser outlet radius, input for sgl_calculate_jc by Jiecheng Zhang
#   NS      =   specific speed based on field units
#   SL      =   Leakage gap width
#   LG      =   Leakage gap length
QBEM_default = {'TE2700': 4500, 'DN1750': 3000, 'GC6100': 7800, 'P100': 11000, 'Flex31': 5000}    

ESP = {'TE2700':
            {
                "R1": 0.017496,     "R2": 0.056054,     "TB": 0.00272,      "TV": 0.00448,      "RD1": 0.056,
                "RD2": 0.017496,    "YI1": 0.012194,    "YI2": 0.007835,    "VOI": 0.000016119, "VOD": 0.000011153,
                "ASF": 0.00176464,  "ASB": 0.00157452,  "AB": 0.001319,     "AV": 0.001516,     "ADF": 0.001482,
                "ADB": 0.000935,    "LI": 0.076,        "LD": 0.08708,      "RLK": 0.056209,    "LG": 0.00806,
                "SL": 0.00005,      "EA": 0.000254,     "ZI": 5,            "ZD": 9,            "B1": 19.5,
                "B2": 24.7,         "NS": 1600,         "DENL": 1000,       "DENG": 11.2,       "DENW": 1000,       "VISL": 0.001,
                "VISG": 0.000018,   "VISW": 0.001,      "ST": 0.073,        "N": 3500,          "SGM": 0.3,
                "QL": 2700,         "QG": 50,           "GVF": 10,          "WC": 0.,           "SN": 1
            },
       'GC6100':
           {
                "R1": 0.027746,     "R2": 0.050013,     "TB": 0.0019862,    "TV": 0.002894,     "RD1": 0.0547,
                "RD2": 0.017517,    "YI1": 0.017399,    "YI2": 0.013716,    "VOI": 1.512E-5,    "VOD": 1.9818E-5,
                "ASF": 8.9654E-4,   "ASB": 9.4143E-4,   "AB": 1.0333E-3,    "AV": 1.769E-3,     "ADF": 2.0486E-3,
                "ADB": 1.0301E-3,   "LI": 0.0529,       "LD": 0.0839,       "RLK": 6.1237E-2,   "LG": 0.0015475,
                "SL": 3.81E-4,      "EA": 0.0003,       "ZI": 7,            "ZD": 8,            "B1": 33.375,
                "B2": 41.387,       "NS": 3220,         "DENL": 1000,       "DENG": 11.2,       "DENW": 1000,       "VISL": 0.001,
                "VISG": 0.000018,   "VISW": 0.001,      "ST": 0.073,        "N": 2400,          "SGM": 0.3,
                "QL": 6100,         "QG": 50,           "GVF": 10,          "WC": 0.,           "SN": 1
           },
        # original R1 = 0.014351   TB = 0.0025896

        'Flex31':
            {
                "R1": 0.018,        "R2": 0.039403,     "TB": 0.0018875,    "TV": 0.0030065,    "RD1": 0.042062,
                "RD2": 0.018841,    "YI1": 0.015341,    "YI2": 0.01046,     "VOI": 0.000010506, "VOD": 0.000010108,
                "ASF": 0,           "ASB": 0,           "AB": 0.0020663,    "AV": 0.0025879,    "ADF": 0, 
                "ADB": 0,           "LI": 0.04252,      "LD": 0.060315,     "RLK": 0.033179,    "LG": 0.005, 
                "SL": 0.000254,     "EA": 0.0003,       "ZI": 6,            "ZD": 8,            "B1": 21.99,        
                "B2": 57.03,        "NS": 2975,         "DENL": 1000,       "DENG": 11.2,       "DENW": 1000,       "VISL": 0.001,      
                "VISG": 0.000018,   "VISW": 0.001,      "ST": 0.073,        "N": 3600,          "SGM": 0.3,
                "QL": 4000,         "QG": 50,           "GVF": 10,          "WC": 0.,           "SN": 1
            },


        'DN1750':
           {
                "R1": 1.9875E-2,    "R2": 3.5599E-2,    "TB": 1.7E-3,       "TV": 3.12E-3,      "RD1": 0.04,
                "RD2": 0.01674,     "YI1": 1.3536E-2,   "YI2": 7.13E-3,     "VOI": 6.283E-6,    "VOD": 7.063E-6,
                "ASF": 6.8159E-04,  "ASB": 6.549E-04,   "AB": 6.9356E-04,   "AV": 7.1277E-04,   "ADF": 1.0605E-03,
                "ADB": 6.6436E-04,  "LI": 0.039,        "LD": 5.185E-02,    "RLK": 0.04,        "LG": 0.01,
                "SL": 0.00005,      "EA": 0.000254,     "ZI": 6,            "ZD": 8,            "B1": 20.3,
                "B2": 36.2,         "NS": 2815,         "DENL": 1000,       "DENG": 11.2,       "VISL": 0.001,
                "VISG": 0.000018,   "VISW": 0.001,      "ST": 0.073,        "N": 3500,          "SGM": 0.3,
                "QL": 1750,         "QG": 50,           "GVF": 10,          "WC": 0.,           "SN": 3
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
                "QL": 9000,         "QG": 50,           "GVF": 10,          "WC": 0.,           "SN": 1
           }
           
       }