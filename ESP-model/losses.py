# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from ESP_Class import SinglePhaseModel

"""
multiple classes that describe different existing loss models in literature
"""

# global variables
G = 9.81
pi = np.pi
E1 = 1e-5
DENW = 997.
psi_to_pa = 1.013e5 / 14.7
bbl_to_m3 = 0.15897


##################################
# Friction loss calculation class
class FrictionLoss(object):
    def __init__(self, inputs, N, DENL, VISL):
        self.R1 = inputs['R1']
        self.R2 = inputs['R2']
        self.RD1 = inputs['RD1']
        self.RD2 = inputs['RD2']
        self.TB = inputs['TB']
        self.TV = inputs['TV']
        self.YI1 = inputs['YI1']
        self.YI2 = inputs['YI2']
        self.VOI = inputs['VOI']
        self.VOD = inputs['VOD']
        self.ASF = inputs['ASF']
        self.ASB = inputs['ASB']
        self.AB = inputs['AB']
        self.AV = inputs['AV']
        self.ADF = inputs['ADF']
        self.ADB = inputs['ADB']
        self.LI = inputs['LI']
        self.LD = inputs['LD']
        self.RLK = inputs['RLK']
        self.LG = inputs['LG']
        self.SL = inputs['SL']
        self.EA = inputs['EA']
        self.ZI = inputs['ZI']
        self.ZD = inputs['ZD']
        self.B1 = inputs['B1'] * (pi / 180.0)
        self.B2 = inputs['B2'] * (pi / 180.0)
        self.DENL = DENL
        self.VISL = VISL
        self.VISW = inputs['VISW']
        self.N = N
        self.NS = inputs['NS']
        self.OMEGA = 2.0 * pi * self.N / 60.0
        a1 = 2 * pi * self.R1 / self.ZI * np.sin(self.B1)
        b1 = self.YI1
        a2 = 2 * pi * self.R2 / self.ZI * np.sin(self.B2)
        b2 = self.YI2
        self.DH1 = 2 * a1 * b1 / (a1 + b1)
        self.DH2 = 2 * a2 * b2 / (a2 + b2)
        self.DH = (self.DH1 + self.DH2) / 2.0
        self.beta = (self.B1 + self.B2) / 2.0
        self.H = (self.YI1 + self.YI2) / 2.0
        self.R = (self.R1 + self.R2) / 2.0

    # Sun(2003) model consider both impeller and diffuser friction
    # In this function, only consider impeller, which is the dominant part
    def sun2003(self, Q):
        fff = self.Churchill_1977(Q)
        Fshape = self.sun2003_shape_effect(Q)
        Fcurve = self.sun2003_curvature_effect(Q)
        Frotation = self.sun2003_rotation_effect(Q)
        friction_factor = fff * Fshape * Fcurve * Frotation
        bm = (self.YI1 + self.YI2) / 2.
        h_friction = friction_factor * Q ** 2 * (self.R2 - self.R1) / (8 * 9.81 * self.DH * np.pi ** 2 * bm ** 2
                                                                       * (np.sin(self.beta)) ** 3 * self.R1 * self.R2)
        return h_friction

    def Churchill_1977(self, Q):
        NRe = self.DH * Q * self.DENL / (2 * np.pi * self.R * self.H * self.VISL * np.sin(self.beta))

        if NRe < 2300:
            return 64.0 / NRe
        else:
            return 8 * (2.457 * np.log(1. / ((7. / NRe) ** 0.9 + 0.27 * (self.EA / self.DH)))) ** (-2)

    def sun2003_shape_effect(self, Q):
        # a channel width, b channel height
        a1 = 2 * np.pi * self.R1 / self.ZI
        b1 = self.YI1
        a2 = 2 * np.pi * self.R2 / self.ZI
        b2 = self.YI2
        ll = (np.min([a1, b1]) / np.max([a1, b1]) + np.min([a2, b2]) / np.max([a2, b2])) / 2.0
        NRe = self.DH * Q * self.DENL / (2 * np.pi * self.R * self.H * self.VISL * np.sin(self.beta))

        if NRe < 2300:
            return 1 / (2. / 3. + 11. / 24. * ll * (2 - ll))
        else:
            return 1 / (2. / 3. + 11. / 24. * ll * (2 - ll)) ** 0.25

    def sun2003_curvature_effect(self, Q):
        # radius of curvature of radial type centrifugal pump used in Sun (2003) dissertation
        RC = 1. / np.sin(self.beta) * (1. / np.tan(self.beta) / self.R -
                                       (self.B1 - self.B2) / (self.R1 - self.R2)) ** (-1)
        rH = self.DH / 2.
        NRe = self.DH * Q * self.DENL / (2 * np.pi * self.R * self.H * self.VISL * np.sin(self.beta))

        if RC / rH >= 860:
            NRec = 2300.
        else:
            NRec = 2e4 * (rH / RC) ** 0.32

        if NRe < NRec:
            if RC / rH >= 860:
                return 1.
            else:
                return 0.266 * NRe ** 0.389 * (rH / RC) ** 0.1945
        else:
            if NRe * (rH / RC) ** 2 >= 300:
                return (NRe * (rH / RC) ** 2) ** 0.05
            elif (NRe * (rH / RC) ** 2 > 0.0034) and (NRe * (rH / RC) ** 2 < 300):
                return 0.092 * (NRe * (rH / RC) ** 2) ** 0.25 + 0.962
            else:
                return 1.

    def sun2003_rotation_effect(self, Q):
        NRe_omega = self.OMEGA * self.DH ** 2 * self.DENL / self.VISL
        NRe = self.DH * Q * self.DENL / (2 * np.pi * self.R * self.H * self.VISL * np.sin(self.beta))

        if NRe_omega >= 28:
            NRe_omegac = 1070 * NRe_omega ** 0.23
        else:
            NRe_omegac = 2300

        if NRe < NRe_omegac:
            KLaminar = NRe * NRe_omega
            if (KLaminar <= 220) and (NRe_omega / NRe < 0.5):
                return 1.
            elif (KLaminar > 220) and (NRe_omega / NRe < 0.5) and (KLaminar < 1e7):
                fr = 0.0883 * KLaminar ** 0.25 * (1. + 11.2 * KLaminar ** (-0.325))
                if fr > 1:
                    return 1.
                else:
                    return fr
            else:
                fr = (0.0672 * NRe_omega ** 0.5) / (1. - 2.11 * NRe_omega ** (-0.5))
                if fr > 1:
                    return 1.
                else:
                    return fr
        else:
            KTurbulent = NRe_omega ** 2 / NRe
            if KTurbulent < 1.:
                return 1.
            elif (KTurbulent) >= 1. and (KTurbulent < 15.):
                return 0.942 + 0.058 * KTurbulent ** 0.282
            else:
                return 0.942 * KTurbulent ** 0.05

    # Thin et al. (2008) friction loss model
    def Thin2008(self, Q):
        a2 = 2 * np.pi * self.R2 / self.ZI
        b2 = self.YI2
        rH = 2 * a2 * b2 / (a2 + b2)
        C1M = Q / ((2.0 * np.pi * self.R1 - self.ZI * self.TB) * self.YI1)
        C2M = Q / ((2.0 * np.pi * self.R2 - self.ZI * self.TB) * self.YI2)
        W1 = C1M / np.sin(self.B1)
        W2 = C2M / np.sin(self.B2)
        h_friction = b2 * (2 * self.R2 - 2 * self.R1) * (W1 + W2) ** 2 / (8 * 9.81 * np.sin(self.B2) * rH)
        return h_friction

    # Bing et al. (2012) friction model
    def Bing2012(self, Q):
        fff = self.Churchill_1977(Q)
        Fshape = self.sun2003_shape_effect(Q)
        Fcurve = self.sun2003_curvature_effect(Q)
        Frotation = self.sun2003_rotation_effect(Q)
        friction_factor = fff * Fshape * Fcurve * Frotation
        C1M = Q / ((2.0 * np.pi * self.R1 - self.ZI * self.TB) * self.YI1)
        C2M = Q / ((2.0 * np.pi * self.R2 - self.ZI * self.TB) * self.YI2)
        W1 = C1M / np.sin(self.B1)
        W2 = C2M / np.sin(self.B2)
        h_friction = friction_factor * self.LI / self.DH * (W1 ** 2 + W2 ** 2) / (4 * 9.81)
        return h_friction


######################################
# Recirculation loss calculation class
class RecirculationLoss(object):
    def __init__(self, inputs, QBEM):
        self.R1 = inputs['R1']
        self.R2 = inputs['R2']
        self.RD1 = inputs['RD1']
        self.RD2 = inputs['RD2']
        self.TB = inputs['TB']
        self.TV = inputs['TV']
        self.YI1 = inputs['YI1']
        self.YI2 = inputs['YI2']
        self.VOI = inputs['VOI']
        self.VOD = inputs['VOD']
        self.ASF = inputs['ASF']
        self.ASB = inputs['ASB']
        self.AB = inputs['AB']
        self.AV = inputs['AV']
        self.ADF = inputs['ADF']
        self.ADB = inputs['ADB']
        self.LI = inputs['LI']
        self.LD = inputs['LD']
        self.RLK = inputs['RLK']
        self.LG = inputs['LG']
        self.SL = inputs['SL']
        self.EA = inputs['EA']
        self.ZI = inputs['ZI']
        self.ZD = inputs['ZD']
        self.B1 = inputs['B1'] * (pi / 180.0)
        self.B2 = inputs['B2'] * (pi / 180.0)
        self.DENL = inputs['DENL']
        self.VISL = inputs['VISL']
        self.VISW = inputs['VISW']
        self.N = inputs['N']
        self.NS = inputs['NS']
        self.QBEM = QBEM * bbl_to_m3 / 24.0 / 3600.0
        self.OMEGA = 2.0 * pi * self.N / 60.0

    def Thin2008(self, Q):
        if Q > 3500 * 0.15897 / 3600 / 24:
            return 0.
        else:
            return 5e-5 * self.OMEGA ** 3 * (2 * self.R1) ** 2 * (1 - Q / (3500 * 0.15897 / 3600 / 24)) ** 2.5

    def Bing2012(self, Q):
        krec = 8e-5
        b2 = self.YI2
        C1M = Q / ((2.0 * np.pi * self.R1 - self.ZI * self.TB) * self.YI1)
        C2M = Q / ((2.0 * np.pi * self.R2 - self.ZI * self.TB) * self.YI2)
        U2 = self.R2 * self.OMEGA
        W1 = C1M / np.sin(self.B1)
        W2 = C2M / np.sin(self.B2)
        HE = (U2 ** 2 - U2 * C2M / np.tan(self.B2)) / 9.81

        Df = 1 + W2 / W1 * (
            0.75 * 9.81 * HE / (U2 ** 2) / (1 / np.pi * (1 - self.R1 / self.R2 + self.R1 / self.R2)) - 1)
        alpha2 = np.arcsin(Q / (2 * np.pi * self.R2 * b2 * (U2 ** 2 + W2 ** 2 - 2 * U2 * W2 * np.cos(self.B2))))
        h_recirculation = krec * np.sinh(3.5 * alpha2 ** 3) * Df ** 2 * U2 ** 2 / 9.81

        return h_recirculation


###############################
# Shock loss calculation class
class ShockLoss(object):
    def __init__(self, inputs, QBEM):
        self.R1 = inputs['R1']
        self.R2 = inputs['R2']
        self.RD1 = inputs['RD1']
        self.RD2 = inputs['RD2']
        self.TB = inputs['TB']
        self.TV = inputs['TV']
        self.YI1 = inputs['YI1']
        self.YI2 = inputs['YI2']
        self.VOI = inputs['VOI']
        self.VOD = inputs['VOD']
        self.ASF = inputs['ASF']
        self.ASB = inputs['ASB']
        self.AB = inputs['AB']
        self.AV = inputs['AV']
        self.ADF = inputs['ADF']
        self.ADB = inputs['ADB']
        self.LI = inputs['LI']
        self.LD = inputs['LD']
        self.RLK = inputs['RLK']
        self.LG = inputs['LG']
        self.SL = inputs['SL']
        self.EA = inputs['EA']
        self.ZI = inputs['ZI']
        self.ZD = inputs['ZD']
        self.B1 = inputs['B1'] * (pi / 180.0)
        self.B2 = inputs['B2'] * (pi / 180.0)
        self.DENL = inputs['DENL']
        self.VISL = inputs['VISL']
        self.VISW = inputs['VISW']
        self.N = inputs['N']
        self.NS = inputs['NS']
        self.QBEM = QBEM * bbl_to_m3 / 24.0 / 3600.0
        self.OMEGA = 2.0 * pi * self.N / 60.0

    def Amaral2007(self, Q):
        if self.VISL > 270e-3:
            k_shock = 10000
        else:
            k_shock = 20000
        return k_shock * (Q - 2700 * 0.15897 / 24 / 3600) ** 2

    def Bing2012(self, Q):
        if self.VISL > 270e-3:
            k_shock = 1.2
        else:
            k_shock = 0.5

        U1 = self.R1 * self.OMEGA
        return k_shock / (2 * 9.81) * ((Q - 2700 * 0.15897 / 3600 / 24) / (2700 * 0.15897 / 3600 / 24) * U1) ** 2


##################################
# Diffuser loss calculation class
class DiffuserLoss(object):
    def __init__(self, inputs, QBEM):
        self.R1 = inputs['R1']
        self.R2 = inputs['R2']
        self.RD1 = inputs['RD1']
        self.RD2 = inputs['RD2']
        self.TB = inputs['TB']
        self.TV = inputs['TV']
        self.YI1 = inputs['YI1']
        self.YI2 = inputs['YI2']
        self.VOI = inputs['VOI']
        self.VOD = inputs['VOD']
        self.ASF = inputs['ASF']
        self.ASB = inputs['ASB']
        self.AB = inputs['AB']
        self.AV = inputs['AV']
        self.ADF = inputs['ADF']
        self.ADB = inputs['ADB']
        self.LI = inputs['LI']
        self.LD = inputs['LD']
        self.RLK = inputs['RLK']
        self.LG = inputs['LG']
        self.SL = inputs['SL']
        self.EA = inputs['EA']
        self.ZI = inputs['ZI']
        self.ZD = inputs['ZD']
        self.B1 = inputs['B1'] * (pi / 180.0)
        self.B2 = inputs['B2'] * (pi / 180.0)
        self.DENL = inputs['DENL']
        self.VISL = inputs['VISL']
        self.VISW = inputs['VISW']
        self.N = inputs['N']
        self.NS = inputs['NS']
        self.QBEM = QBEM * bbl_to_m3 / 24.0 / 3600.0
        self.OMEGA = 2.0 * pi * self.N / 60.0

        R1 = 0.056
        R2 = 0.017496
        self.R = (R1 + R2) / 2.
        B1 = np.pi / 2.
        B2 = 0
        self.beta = (B1 + B2) / 2.
        a1 = 2 * np.pi * R1 / self.ZD
        b1 = 0.00773
        a2 = 2 * np.pi * R2 / self.ZD
        b2 = 0.010147
        self.ADin = a1 * b1
        self.ADout = a2 * b2
        self.bm = (b1 + b2) / 2.
        self.H = self.bm
        self.DH1 = 2 * a1 * b1 / (a1 + b1)
        self.DH2 = 2 * a2 * b2 / (a2 + b2)
        self.DH = (self.DH1 + self.DH2) / 2.0

    def sun2003_diffuser_loss(self, Q):
        R1 = 0.056
        R2 = 0.017496
        F_garmma = self.sun2003_curvature_effect(Q)
        F_betta = self.sun2003_shape_effect(Q)
        f = self.Churchill_1977(Q)
        h_diffuser = F_garmma * F_betta * f * Q ** 2 / (
            8 * 9.81 * self.DH * np.pi ** 2 * self.bm ** 2 * (np.sin(self.beta)) ** 3) * (R1 - R2) / (R1 * R2)
        return h_diffuser

    def Churchill_1977(self, Q):
        NRe = self.DH * Q * self.DENL / (2 * np.pi * self.R * self.H * self.VISL * np.sin(self.beta))

        if NRe < 2300:
            return 64.0 / NRe
        else:
            return 8 * (2.457 * np.log(1. / ((7. / NRe) ** 0.9 + 0.27 * (self.EA / self.DH)))) ** (-2)

    def sun2003_shape_effect(self, Q):
        # a channel width, b channel height
        a1 = 2 * np.pi * self.R1 / self.ZI
        b1 = self.YI1
        a2 = 2 * np.pi * self.R2 / self.ZI
        b2 = self.YI2
        ll = (np.min([a1, b1]) / np.max([a1, b1]) + np.min([a2, b2]) / np.max([a2, b2])) / 2.0
        NRe = self.DH * Q * self.DENL / (2 * np.pi * self.R * self.H * self.VISL * np.sin(self.beta))

        if NRe < 2300:
            return 1 / (2. / 3. + 11. / 24. * ll * (2 - ll))
        else:
            return 1 / (2. / 3. + 11. / 24. * ll * (2 - ll)) ** 0.25

    def sun2003_curvature_effect(self, Q):
        # radius of curvature of radial type centrifugal pump used in Sun (2003) dissertation
        RC = 1. / np.sin(self.beta) * (1. / np.tan(self.beta) / self.R - (self.B1 - self.B2) / (self.R1 - self.R2)) ** (
            -1)
        rH = self.DH / 2.
        NRe = self.DH * Q * self.DENL / (2 * np.pi * self.R * self.H * self.VISL * np.sin(self.beta))

        if RC / rH >= 860:
            NRec = 2300.
        else:
            NRec = 2e4 * (rH / RC) ** 0.32

        if NRe < NRec:
            if RC / rH >= 860:
                return 1.
            else:
                return 0.266 * NRe ** 0.389 * (rH / RC) ** 0.1945
        else:
            if NRe * (rH / RC) ** 2 >= 300:
                return (NRe * (rH / RC) ** 2) ** 0.05
            elif (NRe * (rH / RC) ** 2 > 0.0034) and (NRe * (rH / RC) ** 2 < 300):
                return 0.092 * (NRe * (rH / RC) ** 2) ** 0.25 + 0.962
            else:
                return 1.

    def Bing2012_diffuser_loss(self, Q):
        F_garmma = self.sun2003_curvature_effect(Q)
        F_betta = self.sun2003_shape_effect(Q)
        f = self.Churchill_1977(Q)
        V2 = Q / self.ADin
        V3 = Q / self.ADout
        h_diffuser = F_garmma * F_betta * f * self.LD / self.DH * (V2 ** 2 + V3 ** 2) / (4 * 9.81)
        return h_diffuser


######################################
# Disk friction loss calculation class
class DiskFrictionLoss(FrictionLoss):
    def __init__(self, inputs):
        super(DiskFrictionLoss, self).__init__(inputs)

    def Kruyt2003(self, Q):
        return 0.5 * self.OMEGA ** 3 * self.R2 ** 5 / (9.81 * Q) / 1000

    def Thin2008_disk_loss(self, Q):
        f_disk = 200
        h_disk = f_disk * self.DENL * self.OMEGA ** 3 * self.R2 ** 5 / (1E9 * Q)
        return h_disk

    def Gulich1999(self, Q):
        h = 0.00049
        f_geo = 0.41
        Re = self.DENL * self.OMEGA * self.R2 ** 2 / self.VISL
        f_RS = 0.63 + 0.6 * 0.056 / self.R2
        k_RR = np.pi * self.R2 / (2 * Re * h) + 0.02 / Re ** 0.2 * (1 + h / self.R2) / (1 + h / (2 * self.R2))
        h_disk = k_RR * self.OMEGA ** 3 * self.R2 ** 5 * f_geo * f_RS / (9.81 * Q)
        return h_disk

    def Ladounai2009(self, Q):
        h = 0.05 * self.R2
        f_therm = np.exp(-2e-5 * (self.VISL / self.DENL / 1e-6)) ** 1.34
        Re = self.DENL * self.OMEGA * self.R2 ** 2 / self.VISL
        k_RR = np.pi * self.R2 / (2 * Re * h) + 0.02 / Re ** 0.2 * (1 + h / self.R2) / (1 + h / (2 * self.R2))
        f_geo = 1 - (self.R1 / self.R2) ** 5
        k_DF = k_RR * f_therm
        h_disk = k_DF * self.OMEGA ** 3 * self.R2 ** 5 * f_geo / (9.81 * Q)
        return h_disk


########################################################################################################################
if __name__ == "__main__":
    TE2700 = {"R1": 0.017496, "R2": 0.056054, "TB": 0.00272, "TV": 0.00448,
              "YI1": 0.012194, "YI2": 0.007835, "VOI": 0.000016119, "VOD": 0.000011153,
              "ASF": 0.00176464, "ASB": 0.00157452, "AB": 0.001319, "AV": 0.001516,
              "ADF": 0.001482, "ADB": 0.000935, "LI": 0.076, "LD": 0.08708,
              "RLK": 0.056209, "LG": 0.00806, "SL": 0.00005, "EA": 0.0000254,
              "ZI": 5, "ZD": 9, "B1": 19.5, "B2": 24.7,
              "DENL": 1000, "DENG": 11.2, "VISL": 0.001, "VISG": 0.000018,
              "VISW": 0.001, "ST": 0.073, "N": 3500, "SGM": 0.3,
              "QL": 2700, "QG": 50}

    print(len(TE2700))

    # sgl = SinglePhaseModel(TE2700, 4800)
    # QL, hpsgl, hesgl, hfsgl, htsgl, hdsgl, hresgl, hlksgl, qlksgl = sgl.performance_curve()

    # # friction loss
    # frictionloss = FrictionLoss(TE2700, 0)
    # hf_Sun2003 = []
    # hf_Thin2008 = []
    # hf_Bing2012 = []
    #
    # # recirculation loss
    # recirculationloss = RecirculationLoss(TE2700, 0)
    # hre_Thin2008 = []
    # hre_Bing2012 = []
    #
    # # shock loss
    # shockloss = ShockLoss(TE2700, 0)
    # ht_Amaral2007 = []
    # ht_Bing2012 = []
    #
    # # Diffuser loss
    # diffuserloss = DiffuserLoss(TE2700, 0)
    # hd_Sun2003 = []
    # hd_Bing2012 = []
    #
    # # Disk loss
    # diskloss = DiskFrictionLoss(TE2700, 0)
    # hdisk_Kruyt2003 = []
    # hdisk_Thin2008 = []
    # hdisk_Gulich1999 = []
    # hdisk_Ladounai2009 = []
    #
    # for q in QL:
    #     q = q * 0.15897 / 3600 / 24
    #
    #     # friction
    #     hf_Sun2003.append(frictionloss.sun2003(q) * 9.81 * 1000 / 6891.2)
    #     hf_Thin2008.append(frictionloss.Thin2008(q) * 9.81 * 1000 / 6891.2)
    #     hf_Bing2012.append(frictionloss.Bing2012(q) * 9.81 * 1000 / 6891.2)
    #
    #     # recirculation
    #     hre_Thin2008.append(recirculationloss.Thin2008(q) * 9.81 * 1000 / 6891.2)
    #     hre_Bing2012.append(recirculationloss.Bing2012(q) * 9.81 * 1000 / 6891.2)
    #
    #     # shock loss
    #     ht_Amaral2007.append(shockloss.Amaral2007(q) * 9.81 * 1000 / 6891.2)
    #     ht_Bing2012.append(shockloss.Bing2012(q) * 9.81 * 1000 / 6891.2)
    #
    #     # diffuser loss
    #     hd_Sun2003.append(diffuserloss.sun2003_diffuser_loss(q) * 9.81 * 1000 / 6891.2)
    #     hd_Bing2012.append(diffuserloss.Bing2012_diffuser_loss(q) * 9.81 * 1000 / 6891.2)
    #
    #     # disk loss
    #     hdisk_Kruyt2003.append(diskloss.Kruyt2003(q) * 9.81 * 1000 / 6891.2)
    #     hdisk_Thin2008.append(diskloss.Thin2008_disk_loss(q) * 9.81 * 1000 / 6891.2)
    #     hdisk_Gulich1999.append(diskloss.Gulich1999(q) * 9.81 * 1000 / 6891.2)
    #     hdisk_Ladounai2009.append(diskloss.Ladounai2009(q) * 9.81 * 1000 / 6891.2)
    #
    # df_friction = pd.DataFrame({'QL': QL, 'friction_mechanistic': hfsgl, 'friction_sun2003': hf_Sun2003,
    #                             'friction_Thin2008': hf_Thin2008, 'friction_Bing2012': hf_Bing2012})
    # # figure1_friction = plt.figure()
    # # ax1 = figure1_friction.add_subplot(111)
    # # ax1.plot(df_friction['QL'], df_friction['friction_mechanistic'], 'bo-', markerfacecolor='None',
    # #         label='Mechanistic model')
    # # ax1.plot(df_friction['QL'][4:], df_friction['friction_sun2003'][4:], 'r^-', markerfacecolor='None',
    # #         label='Sun (2003)')
    # # ax1.plot(df_friction['QL'][4:], df_friction['friction_Thin2008'][4:], 'ks-', markerfacecolor='None',
    # #         label='Thin et al. (2008)')
    # # ax1.plot(df_friction['QL'][4:], df_friction['friction_Bing2012'][4:], 'gd-', markerfacecolor='None',
    # #         label='Bing et al. (2012)')
    # #
    # # ax1.set_xlabel(r'$Q_{L}$ (bpd)')
    # # ax1.set_ylabel(r'P (psi)')
    # # ax1.tick_params(direction='in')
    # # ax1.legend(frameon=False)
    # # figure1_friction.savefig('Friction_model.jpg')
    #
    # df_recirculation = pd.DataFrame(
    #     {'QL': QL, 'recirculation_mechanistic': hresgl, 'recirculation_Thin2008': hre_Thin2008,
    #      'recirculation_Bing2012': hre_Bing2012})
    # # figure2_recirculation = plt.figure()
    # # ax2 = figure2_recirculation.add_subplot(111)
    # # ax2.plot(df_recirculation['QL'], df_recirculation['recirculation_mechanistic'], 'bo-', markerfacecolor='None',
    # #                  label='Mechanistic model')
    # # ax2.plot(df_recirculation['QL'], df_recirculation['recirculation_Thin2008'], 'r^-', markerfacecolor='None',
    # #          label='Thin et al. (2008)')
    # # ax2.plot(df_recirculation['QL'], df_recirculation['recirculation_Bing2012'], 'ks-', markerfacecolor='None',
    # #          label='Bing et al. (2012)')
    # #
    # # ax2.set_xlabel(r'$Q_{L}$ (bpd)')
    # # ax2.set_ylabel(r'P (psi)')
    # # ax2.tick_params(direction='in')
    # # ax2.legend(frameon=False)
    # # figure2_recirculation.savefig('Recirculation_model.jpg')
    #
    # df_shock = pd.DataFrame({'QL': QL, 'shock_mechanistic': htsgl, 'shock_Amaral2007': ht_Amaral2007,
    #                          'shock_Bing2012': ht_Bing2012})
    # # figure3_shock = plt.figure()
    # # ax3 = figure3_shock.add_subplot(111)
    # # ax3.plot(df_shock['QL'], df_shock['shock_mechanistic'], 'bo-', markerfacecolor='None',
    # #          label='Mechanistic model')
    # # ax3.plot(df_shock['QL'], df_shock['shock_Amaral2007'], 'r^-', markerfacecolor='None',
    # #          label='Amaral et al. (2007)')
    # # ax3.plot(df_shock['QL'], df_shock['shock_Bing2012'], 'ks-', markerfacecolor='None',
    # #          label='Bing et al. (2012)')
    # #
    # # ax3.set_xlabel(r'$Q_{L}$ (bpd)')
    # # ax3.set_ylabel(r'P (psi)')
    # # ax3.tick_params(direction='in')
    # # ax3.legend(frameon=False)
    # # figure3_shock.savefig('Shock_model.jpg')
    #
    # df_diffuser = pd.DataFrame({'QL': QL, 'diffuser_mechanistic': hdsgl, 'diffuser_Sun2003': hd_Sun2003,
    #                             'diffuser_Bing2012': hd_Bing2012})
    # # figure4_diffuser = plt.figure()
    # # ax4 = figure4_diffuser.add_subplot(111)
    # # ax4.plot(df_diffuser['QL'], df_diffuser['diffuser_mechanistic'], 'bo-', markerfacecolor='None',
    # #          label='Mechanistic model')
    # # ax4.plot(df_diffuser['QL'], df_diffuser['diffuser_Sun2003'], 'r^-', markerfacecolor='None',
    # #          label='Sun (2003)')
    # # ax4.plot(df_diffuser['QL'], df_diffuser['diffuser_Bing2012'], 'ks-', markerfacecolor='None',
    # #          label='Bing et al. (2012)')
    # #
    # # ax4.set_xlabel(r'$Q_{L}$ (bpd)')
    # # ax4.set_ylabel(r'P (psi)')
    # # ax4.tick_params(direction='in')
    # # ax4.legend(frameon=False)
    # # figure4_diffuser.savefig('Diffuser_model.jpg')
    #
    # df_disk = pd.DataFrame({'QL': QL, 'disk_Kruyt2003': hdisk_Kruyt2003, 'disk_Thin2008': hdisk_Thin2008,
    #                         'disk_Gulich1999': hdisk_Gulich1999, 'disk_Ladounani2009': hdisk_Ladounai2009})
    # # figure5_disk = plt.figure()
    # # ax5 = figure5_disk.add_subplot(111)
    # # ax5.plot(df_disk['QL'], df_disk['disk_Gulich1999'], 'ks-', markerfacecolor='None',
    # #          label='Gulich (1999)')
    # # ax5.plot(df_disk['QL'], df_disk['disk_Kruyt2003'], 'bo-', markerfacecolor='None',
    # #          label='Kruyt (2003)')
    # # ax5.plot(df_disk['QL'], df_disk['disk_Thin2008'], 'r^-', markerfacecolor='None',
    # #          label='Thin et al. (2008)')
    # # ax5.plot(df_disk['QL'], df_disk['disk_Ladounani2009'], 'gd-', markerfacecolor='None',
    # #          label='Ladounani(2009)')
    # #
    # # ax5.set_xlabel(r'$Q_{L}$ (bpd)')
    # # ax5.set_ylabel(r'P (psi)')
    # # ax5.tick_params(direction='in')
    # # ax5.legend(frameon=False)
    # # figure5_disk.savefig('Disk_model.jpg')
    #
    # figure6_mechanistic_model = plt.figure()
    # ax6 = figure6_mechanistic_model.add_subplot(111)
    # ax6.plot(QL, hesgl, 'bs-', label='Euler head', markerfacecolor='None')
    # ax6.plot(QL, hfsgl, 'rd-', label='friction loss', markerfacecolor='None')
    # ax6.plot(QL, hresgl, 'gv-', label='recirculation loss', markerfacecolor='None')
    # ax6.plot(QL, htsgl, 'c^-', label='turn loss', markerfacecolor='None')
    # ax6.plot(QL, hlksgl, 'm*-', label='leakage loss', markerfacecolor='None')
    # ax6.plot(QL, hpsgl, 'ko-', label='predicted head', markerfacecolor='None')
    #
    # ax6.set_xlabel(r'$Q_{L}$ (bpd)')
    # ax6.set_ylabel(r'P (psi)')
    # ax6.set_ylim([0, 80])
    # ax6.tick_params(direction='in')
    # ax6.legend(frameon=False)
    #
    # # figure6_mechanistic_model.savefig('Losses_in_mechanistic_model_new.jpg')
    #
    # # plt.show()
