# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from losses import FrictionLoss

"""
The ESP calculation module contains two classes for computing inputs performance under single-phase water/viscous
fluid flow or gas-liquid two-phase flow conditions.

The mechanistic model was originally proposed by Dr. Zhang, TUALP ABM (2013) and later revised by Zhu (2017)

Version:    1st, Aug, 2017
Developer:  Jianjun Zhu
"""

# global variables
G = 9.81
pi = np.pi
E1 = 1e-5
DENW = 997.
psi_to_pa = 1.013e5 / 14.7
bbl_to_m3 = 0.15897


########################################################################################################################
# single-phase mechanistic model
class SinglePhaseModel(object):
    def __init__(self, inputs, QBEM):
        """
        Typical inputs: Q in bpd, N in rpm, angle in degree, all others in SI units
        inputs = {  "R1": 0.017496, "R2": 0.056054, "TB": 0.00272, "TV": 0.00448, "RD1": 0.0219507825, "RD2": 0.054735,
                    "YI1": 0.012194, "YI2": 0.007835, "VOI": 0.000016119, "VOD": 0.000011153, "ASF": 0.00176464,
                    "ASB": 0.00157452, "AB": 0.001319, "AV": 0.001516, "ADF": 0.001482, "ADB": 0.000935, "LI": 0.076,
                    "LD": 0.08708, "RLK": 0.056209, "LG": 0.00806, "SL": 0.00005, "EA": 0.0000254, "ZI": 5, "ZD": 9,
                    "B1": 19.5, "B2": 24.7, "NS": 1600, "DENL": 1000, "DENG": 11.2, "VISL": 0.001, "VISG": 0.000018,
                    "VISW": 0.001, "ST": 0.073, "N": 3500, "SGM": 0.3, "QL": 2700, "QG": 50
                }
        """
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

        # use literature loss models
        self.friction_loss = FrictionLoss(inputs)

    @staticmethod
    def get_fff(re, ed):
        # friction factor used in unified model
        lo = 1000.
        hi = 3000.
        REM = re
        ED = ed
        Fl = 16.0 / REM
        Fh = 0.07716 / (np.log(6.9 / REM + (ED / 3.7)**1.11))**2.0

        if REM < lo:
            return Fl
        elif REM > hi:
            return Fh
        else:
            return (Fh * (REM - lo) + Fl * (hi - REM)) / (hi - lo)

    @staticmethod
    def get_fff_jc(re, ed):
        # friction factor based on Churchill correlation
        REM = re
        AA = (2.457 * np.log(1.0 / ((7.0 / REM) ** 0.9 + 0.27 * ed))) ** 16.0
        BB = (37530.0 / REM) ** 16.0
        fff = 2.0 * ((8.0 / REM) ** 12.0 + 1.0 / (AA + BB) ** 1.5) ** (1.0 / 12.0) / 4
        return fff

    @staticmethod
    def emulsion_zhang():
        pass

    @staticmethod
    def emulsion_banjar():
        pass

    def sgl_calculate_old(self, Q):
        """
        Original development on single-phase model by Zhu and Zhang
        :param Q: flow rate in bpd
        :return: HP, HE, HF, HT, HD, HRE, HLK, QLK
        """
        # Q in bpd
        Q = Q * bbl_to_m3 / 24.0 / 3600.0
        QLK = 0.02 * Q
        HP = 10.
        HE = 0.
        HEE = 0.
        HFI = 0.
        HFD = 0.
        HTI = 0.
        HTD = 0.
        HLK = 0.
        icon = 0
        QLK_new = 1.
        HP_new = 1.
        ABH = np.abs((HP - HP_new) / HP_new)
        ABQ = np.abs((QLK - QLK_new) / HP_new)
        AIW = self.AB + self.ASF + self.ASB
        ADW = self.AV + self.ADF + self.ADB
        SGM = 0.

        while (ABH > E1) or (ABQ > E1):
            C1M = (Q + QLK) / ((2.0 * pi * self.R1 - self.ZI * self.TB) * self.YI1)
            C2M = (Q + QLK) / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            U1 = self.R1 * self.OMEGA
            U2 = self.R2 * self.OMEGA
            W1 = C1M / np.sin(self.B1)
            W2 = C2M / np.sin(self.B2)
            C1 = np.sqrt(C1M**2 + (U1 - C1M / np.tan(self.B1))**2)
            C2 = np.sqrt(C2M**2 + (U2 - C2M / np.tan(self.B2))**2)
            CMB = self.QBEM / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            C2B = np.sqrt(CMB**2 + (U2 - CMB / np.tan(self.B2))**2)

            # Euler head
            # HE=(U2**2-U1**2+W1**2-W2**2+C2**2-C1**2)/(2.0*G)					# with pre-rotation
            HE = (U2**2 - U2 * C2M / np.tan(self.B2)) / G 						# without pre-rotation

            # head loss due to recirculation
            if (Q + QLK) < self.QBEM:
                VSH = U2 * (self.QBEM - (Q + QLK)) / self.QBEM
                C2F = C2B * (Q + QLK) / self.QBEM
                DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
                REC = self.DENL * VSH * DC / self.VISL

                if np.abs(self.VISL - self.VISW) > 0.001:
                    SGM = (self.VISW / self.VISL)**0.5 / (1.0 + 0.2 * REC**0.2)  # consider viscosity

                C2P = (C2**2 + C2F**2 - VSH**2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F)
                HEE = HE + (C2E**2 - C2**2) / (2.0 * G)
            else:
                VSH = U2 * (Q + QLK - self.QBEM) / self.QBEM
                C2F = C2B * (Q + QLK) / self.QBEM
                C2E = (C2**2 + C2F**2 - VSH**2) / (2.0 * C2F)
                # C2E = C2F + SGM * (C2P - C2F) * (Q + QLK - self.QBEM) / self.QBEM     # new development by Dr Zhang
                HEE = HE + (C2E**2 - C2**2) / (2.0 * G)

            # friction loss
            AI = self.VOI / self.LI
            AD = self.VOD / self.LD
            VI = (Q + QLK) / self.ZI / AI
            VD = Q / self.ZD / AD
            DI = 4.0 * self.VOI / AIW
            DD = 4.0 * self.VOD / ADW
            REI = self.DENL * (W1 + W2) * DI / self.VISL / 2.0
            RED = self.DENL * VD * DD / self.VISL
            EDI = self.EA / DI
            EDD = self.EA / DD
            FFI = self.get_fff(REI, EDI)
            FFD = self.get_fff(RED, EDD)
            HFI = 5 * 4.0 * FFI * (W1 + W2)**2 * self.LI / (8.0 * G * DI)
            HFD = 5 * 4.0 * FFD * VD**2 * self.LD / (2.0 * G * DD)

            # turn loss
            FTI = 2.5
            FTD = 2.5
            HTI = FTI * VI**2 / (2.0 * G)
            HTD = FTD * VD**2 / (2.0 * G)

            # new pump head
            HP_new = HEE - HFI - HFD - HTI - HTD

            # calculate leakage
            UL = self.RLK * self.OMEGA
            HIO = HEE - HFI - HTI
            HLK = HIO - (U2**2 - UL**2) / (8.0 * G)
            if HLK >= 0:
                VL = np.abs(QLK) / (2.0 * pi * self.RLK * self.SL)
                REL = self.DENL * VL * self.SL / self.VISL
                EDL = 0.0
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * G * HLK / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLK_new = 2.0 * pi * self.RLK * self.SL * VL
            else:
                VL = np.abs(QLK / (2.0 * pi * self.RLK * self.SL))
                REL = self.DENL * VL * self.SL / self.VISL
                EDL = 0.
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * G * np.abs(HLK) / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLK_new = -2.0 * pi * self.RLK * self.SL * VL

            ABQ = np.abs((QLK_new - QLK) / QLK_new)
            QLK = QLK_new
            ABH = np.abs((HP_new - HP) / HP_new)
            HP = HP_new

            if icon > 500:
                break
            else:
                icon += 1

        # return pressure in psi, flow rate in bpd
        HP = HP * G * DENW / psi_to_pa
        HE = HE * G * DENW / psi_to_pa
        HEE = HEE * G * DENW / psi_to_pa
        HF = (HFI + HFD) * G * DENW / psi_to_pa
        HT = (HTI + HTD) * G * DENW / psi_to_pa
        HD = (HFD + HTD) * G * DENW / psi_to_pa
        HRE = np.abs(HE - HEE)
        HLK = HLK * G * DENW / psi_to_pa
        QLK = QLK * 24.0 * 3600.0 / bbl_to_m3
        return HP, HE, HF, HT, HD, HRE, HLK, QLK

    def sgl_calculate_new(self, Q):
        """
        Dr Zhang new update on previous single-phase model with new consideration on recirculation flow loss
        :param Q: flow rate in bpd
        :return: HP, HE, HF, HT, HD, HRE, HLK, QLK
        """
        # Q in bpd
        Q = Q * bbl_to_m3 / 24.0 / 3600.0
        QLK = 0.02 * Q
        HP = 10.
        HE = 0.
        HEE = 0.
        HFI = 0.
        HFD = 0.
        HTI = 0.
        HTD = 0.
        HLK = 0.
        icon = 0
        QLK_new = 1.
        HP_new = 1.
        ABH = np.abs((HP - HP_new) / HP_new)
        ABQ = np.abs((QLK - QLK_new) / HP_new)
        AIW = self.AB + self.ASF + self.ASB
        ADW = self.AV + self.ADF + self.ADB
        SGM = 0.

        while (ABH > E1) or (ABQ > E1):
            C1M = (Q + QLK) / ((2.0 * pi * self.R1 - self.ZI * self.TB) * self.YI1)
            C2M = (Q + QLK) / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            U1 = self.R1 * self.OMEGA
            U2 = self.R2 * self.OMEGA
            W1 = C1M / np.sin(self.B1)
            W2 = C2M / np.sin(self.B2)
            C1 = np.sqrt(C1M ** 2 + (U1 - C1M / np.tan(self.B1)) ** 2)
            C2 = np.sqrt(C2M ** 2 + (U2 - C2M / np.tan(self.B2)) ** 2)
            CMB = self.QBEM / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            C2B = np.sqrt(CMB ** 2 + (U2 - CMB / np.tan(self.B2)) ** 2)

            # Euler head
            # HE=(U2**2-U1**2+W1**2-W2**2+C2**2-C1**2)/(2.0*G)					        # with pre-rotation
            HE = (U2 ** 2 - U2 * C2M / np.tan(self.B2)) / G                             # without pre-rotation

            # head loss due to recirculation
            if (Q + QLK) < self.QBEM:
                VSH = U2 * (self.QBEM - (Q + QLK)) / self.QBEM
                C2F = C2B * (Q + QLK) / self.QBEM
                DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
                REC = self.DENL * VSH * DC / self.VISL
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F
                HEE = HE + (C2E ** 2 - C2 ** 2) / (2.0 * G)
            else:
                VSH = U2 * (Q + QLK - self.QBEM) / self.QBEM
                C2F = C2B * (Q + QLK) / self.QBEM
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F) * (Q + QLK - self.QBEM) / self.QBEM     # new development by Dr Zhang
                HEE = HE + (C2E ** 2 - C2 ** 2) / (2.0 * G)

            # friction loss
            AI = self.VOI / self.LI
            AD = self.VOD / self.LD
            VI = (Q + QLK) / self.ZI / AI
            VD = Q / self.ZD / AD
            DI = 4.0 * self.VOI / AIW
            DD = 4.0 * self.VOD / ADW
            REI = self.DENL * (W1 + W2) * DI / self.VISL / 2.0
            RED = self.DENL * VD * DD / self.VISL
            EDI = self.EA / DI
            EDD = self.EA / DD
            FFI = self.get_fff(REI, EDI)
            FFD = self.get_fff(RED, EDD)
            HFI = 2.5 * 4.0 * FFI * (W1 + W2) ** 2 * self.LI / (8.0 * G * DI)
            HFD = 2.5 * 4.0 * FFD * VD ** 2 * self.LD / (2.0 * G * DD)

            # turn loss
            FTI = 3.0
            FTD = 3.0
            HTI = FTI * VI ** 2 / (2.0 * G)
            HTD = FTD * VD ** 2 / (2.0 * G)

            # new pump head
            HP_new = HEE - HFI - HFD - HTI - HTD

            # calculate leakage
            UL = self.RLK * self.OMEGA
            HIO = HEE - HFI - HTI
            HLK = HIO - (U2 ** 2 - UL ** 2) / (8.0 * G)
            if HLK >= 0:
                VL = np.abs(QLK) / (2.0 * pi * self.RLK * self.SL)
                REL = self.DENL * VL * self.SL / self.VISL
                EDL = 0.0
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * G * HLK / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLK_new = 2.0 * pi * self.RLK * self.SL * VL
            else:
                VL = np.abs(QLK / (2.0 * pi * self.RLK * self.SL))
                REL = self.DENL * VL * self.SL / self.VISL
                EDL = 0.
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * G * np.abs(HLK) / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLK_new = -2.0 * pi * self.RLK * self.SL * VL

            ABQ = np.abs((QLK_new - QLK) / QLK_new)
            QLK = QLK_new
            ABH = np.abs((HP_new - HP) / HP_new)
            HP = HP_new

            if icon > 500:
                break
            else:
                icon += 1

        # return pressure in psi, flow rate in bpd
        HP = HP * G * DENW / psi_to_pa
        HE = HE * G * DENW / psi_to_pa
        HEE = HEE * G * DENW / psi_to_pa
        HF = (HFI + HFD) * G * DENW / psi_to_pa
        HT = (HTI + HTD) * G * DENW / psi_to_pa
        HD = (HFD + HTD) * G * DENW / psi_to_pa
        HRE = np.abs(HE - HEE)
        HLK = HLK * G * DENW / psi_to_pa
        QLK = QLK * 24.0 * 3600.0 / bbl_to_m3
        return HP, HE, HF, HT, HD, HRE, HLK, QLK

    def sgl_calculate_jc(self, Q):
        """
        Jiecheng Zhang (2017) thesis development
        :param Q: flow rate in bpd
        :return: HP, HE, HF, HT, HD, HRE, HLK, QLK
        """
        # Q in bpd
        Q = Q * bbl_to_m3 / 24.0 / 3600.0
        QLK = 0.02 * Q
        HP = 10.
        HE = 0.
        HEE = 0.
        HFI = 0.
        HFD = 0.
        HTI = 0.
        HTD = 0.
        HLK = 0.
        icon = 0
        QLK_new = 1.
        HP_new = 1.
        ABH = np.abs((HP - HP_new) / HP_new)
        ABQ = np.abs((QLK - QLK_new) / HP_new)
        AIW = self.AB + self.ASF + self.ASB
        ADW = self.AV + self.ADF + self.ADB

        # slip factor compared to Wiesner (1961)
        SGMU = 1 - np.sqrt(np.sin(self.B2)) / self.ZI ** (1.5 * (3448 / self.NS) ** 0.4)

        # new QBEM due to liquid viscosity
        QBEM = self.QBEM * (self.N / 3600) * (self.VISL / self.VISW) ** (0.01 * (3448 / self.NS))

        for icon in range(1, 100):
            C1M = (Q + QLK) / ((2.0 * pi * self.R1 - self.ZI * self.TB) * self.YI1)
            C2M = (Q + QLK) / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            U1 = self.R1 * self.OMEGA
            U2 = self.R2 * self.OMEGA
            W1 = C1M / np.sin(self.B1)
            W2 = C2M / np.sin(self.B2)
            C1 = np.sqrt(C1M ** 2 + (U1 - C1M / np.tan(self.B1)) ** 2)
            C2 = np.sqrt(C2M ** 2 + (U2 - C2M / np.tan(self.B2)) ** 2)
            CMB = QBEM / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            C2B = np.sqrt(CMB ** 2 + (U2 - CMB / np.tan(self.B2)) ** 2)
            HE = (U2 ** 2 * SGMU - U2 * C2M / np.tan(self.B2)) / G

            if (Q+QLK) <= QBEM:
                VSH = U2 * (QBEM - (Q + QLK)) / QBEM
                C2F = C2B * (Q + QLK) / QBEM
                DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
                REC = self.DENL * VSH * DC / self.VISL
                SGM = (self.VISW / self.VISL) ** 0.1 / (10.0 + 0.02 * REC ** 0.2)   # SGM: shear factor due to viscosity
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F)
                HEE = HE + (C2E ** 2 - C2 ** 2) / (2.0 * G)
            else:
                VSH = U2 * (Q + QLK - QBEM) / QBEM
                C2F = C2B * (Q + QLK) / QBEM
                DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
                REC = self.DENL * VSH * DC / self.VISL
                SGM = (self.VISW / self.VISL) ** 0.1 / (10.0 + 0.02 * REC ** 0.2)
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F) * (Q + QLK - QBEM) / QBEM
                HEE = HE + (C2E ** 2 - C2 ** 2) / (2.0 * G)

            AI = self.VOI / self.LI
            AD = self.VOD / self.LD
            VI = (Q + QLK) / self.ZI / AI
            VD = Q / self.ZD / AD
            DI = 4.0 * self.VOI / AIW
            DD = 4.0 * self.VOD / ADW
            REI = self.DENL * VI * DI / self.VISL
            RED = self.DENL * VD * DD / self.VISL / 1.0
            EDI = self.EA / DI
            EDD = self.EA / DD
            FFI = self.get_fff_jc(REI, EDI)
            FFD = self.get_fff_jc(RED, EDD)

            FFIC = 10
            FFDC = FFIC / ((self.R2 - self.R1) / self.LI) * (abs(self.RD2 - self.RD1) / self.LD)
            HFI = 4.0 * FFI * (W1 + W2) ** 2 * self.LI / (8.0 * G * DI) * FFIC
            HFD = 1 * 4.0 * FFD * VD ** 2 * self.LD / (2.0 * G * DD) * FFDC

            FTI = 6
            FTD = FTI / ((self.R2 - self.R1) / self.LI) * (abs(self.RD2 - self.RD1) / self.LD)
            HTI = FTI * VI ** 2 / (2.0 * G)
            HTD = FTD * VD ** 2 / (2.0 * G)
            HPN = HEE - HFI - HFD - HTI - HTD

            UL = self.RLK * self.OMEGA
            HIO = HEE - HFI - HTI
            HLK = HIO - (U2 ** 2 - UL ** 2) / (8.0 * G)

            if HLK >= 0.:
                VL = abs(QLK) / (2.0 * pi * self.RLK * self.SL)
                REL = self.DENL * VL * self.SL / self.VISL
                EDL = 0.0
                FFL = self.get_fff_jc(REL, EDL)
                VL = np.sqrt(2.0 * G * HLK / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLKN = 2.0 * pi * self.RLK * self.SL * VL
            else:
                VL = abs(QLK / (2.0 * pi * self.RLK * self.SL))
                REL = self.DENL * VL * self.SL / self.VISL
                EDL = 0.0
                FFL = self.get_fff_jc(REL, EDL)
                VL = np.sqrt(2.0 * G * abs(HLK) / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLKN = -2.0 * pi * self.RLK * self.SL * VL

            ABL = abs((QLKN - QLK) / QLKN)
            QLK = QLKN
            ABC = abs((HPN - HP) / HPN)
            HP = HPN

            if (ABC < E1) and (ABL < E1):
                break

        # return pressure in psi, flow rate in bpd
        HP = HP * G * DENW / psi_to_pa
        HE = HE * G * DENW / psi_to_pa
        HEE = HEE * G * DENW / psi_to_pa
        HF = (HFI + HFD) * G * DENW / psi_to_pa
        HT = (HTI + HTD) * G * DENW / psi_to_pa
        HD = (HFD + HTD) * G * DENW / psi_to_pa
        HRE = np.abs(HE - HEE)
        HLK = HLK * G * DENW / psi_to_pa
        QLK = QLK * 24.0 * 3600.0 / bbl_to_m3
        return HP, HE, HF, HT, HD, HRE, HLK, QLK

    def sgl_calculate_2018(self, Q):
        """
        Based on previous single-phase model (Jiecheng Zhang) with new development on pressure losses
        :param Q: flow rate in bpd
        :return: HP, HE, HF, HT, HD, HRE, HLK, QLK
        """
        # Q in bpd
        Q = Q * bbl_to_m3 / 24.0 / 3600.0
        QLK = 0.02 * Q
        HP = 10.
        HE = 0.
        HEE = 0.
        HFI = 0.
        HFD = 0.
        HTI = 0.
        HTD = 0.
        HLK = 0.
        icon = 0
        QLK_new = 1.
        HP_new = 1.
        ABH = np.abs((HP - HP_new) / HP_new)
        ABQ = np.abs((QLK - QLK_new) / HP_new)
        AIW = self.AB + self.ASF + self.ASB
        ADW = self.AV + self.ADF + self.ADB

        # slip factor compared to Wiesner (1961)
        SGMU = 1 - np.sqrt(np.sin(self.B2)) / self.ZI ** (1.5 * (3448 / self.NS) ** 0.4)
        SGM = 0.

        # new QBEM due to liquid viscosity
        QBEM = self.QBEM * (self.N / 3600) * (self.VISL / self.VISW) ** (0.01 * (3448 / self.NS))

        while (ABH > E1) or (ABQ > E1):
            C1M = (Q + QLK) / ((2.0 * pi * self.R1 - self.ZI * self.TB) * self.YI1)
            C2M = (Q + QLK) / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            U1 = self.R1 * self.OMEGA
            U2 = self.R2 * self.OMEGA
            W1 = C1M / np.sin(self.B1)
            W2 = C2M / np.sin(self.B2)
            C1 = np.sqrt(C1M ** 2 + (U1 - C1M / np.tan(self.B1)) ** 2)
            C2 = np.sqrt(C2M ** 2 + (U2 - C2M / np.tan(self.B2)) ** 2)
            CMB = QBEM / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            C2B = np.sqrt(CMB ** 2 + (U2 - CMB / np.tan(self.B2)) ** 2)

            # Euler head
            # HE=(U2**2-U1**2+W1**2-W2**2+C2**2-C1**2)/(2.0*G)
            # HE = (U2 ** 2 - U2 * C2M / np.tan(self.B2)) / G
            HE = (U2 ** 2 * SGMU - U2 * C2M / np.tan(self.B2)) / G

            # head loss due to recirculation
            if (Q+QLK) <= QBEM:
                VSH = U2 * (QBEM - (Q + QLK)) / QBEM
                C2F = C2B * (Q + QLK) / QBEM
                DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
                REC = self.DENL * VSH * DC / self.VISL
                
                # SGM: shear factor due to viscosity
                SGM = (self.VISW / self.VISL) ** 0.1 / (10.0 + 0.02 * REC ** 0.2)
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F)
                HEE = HE + (C2E ** 2 - C2 ** 2) / (2.0 * G)
            else:
                VSH = U2 * (Q + QLK - QBEM) / QBEM
                C2F = C2B * (Q + QLK) / QBEM
                DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
                REC = self.DENL * VSH * DC / self.VISL
                SGM = (self.VISW / self.VISL) ** 0.1 / (10.0 + 0.02 * REC ** 0.2)
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F) * (Q + QLK - QBEM) / QBEM
                HEE = HE + (C2E ** 2 - C2 ** 2) / (2.0 * G)

            # friction loss
            AI = self.VOI / self.LI
            AD = self.VOD / self.LD
            VI = (Q + QLK) / self.ZI / AI
            VD = Q / self.ZD / AD
            DI = 4.0 * self.VOI / AIW
            DD = 4.0 * self.VOD / ADW
            REI = self.DENL * (W1 + W2) * DI / self.VISL / 2.0
            RED = self.DENL * VD * DD / self.VISL
            EDI = self.EA / DI
            EDD = self.EA / DD
            FFI = self.get_fff(REI, EDI)
            FFD = self.get_fff(RED, EDD)
            # HFI = 2.5 * 4.0 * FFI * (W1 + W2) ** 2 * self.LI / (8.0 * G * DI)
            HFI = self.friction_loss.sun2003(Q)                         # use Sun(2003) friction model
            HFD = 2.5 * 4.0 * FFD * VD ** 2 * self.LD / (2.0 * G * DD)

            # turn loss
            FTI = 3.0
            FTD = 3.0
            HTI = FTI * VI ** 2 / (2.0 * G)
            HTD = FTD * VD ** 2 / (2.0 * G)

            # new pump head
            HP_new = HEE - HFI - HFD - HTI - HTD

            # calculate leakage
            UL = self.RLK * self.OMEGA
            HIO = HEE - HFI - HTI
            HLK = HIO - (U2 ** 2 - UL ** 2) / (8.0 * G)
            if HLK >= 0:
                VL = np.abs(QLK) / (2.0 * pi * self.RLK * self.SL)
                REL = self.DENL * VL * self.SL / self.VISL
                EDL = 0.0
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * G * HLK / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLK_new = 2.0 * pi * self.RLK * self.SL * VL
            else:
                VL = np.abs(QLK / (2.0 * pi * self.RLK * self.SL))
                REL = self.DENL * VL * self.SL / self.VISL
                EDL = 0.
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * G * np.abs(HLK) / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLK_new = -2.0 * pi * self.RLK * self.SL * VL

            ABQ = np.abs((QLK_new - QLK) / QLK_new)
            QLK = QLK_new
            ABH = np.abs((HP_new - HP) / HP_new)
            HP = HP_new

            if icon > 500:
                break
            else:
                icon += 1

        # return pressure in psi, flow rate in bpd
        HP = HP * G * DENW / psi_to_pa
        HE = HE * G * DENW / psi_to_pa
        HEE = HEE * G * DENW / psi_to_pa
        HF = (HFI + HFD) * G * DENW / psi_to_pa
        HT = (HTI + HTD) * G * DENW / psi_to_pa
        HD = (HFD + HTD) * G * DENW / psi_to_pa
        HRE = np.abs(HE - HEE)
        HLK = HLK * G * DENW / psi_to_pa
        QLK = QLK * 24.0 * 3600.0 / bbl_to_m3
        return HP, HE, HF, HT, HD, HRE, HLK, QLK

    def performance_curve(self):
        QL = np.arange(500)[1:] * 150.0
        hpsgl = []
        hesgl = []
        hfsgl = []
        htsgl = []
        hdsgl = []
        hresgl = []
        hlksgl = []
        qlksgl = []

        for ql in QL:
            # q in bpd
            # HP, HE, HF, HT, HD, HRE, HLK, QLK = self.SglCalculate_old(ql)
            HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_2018(ql)

            if HP > 0:
                hpsgl.append(HP)
                hesgl.append(HE)
                hfsgl.append(HF)
                htsgl.append(HT)
                hdsgl.append(HD)
                hresgl.append(HRE)
                hlksgl.append(HLK)
                qlksgl.append(QLK)
            else:
                break
        return QL[:len(hpsgl)], hpsgl, hesgl, hfsgl, htsgl, hdsgl, hresgl, hlksgl, qlksgl


########################################################################################################################
# two-phase mechanistic model
class GasLiquidModel(SinglePhaseModel):
    def __init__(self, inputs, QBEM):
        super(GasLiquidModel, self).__init__(inputs, QBEM)
        self.DENG = inputs['DENG']
        self.VISG = inputs['VISG']
        self.ST = inputs['ST']
        self.QL = inputs['QL'] * bbl_to_m3 / 24.0 / 3600.0
        self.QG = inputs['QG'] * bbl_to_m3 / 24.0 / 3600.0
        self.RI = (self.R1 + self.R2) / 2.0
        self.YI = (self.YI1 + self.YI2) / 2.0

    def get_lambda_c1(self, Q, HP):
        # Q in bpd, HP in psi
        Q = Q * bbl_to_m3 / 24.0 / 3600.0
        HP = HP * psi_to_pa

        # critical bubble size in turbulent flow (Barnea 1982)
        DBCrit = 2.0 * (0.4 * 0.073 / (self.DENL - self.DENG) / (self.OMEGA**2 * self.RI))**0.5
        LambdaC1 = DBCrit / (6.034 / 0.6 * (self.ST / self.DENL)**0.6 * (HP * Q / (self.DENL * self.ZI * self.VOI))**(-0.4) * (self.DENL / self.DENG)**0.2)
        return LambdaC1

    def get_lambda_c2(self, Q, HP):
        # Q in bpd, HP in psi
        Q = Q * bbl_to_m3 / 24.0 / 3600.0
        HP = HP * psi_to_pa
        alphaG_crit = pi / 6.0 - (pi / 6.0 - 0.25) * np.exp(-self.N / 3600.0)
        VSR = 0.1
        LambdaC2 = 1.0
        ABV = 1.0
        icon = 0

        while ABV > 0.001:
            LambdaC2 -= 1e-5
            DB = 6.034 * LambdaC2 * (self.ST / self.DENL)**0.6 * (HP * Q / (self.DENL * self.ZI * self.VOI))**(-0.4) * (self.DENL / self.DENG)**0.2
            DBMAX = DB / 0.6
            REB = self.DENL * VSR * DBMAX / self.VISL
            SR = DBMAX * self.OMEGA / VSR

            if REB <= 50:
                CD = 24.0 / REB * (1.0 + 0.15 * REB**0.687) * (1.0 + 0.3 * SR**2.5)
            else:
                CD = 24.0 / REB * (1.0 + 0.15 * REB**0.687) * (1.0 + 0.55 * SR**2.0)

            VSR = np.sqrt(4.0 * DBMAX * (self.DENL - self.DENG) * self.RI / (3.0 * CD * self.DENL)) * self.OMEGA
            RS = VSR * (2.0 * pi * self.RI - self.ZI * self.TB) * self.YI / Q
            alphaG = (RS - 1.0 + np.sqrt((1.0 - RS)**2 + 4.0 * RS * LambdaC2)) / (2.0 * RS)
            ABV = np.abs((alphaG - alphaG_crit) / alphaG_crit)

            if LambdaC2 < 0:
                LambdaC2 = -LambdaC2
                return LambdaC2

            if icon > 1e5:
                return LambdaC2
            else:
                icon += 1
        return LambdaC2

    def get_lambda_c3(self, Q):
        # Q in bpd
        Q = Q * bbl_to_m3 / 24.0 / 3600.0
        LambdaC3 = 0.
        ANG = -90.0 * pi / 180.0
        FIC = 0.0142
        BI = (self.B1 + self.B2) / 2.0
        AI = self.VOI / self.LI
        AIW = self.AB + self.ASF + self.ASB
        DI = 4.0 * self.VOI / AIW
        EDI = self.EA / DI
        GC = self.OMEGA**2 * self.RI * np.sin(BI)
        VSL = Q / self.ZI / AI
        CS = (32.0 * np.cos(ANG)**2 + 16.0 * np.sin(ANG)**2) * DI
        CC = 1.25 - 0.5 * np.abs(np.sin(ANG))

        if CS > self.LI:
            CS = self.LI

        V24 = 19.0 * VSL / 6.0

        # guess a VSG
        VSG = VSL
        VM = VSL + VSG
        HLS = 1.0 / (1.0 + np.abs(VM / 8.66)**1.39)

        if HLS < 0.24:
            HLS = 0.24

        FE = 0.0
        HLF = VSL / VM
        VC = VSG / (1.0 - HLF)
        VCN = 0
        VF = VSL / HLF
        FI = FIC
        REMX = 5000.0
        RESL = self.DENL * VSL * DI / self.VISL
        ABCD = V24**2

        ABU = np.abs((VCN - VC) / VC)
        icon = 0

        while ABU > E1:
            if icon > 1000:
                break
            else:
                icon += 1

            # Entrainement fraction based on Oliemans et al's (1986)
            WEB = self.DENG * VSG * VSG * DI / self.ST
            FRO = np.sqrt(G * DI) / VSG
            RESG = self.DENG * VSG * DI / self.VISG
            CCC = 0.003 * WEB**1.8 * FRO**0.92 * RESL**0.7 * (self.DENL / self.DENG)**0.38 * (self.VISL / self.VISG)**0.97 / RESG**1.24
            FEN = CCC / (1.0 + CCC)

            if FEN > 0.9:
                FEN = 0.9
            FE = (FEN + 9.0 * FE) / 10.0

            # Translational velocity based on Nicklin (1962), Bendiksen (1984) and Zhang et al. (2000)
            if REMX < 2000.0:
                VAV = 2.0 * VM
            elif REMX > 2000.0:
                VAV = 1.3 * VM
            else:
                VAV = (2.0 - 0.7 * (REMX - 2000.0) / 2000.0) * VM

            VT = VAV + (0.54 * np.cos(ANG) + 0.35 * np.sin(ANG)) * np.sqrt(GC * DI * np.abs(self.DENL - self.DENG) / self.DENL)

            if VT < 0:
                VT = np.abs(VT)

            HLFN = ((HLS * (VT - VM) + VSL) * (VSG + VSL * FE) - VT * VSL * FE) / (VT * VSG)

            if HLFN < 0.0:
                HLFN = np.abs(HLFN)

            if HLFN > 1.0:
                HLFN = 1.0 / HLFN

            HLF = HLFN
            HLC = (1.0 - HLF) * VSL * FE / (VM - VSL * (1.0 - FE))

            if HLC < 0.0:
                HLC = 0.0

            # Taylor bubble geometries
            DeltaL = DI / 2.0 * (1. - np.sqrt(1.0 - HLF))
            AC = pi * (DI - 2.0 * DeltaL)**2 / 4.0
            AF = pi * DeltaL * (DI - 1.0 * DeltaL)
            SI = pi * (DI - 2.0 * DeltaL)
            SF = pi * DI
            DF = 4.0 * DeltaL * (DI - DeltaL) / DI
            DC = DI - 2.0 * DeltaL
            THF = DeltaL
            VFN = VSL * (1.0 - FE) / HLF
            VF = (VFN + 9.0 * VF) / 10.0

            # Reynolds number to get friction factor
            DENC = (self.DENL * HLC + self.DENG * (1.0 - HLF - HLC)) / (1.0 - HLF)
            REF = np.abs(self.DENL * VF * DF / self.VISL)
            REC1 = np.abs(self.DENG * VC * DC / self.VISG)
            FF = self.get_fff(REF, EDI)
            FIM = self.get_fff(REC1, THF / DI)

            # Interfacial friction factor based on Ambrosini et al. (1991)
            REG = np.abs(VC * self.DENG * DI / self.VISG)
            WED = self.DENG * VC * VC * DI / self.ST
            FS = 0.046 / REG**0.2
            SHI = np.abs(FI * self.DENG * (VC - VF)**2 / 2.0)
            THFO = THF * np.sqrt(np.abs(SHI * self.DENG)) / self.VISG
            FRA = FS * (1.0 + 13.8 * (THFO - 200.0 * np.sqrt(self.DENG / self.DENL)) * WED**0.2 / REG**0.6)
            FRA1 = self.get_fff(REC1, THF / DI)

            if FRA > FRA1:
                FRA = FRA1

            FIN = FRA

            if FIN > 0.5:
                FIN = 0.5

            FI = (FIN + 9.0 * FI) / 10.0
            ABCDN = (-(self.DENL - DENC) * GC + SF * FF * self.DENL * VF * VF / (2.0 * AF)) * 2 / (SI * FI * self.DENG * (1.0 / AF + 1.0 / AC))
            ABCD = (ABCDN + 9.0 * ABCD) / 10.0

            if ABCD < 0.:
                VCN = VC * 0.9
            else:
                VCN = np.sqrt(ABCD) + VF

            if VCN < V24:
                VCN = V24

            ABU = np.abs((VCN - VC) / VC)
            VC = VCN
            VSG = VC * (1.0 - HLF) - VSL * FE

            if VSG < 0.0:
                VSG = -VSG

            VM = VSL + VSG
            DPEX = (self.DENL * (VM - VF) * (VT - VF) * HLF + DENC * (VM - VC) * (VT - VC) * (1.0 - HLF)) * DI / CS / 4.0
            REMX = np.abs(DI * VM * self.DENL / self.VISL)
            FM = self.get_fff(REMX, EDI)
            DPSL = FM * self.DENL * VM * VM / 2.0
            DPAL = DPSL + DPEX

            if REMX < 5000.0:
                DPAL = DPAL * np.sqrt(REMX / 5000.0)

            AD = DPAL / (3.16 * CC * np.sqrt(self.ST * np.abs(self.DENL - self.DENG) * GC))
            HLSN = 1.0 / (1.0 + AD)

            if HLSN < 0.24:
                HLSN = 0.24

            HLS = HLSN
            LambdaC3 = VSG / (VSG + VSL)

        return LambdaC3

    def bubble_flow(self, HP, GF, Q, QLK):
        """
        Bubble flow model based on Zhu (2017)
        :param HP: single-phase pump pressure increment in psi
        :param GF: gas volumetric fraction at pump inlet
        :param Q: liquid flow rate in bpd
        :param QLK: leakage flow rate in bpd
        :return: GV (in-situ gas void fraction)
        """
        Q = Q * bbl_to_m3 / 24.0 / 3600.0       # bpd to m^3/s
        QLK = QLK * bbl_to_m3 / 24.0 / 3600.0
        HP = HP * psi_to_pa                     # psi to pa
        GV = GF
        ABV = 1.
        VSR = 0.1
        counter = 0

        while ABV > E1:
            counter += 1

            if counter > 10000:
                return GV

            DB = 6.034 * GF * (self.ST / self.DENL) ** 0.6 * \
                 (HP * (Q + QLK) / (self.DENL * self.ZI * self.VOI)) ** (-0.4) * (self.DENL / self.DENG) ** 0.2

            DBMAX = DB / 0.6
            REB = self.DENL * VSR * DBMAX / self.VISL
            SR = DBMAX * self.OMEGA / VSR

            if REB < 50.0:
                CD = 24.0 / REB * (1.0 + 0.15 * REB ** 0.687) * (1.0 + 0.3 * SR ** 2.5)
            else:
                CD = 24.0 / REB * (1.0 + 0.15 * REB ** 0.687) * (1.0 + 0.55 * SR ** 2.0)

            VSRN = np.sqrt(4.0 * DBMAX * (self.DENL - self.DENG) * self.RI / (3.0 * CD * self.DENL)) * self.OMEGA
            ABV = np.abs((VSRN - VSR) / VSR)
            VSR = VSRN
            RS = VSR * (2.0 * pi * self.RI - self.ZI * self.TB) * self.YI / (Q + QLK)

            if GV < 0.0:
                GV = 0.0
            else:
                GV = (RS - 1.0 + np.sqrt((1.0 - RS) ** 2 + 4.0 * RS * GF)) / (2.0 * RS)
        return GV

    def bubble_flow_sun(self, HP, GF, Q, QLK):
        """
        Use Sun (2003) drag force coefficient
        :param HP: single-phase pump pressure increment in psi
        :param GF: gas volumetric fraction at pump inlet
        :param Q: liquid flow rate in bpd
        :param QLK: leakage flow rate in bpd
        :return: GV (in-situ gas void fraction)
        """
        pass

    def bubble_flow_barrios(self, HP, GF, Q, QLK):
        """
        Use Barrios (2007) drag force coefficient, which is the original form Dr Zhang's development
        :param HP: single-phase pump pressure increment in psi
        :param GF: gas volumetric fraction at pump inlet
        :param Q: liquid flow rate in bpd
        :param QLK: leakage flow rate in bpd
        :return: GV (in-situ gas void fraction)
        """
        pass

    def intermittent_flow(self, Q, QLK, QG):
        ANG = -pi/2.0
        FIC = 0.0142
        Q = Q * bbl_to_m3 / 24.0 / 3600.0
        QLK = QLK * bbl_to_m3 / 24.0 / 3600.0
        QG = QG * bbl_to_m3 / 24.0 / 3600.0
        BI = (self.B1 + self.B2) / 2.0
        AI = self.VOI / self.LI
        AIW = self.AB + self.ASF + self.ASB
        DI = 4.0 * self.VOI / AIW
        EDI = self.EA / DI
        GC = self.OMEGA ** 2 * self.RI * np.sin(BI)

        VSL = (Q + QLK) / self.ZI / AI
        VSG = QG / self.ZI / AI
        VM = VSL + VSG
        HLS = 1.0 / (1.0 + np.abs(VM / 8.66) ** 1.39)
        GV = 0.25                                         # a first guess on gas void fraciton in slug unit

        # Translational velocity according to Nicklin (1962), Bendiksen (1984) and Zhang et al. (2000)
        REMX = np.abs(self.DENL * DI * VM / self.VISL)
        if REMX < 2000.0:
            VAV = 2.0 * VM
        elif REMX > 4000.0:
            VAV = 1.3 * VM
        else:
            VAV = (2.0 - 0.7 * (REMX - 2000.0) / 2000.0) * VM
        VT = VAV + (0.54 * np.cos(ANG) + 0.35 * np.sin(ANG)) * np.sqrt(GC * DI * np.abs(self.DENL - self.DENG) / self.DENL)

        if VT < 0:
            VT = -VT

        VMT = (VM + VT) / 2.0
        CS = 16 * DI
        CE = 1.25 - 0.5 * np.abs(np.sin(ANG))
        HLF = (VT - VM) / VT

        if HLF < 0:
            HLF = -HLF

        CU = CS * VM / VSL
        CF = CU - CS
        VF = VSL / 2.0
        VFN1 = VSL / 2.0
        VC = VM
        ABCD = 0.0
        FF = 0.1
        FI = FIC
        ABU = ABV = 1.
        counter = 0

        while (ABU > 0.001) or (ABV > 0.001):
            counter += 1
            if counter > 1500:
                return GV

            # Overall liquid holdup
            HL = (HLS * (VT - VM) + VSL) / VT
            HLF = HLS * (VT - VM) / (VT - VF)

            if HLF > HL:
                HLF = 0.99 * HL

            VCN = (VM - HLF * VF) / (1.0 - HLF)

            if VCN < 0.0:
                VCN = -VCN

            VC = VCN
            CFN = (HLS - HL) * CS / (HL - HLF)
            ABU = np.abs((CFN - CF) / CF)
            CF = (CFN + 9.0 * CF) / 10.0
            RSU = CS / (CS + CF)

            # Taylor bubble geometries
            DeltaL = DI / 2.0 * (1 - np.sqrt(1.0 - HLF))
            AC = pi * (DI - 2.0 * DeltaL) ** 2 / 4.0
            AF = pi * DeltaL * (DI - 1.0 * DeltaL)
            SI = pi * (DI - 2.0 * DeltaL)
            SF = pi * DI
            DF = 4.0 * DeltaL * (DI - DeltaL) / DI
            DC = DI - 2.0 * DeltaL
            THF = DeltaL

            # Slug liquid holdup
            DPEX = (self.DENL * (VM - VF) * (VT - VF) * HLF + self.DENG * (VM - VC) * (VT - VC) * (1.0 - HLF)) * DI / CS / 4.0
            REM = np.abs(self.DENL * VM * DI / self.VISL)
            FM = self.get_fff(REM, EDI)
            DPSL = FM * self.DENL * VM * VM / 2.0
            DPAL = DPSL + DPEX

            if REM < 5000.0:
                DPAL = DPAL * (REM / 5000.0)

            AD = DPAL / (3.16 * CE * np.sqrt(self.ST * np.abs(self.DENL - self.DENG) * GC))
            HLSN = 1.0 / (1.0 + AD)

            if HLSN < 0.24:
                HLSN=0.24

            HLS = HLSN

            if (REM < 1500.0) and (HLS < 0.6):
                HLS = 0.6

            if (VSL / VM > HLS) or (np.abs(CF) < DI):
                return GV

            # Reynolds numbers
            REF = np.abs(self.DENL * VF * DF / self.VISL)
            REC1 = np.abs(self.DENG * VC * DC / self.VISG)

            # Friction factors
            FFN = self.get_fff(REF, EDI)
            FIM = self.get_fff(REC1, THF / DI)
            FF = (FFN + 9.0 * FF) / 10.0

            # Interfacial friction factor (annular) according to Ambrosini et al. (1991)
            REG = np.abs(VC * self.DENG * DI / self.VISG)
            WED = self.DENG * VC * VC * DI / self.ST
            FS = 0.046 / REG ** 0.2
            SHI = np.abs(FI * self.DENG * (VC - VF) ** 2 / 2.0)
            THFO = THF * np.sqrt(np.abs(SHI * self.DENG)) / self.VISG
            FRA = FS * (1.0 + 13.8 * (THFO - 200.0 * np.sqrt(self.DENG / self.DENL)) * WED ** 0.2 / REG ** 0.6)
            FRA1 = self.get_fff(REC1, THF / DI)

            if FRA > FRA1:
                FRA=FRA1

            FIN = FRA

            if FIN > 0.5:
                FIN=0.5

            FI = (FIN + 9.0 * FI) / 10.0

            # Calculate film length CF using the combined momentum equation
            FSL = RSU * (self.DENL * (VM - VF) * (VT - VF) - self.DENG * (VM - VC) * (VT - VC)) / CF
            ABCDN = (FSL - SI * FI * self.DENG * (VC - VF) * np.abs(VC - VF) / 2.0 * (1.0 / AF + 1.0 / AC) + (self.DENL - self.DENG) * GC) * 2.0 * AF / (SF * FF * self.DENL)
            ABCD = (ABCDN + 19.0 * ABCD) / 20.0

            if ABCD > 0.0:
                VFN = np.sqrt(ABCD)
                if VFN > VM:
                    VFN = VM
            else:
                VFN = -np.sqrt(-ABCD)
                if VFN < - VM:
                    VFN = -VM

            ABV = abs((VFN - VF) / VF)
            VF = (VFN + 9.0 * VF) / 10.0
        return GV

    def segregated_flow(self, Q, QLK, QG):
        Q = Q * bbl_to_m3 / 24.0 / 3600.0
        QLK = QLK * bbl_to_m3 / 24.0 / 3600.0
        QG = QG * bbl_to_m3 / 24.0 / 3600.0
        FIC = 0.0142
        FI = FIC
        BI = (self.B1 + self.B2) / 2.0
        AI = self.VOI / self.LI
        AIW = self.AB + self.ASF + self.ASB
        DI = 4.0 * self.VOI / AIW
        EDI = self.EA / DI
        GC = self.OMEGA ** 2 * self.RI * np.sin(BI)

        VSL = (Q + QLK) / self.ZI / AI
        VSG = QG / self.ZI / AI
        VF = VC = 0.
        HLF = 1.0
        GV = 0.5

        while HLF > 0:
            HLF = HLF - 0.001

            # Entrainement fraction basd on Oliemans et al's (1986)
            WEB = self.DENG * VSG * VSG * DI / self.ST
            FRO = np.sqrt(G * DI) / VSG
            RESG = self.DENG * VSG * DI / self.VISG
            RESL = self.DENL * VSL * DI / self.VISL
            CCC = 0.003 * WEB ** 1.8 * FRO ** 0.92 * RESL ** 0.7 * (self.DENL / self.DENG) ** 0.38 * (self.VISL / self.VISG) ** 0.97 / RESG ** 1.24
            FEN = CCC / (1.0 + CCC)
            if FEN > 0.9:
                FEN = 0.9
            FE = FEN

            # Film geometries
            DeltaL = DI / 2.0 * (1 - np.sqrt(1.0 - HLF))
            AC = pi * (DI - 2.0 * DeltaL) ** 2 / 4.0
            AF = pi * DeltaL * (DI - 1.0 * DeltaL)
            SI = pi * (DI - 2.0 * DeltaL)
            SF = pi * DI
            DF = 4.0 * DeltaL * (DI - DeltaL) / DI
            DC = DI - 2.0 * DeltaL
            THF = DeltaL

            VFN = VSL * (1.0 - FE) / HLF
            VF = (VFN + 9.0 * VF) / 10.0
            VCN = (VSG + VSL * FE) / (1 - HLF)
            VC = (VCN + 9.0 * VC) / 10.0

            ALPHAC = VSG / (VSG + VSL * FE)
            DENC = self.DENG * ALPHAC + self.DENL * (1.0 - ALPHAC)
            VISC = self.VISG * ALPHAC + self.VISL * (1.0 - ALPHAC)
            REF = np.abs(self.DENL * VF * DF / self.VISL)
            REC1 = np.abs(self.DENG * VC * DC / self.VISC)

            FF = self.get_fff(REF, EDI)
            FIM = self.get_fff(REC1, THF/DI)

            # Interfacial friction factor (annular) according to Ambrosini et al. (1991)
            REG = np.abs(VC * self.DENG * DI / self.VISG)
            WED = self.DENG * VC * VC * DI / self.ST
            FS = 0.046 / REG ** 0.2
            SHI = np.abs(FI * self.DENG * (VC - VF) ** 2 / 2.0)
            THFO = THF * np.sqrt(np.abs(SHI * self.DENG)) / self.VISG
            FRA = FS * (1.0 + 13.8 * (THFO - 200.0 * np.sqrt(self.DENG / self.DENL)) * WED ** 0.2 / REG ** 0.6)
            FRA1 = self.get_fff(REC1, THF / DI)

            if FRA > FRA1:
                FRA=FRA1

            FIN = FRA

            if FIN > 0.5:
                FIN=0.5

            FI = (FIN + 9.0 * FI) / 10.0
            ABCD = (self.DENL - DENC) * GC - SF * FF * self.DENL * VF * VF / (2.0 * AF) + SI * FI * DENC * (VC - VF) ** 2.0 * (1.0 / 2.0 / AF + 1.0 / 2.0 / AC)

            if np.abs(ABCD) < 0.001:
                GV = ALPHAC * (1.0 - 2.0 * DeltaL / DI) ** 2.0
                break
        return GV

    def gl_calculate(self, qg, ql):
        ql = ql * bbl_to_m3 / 24.0 / 3600.0
        qg = qg * bbl_to_m3 / 24.0 / 3600.0
        GF = ql / (ql + qg)
        pass

    def surging_curve(self):
        pass

    def mapping_curve(self):
        pass

    def flow_pattern(self):
        QL = np.arange(500)[1:] * 50.0
        QSG1 = []
        QSG2 = []
        QSG3 = []
        for ql in QL:
            # ql in bpd
            HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_old(ql)
            if HP > 0:
                LambdaC1 = self.get_lambda_c1(ql + QLK, HP)
                LambdaC2 = self.get_lambda_c2(ql + QLK, HP)
                LambdaC3 = self.get_lambda_c3(ql + QLK)

                if LambdaC1 > LambdaC2: LambdaC2 = LambdaC1

                QSG1.append(LambdaC1 / (1.0 - LambdaC1) * ql)
                QSG1.append(LambdaC2 / (1.0 - LambdaC2) * ql)
                QSG1.append(LambdaC3 / (1.0 - LambdaC3) * ql)
            else:
                break
        return QL[:len(QSG1)], QSG1, QSG2, QSG3


########################################################################################################################
if __name__ == "__main__":

    TE2700 = {"R1": 0.017496, "R2": 0.056054, "TB": 0.00272, "TV": 0.00448, "RD1": 0.0219507825, "RD2": 0.054735,
              "YI1": 0.012194, "YI2": 0.007835, "VOI": 0.000016119, "VOD": 0.000011153, "ASF": 0.00176464,
              "ASB": 0.00157452, "AB": 0.001319, "AV": 0.001516, "ADF": 0.001482, "ADB": 0.000935, "LI": 0.076,
              "LD": 0.08708, "RLK": 0.056209, "LG": 0.00806, "SL": 0.00005, "EA": 0.0000254, "ZI": 5, "ZD": 9,
              "B1": 19.5, "B2": 24.7, "NS": 1600, "DENL": 1000, "DENG": 11.2, "VISL": 0.001, "VISG": 0.000018,
              "VISW": 0.001, "ST": 0.073, "N": 3500, "SGM": 0.3, "QL": 2700, "QG": 50}

    # sgl = SinglePhaseModel(TE2700, 4800)
    # QL, hpsgl, hesgl, hfsgl, htsgl, hresgl, qlksgl = sgl.PerformanceCurve()

    # figure1 = plt.figure()
    # plt.plot(QL, hpsgl)
    # plt.show()

    gl = GasLiquidModel(TE2700, 5200)
    HP, HE, HF, HT, HD, HRE, HLK, QLK = gl.sgl_calculate_new(1000)
