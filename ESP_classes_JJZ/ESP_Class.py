# -*- coding: utf-8 -*-

"""
The ESP calculation module contains two classes for computing inputs performance under single-phase water/viscous
fluid flow or gas-liquid two-phase flow conditions.

The two-phase model was originally proposed by Dr. Zhang, TUALP ABM (2013) and later revised by Zhu (2017). Only the 
simplified version is programmed.

Version:    1st, Aug, 2017
Developer:  Jianjun Zhu
"""

import numpy as np
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
import matplotlib as mpl
from losses import FrictionLoss
from ESP_inputs import ESP

# global variables
G = 9.81
pi = np.pi
E1 = 1e-5
DENW = 997.                 # water density
VISW = 1e-3                 # water viscosity
psi_to_pa = 1.013e5 / 14.7
psi_to_ft = 2.3066587368787
bbl_to_m3 = 0.15897
colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
symbols = ['o', 's', '^', '*', 'd', 'p', 'v', 'D', '<', '>']

# customize matplotlib
# mpl.rcParams['xtick.labelsize'] = 18
# mpl.rcParams['ytick.labelsize'] = 18
# mpl.rcParams['font.size'] = 18
# mpl.rcParams['xtick.direction'] = 'in'
# mpl.rcParams['ytick.direction'] = 'in'
# mpl.rcParams['xtick.top'] = 'True'
# mpl.rcParams['ytick.right'] = 'True'
# mpl.rcParams['figure.figsize'] = (5, 5)
#
plt.style.use('seaborn-ticks')
mpl.rcParams['figure.figsize'] = (4, 3)
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['markers.fillstyle'] = 'none'
mpl.rcParams['lines.markersize'] = 5


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
    conn.commit()
    conn.close()


def gasdensity(p, t, h):
    """
    gas density based on CIPM-81 (Davis, 1972) correlations
    :param p: pressure in psig
    :param t: temperature in Fahrenheit
    :param h: humidity in %
    :return: gas density in kg/m3
    """
    A = 1.2811805e-5
    B = -1.9509874e-2
    C = 34.04926034
    D = -6.3536311e3
    alpha = 1.00062
    beta = 3.14e-8
    gamma = 5.6e-7
    a0 = 1.62419e-6
    a1 = -2.8969e-8
    a2 = 1.0880e-10
    b0 = 5.757e-6
    b1 = -2.589e-8
    c0 = 1.9297e-4
    c1 = -2.285e-6
    d = 1.73e-11
    e = -1.034e-8
    R = 8.31441
    Ma = 28.9635  # air molecular weight, g/mol
    Mv = 18  # water molecular weight, g/mol

    Pabs = (p + 14.7) * 6894.76
    Tt = (t - 32) / 1.8
    Tabs = Tt + 273.15
    psv = 1.0 * np.exp(A * (Tabs) ** 2.0 + B * Tabs + C + D / Tabs)
    f = alpha + beta * Pabs + gamma * (Tt) ** 2
    xv = h / 100.0 * f * psv / Pabs
    Z = 1.0 - Pabs / Tabs * (a0 + a1 * Tt + a2 * (Tt) ** 2.0 + (b0 + b1 * Tt) * xv + (c0 + c1 * Tt) * (xv) ** 2) + \
        (Pabs / Tabs) ** 2.0 * (d + e * (xv) ** 2.0)
    return Pabs * Ma / 1000.0 / (Z * R * Tabs) * (1.0 - xv * (1.0 - Mv / Ma))


#################################
# single-phase mechanistic model
class SinglePhaseModel(object):
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
        self.SGM = inputs['SGM']
        self.ST = inputs['ST']
        self.NS = inputs['NS']
        self.WC = inputs['WC']
        self.SN = inputs['SN']

        self.QBEM = QBEM * bbl_to_m3 / 24.0 / 3600.0
        self.OMEGA = 2.0 * pi * self.N / 60.0

        # use literature loss models
        self.friction_loss = FrictionLoss(inputs, self.N, self.DENL, self.VISL)

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
    def emulsion(VOI, R2, VISO, VISW, DENO, DENW, WC, ST, N, Q, SN, mod='tualp'):
        """
        The model is based on Brinkman (1952) correlation and Zhang (2017, Fall)
        :param VOI: volume of impeller (m3)
        :param R2:  impeller outer radius (m)
        :param VISO:viscosity of oil (kg/m-s)
        :param VISW:viscosity of water (kg/m-s)
        :param DENO:density of oil (kg/m3)
        :param DENW:density of water (kg/m3)
        :param WC:  water cut (%)
        :param ST:  surface tension (N/m)
        :param N:   rotational speed (rpm)
        :param Q:   flow rate (m3/s)
        :param SN:  stage number (-)
        :param mod: select different model type: tualp, banjar, or zhu
        :return: miu in Pas
        """
        E = 3.  # exponential index
        f = N / 60.
        WC = WC / 100.
        miu_tilda = VISO / VISW
        phi_OI = miu_tilda ** (1. / E) / (1. + miu_tilda ** (1. / E))
        phi_WI = 1. - phi_OI
        phi_OE = 1. - (VISW / VISO) ** (1. / E)
        i = 0.

        get_C = lambda SN, WE, RE, ST: SN ** 0.01 * WE ** 0.1 * RE ** 0.1 / (2.5 * ST ** 0.2)

        if mod == "tualp":
            get_C = lambda SN, WE, RE, ST: (SN * WE * RE) ** 0.15 / (10 * ST ** 0.5)
        elif mod == "banjar":
            get_C = lambda SN, WE, RE, ST: (SN * WE * RE) ** 0.1 / (10 * ST ** 0.2)
        elif mod == "zhu":
            # get_C = lambda SN, WE, RE, ST: (SN * WE) ** 0.1 * RE ** 0.1 / (2.5 * ST ** 0.2)
            get_C = lambda SN, WE, RE, ST: SN ** 0.01 * WE ** 0.1 * RE ** 0.1 / (2.5 * ST ** 0.2)

        # find the inversion point
        St = f * VOI / Q
        for i in np.arange(1000) / 1000.:
            if i == 0:
                continue
            rouA = i * DENW + (1 - i) * DENO
            We = rouA * Q ** 2 / (ST * VOI)
            miu_M = VISW / (1 - (1 - i) * phi_OE) ** E

            # assume oil in water
            Re = rouA * Q / (VISW * 2 * R2)
            C = get_C(SN, We, Re, St)
            miu_E = VISW / (1 - (1 - i)) ** E
            miu_A_OW = C * (miu_E - miu_M) + miu_M

            # assume water in oil
            Re = rouA * Q / (VISO * 2 * R2)
            C = get_C(SN, We, Re, St)
            miu_E = VISO / (1 - i) ** E
            miu_A_WO = C * (miu_E - miu_M) + miu_M

            if np.abs(miu_A_OW - miu_A_WO) / miu_A_OW < 0.01:
                break

        if WC > i:
            # oil in water
            rouA = WC * DENW + (1 - WC) * DENO
            We = rouA * Q ** 2 / (ST * VOI)
            miu_M = VISW / (1 - (1 - WC) * phi_OE) ** E

            Re = rouA * Q / (VISW * 2 * R2)
            C = get_C(SN, We, Re, St)
            miu_E = VISW / (1 - (1 - WC)) ** E
            miu_A = C * (miu_E - miu_M) + miu_M
        else:
            # water in oil
            rouA = WC * DENW + (1 - WC) * DENO
            We = rouA * Q ** 2 / (ST * VOI)
            miu_M = VISW / (1 - (1 - WC) * phi_OE) ** E

            Re = rouA * Q / (VISO * 2 * R2)
            C = get_C(SN, We, Re, St)
            miu_E = VISO / (1 - WC) ** E
            miu_A = C * (miu_E - miu_M) + miu_M
        return miu_A

    def sgl_calculate_old(self, Q):
        """
        Original development on single-phase model by Zhu and Zhang
        :param Q: flow rate in bpd
        :return: HP, HE, HF, HT, HD, HRE, HLK, QLK in filed units
        """
        # Q in bpd
        Q = Q * bbl_to_m3 / 24.0 / 3600.0
        QLK = 0.02 * Q
        HP = 10.
        HE, HEE = 0., 0
        HFI, HFD = 0., 0
        HTI, HTD = 0., 0
        HLK = 0.
        icon = 0
        HP_new, QLK_new = 1., 1.

        ABH = np.abs((HP - HP_new) / HP_new)
        ABQ = np.abs((QLK - QLK_new) / QLK_new)
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

    def sgl_calculate_new(self, Q, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC):
        """
        Dr Zhang new update on previous single-phase model with new consideration on recirculation flow loss
        :param Q:   flow rate in m3/s
        :param QEM: best match flow rate in m3/s
        :param DENL: liquid density in kg/m3
        :param DENW: water density in kg/m3
        :param N:   rotational speed in rpm
        :param NS:  specific speed based on field units
        :param SGM: tuning factor
        :param SN:  stage number
        :param ST:  surface tension Nm
        :param VISL: liquid viscosity in Pas
        :param VISW: water viscosity in Pas
        :param WC:  water cut in %
        :return: HP, HE, HF, HT, HD, HRE, HLK, QLK in filed units
        """

        QLK = 0.02 * Q
        HP = 10.
        HE, HEE = 0., 0
        HFI, HFD = 0., 0
        HTI, HTD = 0., 0
        HLK = 0.
        icon = 0
        HP_new, QLK_new = 1., 1.

        ABH = np.abs((HP - HP_new) / HP_new)
        ABQ = np.abs((QLK - QLK_new) / QLK_new)
        AIW = self.AB + self.ASF + self.ASB
        ADW = self.AV + self.ADF + self.ADB

        # needs change
        SGM = self.SGM
        OMEGA = 2.0 * pi * N / 60.0

        # new QBEM due to rotational speed
        QBEM = QBEM * (N / 3500)

        # check if emulsion occurs
        if WC > 0.:
            VISL = self.emulsion_zhang(self.VOI, self.R2, VISL, VISW, DENL, DENW, WC, ST, N, Q, 3)

        while (ABH > E1) or (ABQ > E1):

            C1M = (Q + QLK) / ((2.0 * pi * self.R1 - self.ZI * self.TB) * self.YI1)
            C2M = (Q + QLK) / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            U1 = self.R1 * OMEGA
            U2 = self.R2 * OMEGA
            W1 = C1M / np.sin(self.B1)
            W2 = C2M / np.sin(self.B2)
            C1 = np.sqrt(C1M ** 2 + (U1 - C1M / np.tan(self.B1)) ** 2)
            C2 = np.sqrt(C2M ** 2 + (U2 - C2M / np.tan(self.B2)) ** 2)
            CMB = QBEM / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            C2B = np.sqrt(CMB ** 2 + (U2 - CMB / np.tan(self.B2)) ** 2)

            # Euler head
            # HE=(U2**2-U1**2+W1**2-W2**2+C2**2-C1**2)/(2.0*G)					        # with pre-rotation
            HE = (U2 ** 2 - U2 * C2M / np.tan(self.B2)) / G                             # without pre-rotation

            # head loss due to recirculation
            if (Q + QLK) < QBEM:
                VSH = U2 * (QBEM - (Q + QLK)) / QBEM
                C2F = C2B * (Q + QLK) / QBEM
                DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
                REC = DENL * VSH * DC / VISL
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F
                HEE = HE + (C2E ** 2 - C2 ** 2) / (2.0 * G)
            else:
                VSH = U2 * (Q + QLK - QBEM) / QBEM
                C2F = C2B * (Q + QLK) / QBEM
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F) * (Q + QLK - QBEM) / QBEM     # new development by Dr Zhang
                HEE = HE + (C2E ** 2 - C2 ** 2) / (2.0 * G)

            # friction loss
            AI = self.VOI / self.LI
            AD = self.VOD / self.LD
            VI = (Q + QLK) / self.ZI / AI
            VD = Q / self.ZD / AD
            DI = 4.0 * self.VOI / AIW
            DD = 4.0 * self.VOD / ADW
            REI = DENL * (W1 + W2) * DI / VISL / 2.0
            RED = DENL * VD * DD / VISL
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
            UL = self.RLK * OMEGA
            HIO = HEE - HFI - HTI
            HLK = HIO - (U2 ** 2 - UL ** 2) / (8.0 * G)
            if HLK >= 0:
                VL = np.abs(QLK) / (2.0 * pi * self.RLK * self.SL)
                REL = DENL * VL * self.SL / VISL
                EDL = 0.0
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * G * HLK / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLK_new = 2.0 * pi * self.RLK * self.SL * VL
            else:
                VL = np.abs(QLK / (2.0 * pi * self.RLK * self.SL))
                REL = DENL * VL * self.SL / VISL
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
        HP = HP * G * DENL / psi_to_pa
        HE = HE * G * DENL / psi_to_pa
        HEE = HEE * G * DENL / psi_to_pa
        HF = (HFI + HFD) * G * DENL / psi_to_pa
        HT = (HTI + HTD) * G * DENL / psi_to_pa
        HD = (HFD + HTD) * G * DENL / psi_to_pa
        HRE = np.abs(HE - HEE)
        HLK = HLK * G * DENL / psi_to_pa
        QLK = QLK * 24.0 * 3600.0 / bbl_to_m3
        return HP, HE, HF, HT, HD, HRE, HLK, QLK

    def sgl_calculate_jc(self, Q, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, W):
        """
        Jiecheng Zhang (2017) thesis development
        :param Q:   flow rate in m3/s
        :param QEM: best match flow rate in m3/s
        :param DENL: liquid density in kg/m3
        :param DENW: water density in kg/m3
        :param N:   rotational speed in rpm
        :param NS:  specific speed based on field units
        :param SGM: tuning factor
        :param SN:  stage number
        :param ST:  surface tension Nm
        :param VISL: liquid viscosity in Pas
        :param VISM: water viscosity in Pas
        :param WC:  water cut in %
        :return: HP, HE, HF, HT, HD, HRE, HLK, QLK in filed units
        """
        QLK = 0.02 * Q
        HP = 10.
        HE, HEE = 0., 0
        HFI, HFD = 0., 0
        HTI, HTD = 0., 0
        HLK = 0.

        AIW = self.AB + self.ASF + self.ASB
        ADW = self.AV + self.ADF + self.ADB

        # slip factor compared to Wiesner (1961)
        SGMU = 1 - np.sqrt(np.sin(self.B2)) / self.ZI ** (1.5 * (3448 / NS) ** 0.4)
        OMEGA = 2.0 * pi * N / 60.0

        # new QBEM due to liquid viscosity
        QBEM = QBEM * (N / 3600) * (VISL / VISW) ** (0.01 * (3448 / NS))

        for icon in range(1, 100):
            C1M = (Q + QLK) / ((2.0 * pi * self.R1 - self.ZI * self.TB) * self.YI1)
            C2M = (Q + QLK) / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            U1 = self.R1 * OMEGA
            U2 = self.R2 * OMEGA
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
                REC = DENL * VSH * DC / VISL
                SGM = (VISW / VISL) ** 0.1 / (10.0 + 0.02 * REC ** 0.2)   # SGM: shear factor due to viscosity
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F)
                HEE = HE + (C2E ** 2 - C2 ** 2) / (2.0 * G)
            else:
                VSH = U2 * (Q + QLK - QBEM) / QBEM
                C2F = C2B * (Q + QLK) / QBEM
                DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
                REC = DENL * VSH * DC / VISL
                SGM = (VISW / VISL) ** 0.1 / (10.0 + 0.02 * REC ** 0.2)
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F) * (Q + QLK - QBEM) / QBEM
                HEE = HE + (C2E ** 2 - C2 ** 2) / (2.0 * G)

            AI = self.VOI / self.LI
            AD = self.VOD / self.LD
            VI = (Q + QLK) / self.ZI / AI
            VD = Q / self.ZD / AD
            DI = 4.0 * self.VOI / AIW
            DD = 4.0 * self.VOD / ADW
            REI = DENL * VI * DI / VISL
            RED = DENL * VD * DD / VISL / 1.0
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

            UL = self.RLK * OMEGA
            HIO = HEE - HFI - HTI
            HLK = HIO - (U2 ** 2 - UL ** 2) / (8.0 * G)

            if HLK >= 0.:
                VL = abs(QLK) / (2.0 * pi * self.RLK * self.SL)
                REL = DENL * VL * self.SL / VISL
                EDL = 0.0
                FFL = self.get_fff_jc(REL, EDL)
                VL = np.sqrt(2.0 * G * HLK / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLKN = 2.0 * pi * self.RLK * self.SL * VL
            else:
                VL = abs(QLK / (2.0 * pi * self.RLK * self.SL))
                REL = DENL * VL * self.SL / VISL
                EDL = 0.0
                FFL = self.get_fff_jc(REL, EDL)
                VL = np.sqrt(2.0 * G * abs(HLK) / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLKN = -2.0 * pi * self.RLK * self.SL * VL

            ABL = abs((QLKN - QLK) / QLKN)
            QLK = QLKN
            ABC = abs((HPN - HP) / HPN)
            HP = HPN

            if (ABC < E1) and (ABL < E1): break

        # return pressure in psi, flow rate in bpd
        HP = HP * G * DENL / psi_to_pa
        HE = HE * G * DENL / psi_to_pa
        HEE = HEE * G * DENL / psi_to_pa
        HF = (HFI + HFD) * G * DENL / psi_to_pa
        HT = (HTI + HTD) * G * DENL / psi_to_pa
        HD = (HFD + HTD) * G * DENL / psi_to_pa
        HRE = np.abs(HE - HEE)
        HLK = HLK * G * DENL / psi_to_pa
        QLK = QLK * 24.0 * 3600.0 / bbl_to_m3
        return HP, HE, HF, HT, HD, HRE, HLK, QLK

    def sgl_calculate_2018(self, Q, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC):
        """
        Based on previous single-phase model (Jiecheng Zhang) with new development on pressure losses
        :param Q:   flow rate in m3/s
        :param QEM: best match flow rate in m3/s
        :param DENL: liquid density in kg/m3
        :param DENW: water density in kg/m3
        :param N:   rotational speed in rpm
        :param NS:  specific speed based on field units
        :param SGM: tuning factor
        :param SN:  stage number
        :param ST:  surface tension Nm
        :param VISL: liquid viscosity in Pas
        :param VISM: water viscosity in Pas
        :param WC:  water cut in %
        :return: HP, HE, HF, HT, HD, HRE, HLK, QLK in filed units
        """
        QLK = 0.02 * Q
        HP = 10.
        HE, HEE = 0., 0
        HFI, HFD = 0., 0
        HTI, HTD = 0., 0
        HLK = 0.
        icon = 0
        HP_new, QLK_new = 1., 1.
        OMEGA = 2.0 * pi * N / 60.0

        ABH = np.abs((HP - HP_new) / HP_new)
        ABQ = np.abs((QLK - QLK_new) / HP_new)
        AIW = self.AB + self.ASF + self.ASB
        ADW = self.AV + self.ADF + self.ADB

        # slip factor compared to Wiesner (1961)
        SGMU = 1 - np.sqrt(np.sin(self.B2)) / self.ZI ** (1.5 * (3448 / NS) ** 0.4)

        # new QBEM due to liquid viscosity
        QBEM = QBEM * (N / 3600) * (VISL / VISW) ** (0.01 * (3448 / NS)**4)

        # check if emulsion occurs
        if WC > 0.:
            VISL = self.emulsion(self.VOI, self.R2, VISL, VISW, DENL, DENW, WC, ST, N, Q, SN, "zhu")

        while (ABH > E1) or (ABQ > E1):

            C1M = (Q + QLK) / ((2.0 * pi * self.R1 - self.ZI * self.TB) * self.YI1)
            C2M = (Q + QLK) / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            U1 = self.R1 * OMEGA
            U2 = self.R2 * OMEGA
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
                REC = DENL * VSH * DC / VISL
                
                # SGM: shear factor due to viscosity
                SGM = (VISW / VISL) ** 0.1 / (10.0 + 0.02 * REC ** 0.25)
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F)
                HEE = HE + (C2E ** 2 - C2 ** 2) / (2.0 * G)
            else:
                VSH = U2 * (Q + QLK - QBEM) / QBEM
                C2F = C2B * (Q + QLK) / QBEM
                DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
                REC = DENL * VSH * DC / VISL
                SGM = (VISW / VISL) ** 0.1 / (10.0 + 0.02 * REC ** 0.2)
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
            REI = DENL * (W1 + W2) * DI / VISL / 2.0
            RED = DENL * VD * DD / VISL
            EDI = self.EA / DI
            EDD = self.EA / DD
            FFI = self.get_fff(REI, EDI)
            FFD = self.get_fff(RED, EDD)
            # HFI = 2.5 * 4.0 * FFI * (W1 + W2) ** 2 * self.LI / (8.0 * G * DI)
            HFI = self.friction_loss.sun2003(Q+QLK)                         # use Sun(2003) friction model
            HFD = 2.5 * 4.0 * FFD * VD ** 2 * self.LD / (2.0 * G * DD)

            # turn loss
            FTI = 3.0
            FTD = 3.0
            HTI = FTI * VI ** 2 / (2.0 * G)
            HTD = FTD * VD ** 2 / (2.0 * G)

            # new pump head
            HP_new = HEE - HFI - HFD - HTI - HTD

            # calculate leakage
            UL = self.RLK * OMEGA
            HIO = HEE - HFI - HTI
            HLK = HIO - (U2 ** 2 - UL ** 2) / (8.0 * G)
            if HLK >= 0:
                VL = np.abs(QLK) / (2.0 * pi * self.RLK * self.SL)
                REL = DENL * VL * self.SL / VISL
                EDL = 0.0
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * G * HLK / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLK_new = 2.0 * pi * self.RLK * self.SL * VL
            else:
                VL = np.abs(QLK / (2.0 * pi * self.RLK * self.SL))
                REL = DENL * VL * self.SL / VISL
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
        HP = HP * G * DENL / psi_to_pa
        HE = HE * G * DENL / psi_to_pa
        HEE = HEE * G * DENL / psi_to_pa
        HF = (HFI + HFD) * G * DENL / psi_to_pa
        HT = (HTI + HTD) * G * DENL / psi_to_pa
        HD = (HFD + HTD) * G * DENL / psi_to_pa
        HRE = np.abs(HE - HEE)
        HLK = HLK * G * DENL / psi_to_pa
        QLK = QLK * 24.0 * 3600.0 / bbl_to_m3
        return HP, HE, HF, HT, HD, HRE, HLK, QLK

    def performance_curve(self):
        QL = np.arange(5000)[1:] * 50.0
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
            ql = ql * bbl_to_m3 / 24.0 / 3600.0
            HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_2018(ql, self.QBEM, self.DENL, DENW, self.N, self.NS,
                                                                        self.SGM, self.SN, self.ST, self.VISL, VISW,
                                                                        self.WC)

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


########################################
# two-phase mechanistic model class
class GasLiquidModel(SinglePhaseModel):
    def __init__(self, inputs, QBEM):
        super(GasLiquidModel, self).__init__(inputs, QBEM)
        self.DENG = inputs['DENG']
        self.VISG = inputs['VISG']
        self.QL = inputs['QL']
        self.QG = inputs['QG']
        self.RI = (self.R1 + self.R2) / 2.0
        self.YI = (self.YI1 + self.YI2) / 2.0

    def get_lambda_c1(self, Q, HP):
        # critical bubble size in turbulent flow (Barnea 1982)
        DBCrit = 2.0 * (0.4 * 0.073 / (self.DENL - self.DENG) / (self.OMEGA**2 * self.RI))**0.5
        LambdaC1 = DBCrit / (6.034 / 0.6 * (self.ST / self.DENL)**0.6 * (HP * Q / (self.DENL * self.ZI * self.VOI)) **
                             (-0.4) * (self.DENL / self.DENG)**0.2)
        return LambdaC1

    def get_lambda_c2(self, Q, HP):
        alphaG_crit = pi / 6.0 - (pi / 6.0 - 0.25) * np.exp(-self.N / 3600.0)
        VSR = 0.1
        LambdaC2 = 1.0
        ABV = 1.0
        icon = 0

        while ABV > 0.001:
            LambdaC2 -= 1e-5
            DB = 6.034 * LambdaC2 * (self.ST / self.DENL)**0.6 * (HP * Q / (self.DENL * self.ZI * self.VOI))**(-0.4) * \
                 (self.DENL / self.DENG)**0.2
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

    @staticmethod
    def gv_zhu(DENL, DENG, GF, HP, N, Q, QLK, RI, ST, TB, VISL, VOI, YI, ZI):
        """
        Bubble flow model based on Zhu (2017)
        :param DENL: liquid density in m3/s
        :param DENG: gas density in m3/s
        :param GF: gas volumetric fraction at pump inlet
        :param HP: single-phase pump pressure increment in pa
        :param N: rotatonal speed in rpm
        :param Q: liquid flow rate in m3/s
        :param QLK: leakage flow rate in m3/s
        :param RI: radius of impeller in m
        :param ST: interfacial tension in N/m
        :param TB: impeller blade thickness in m
        :param VISL: liquid viscosity in Pas
        :param VOI: impeller volume in m3
        :param YI: impeller outlet height
        :param ZI: impeller blade number
        :return: GV (in-situ gas void fraction), CD_to_rb(CD over rb)
        """
        GV, DBMAX = GF, 0.3
        CD = 0.
        ABV = 1.
        VSR = 0.1
        counter = 0
        OMEGA = 2.0 * pi * N / 60.0

        while ABV > E1:
            counter += 1
            if counter > 10000:
                return GV, CD/DBMAX

            DB = 6.034 * GF * (ST / DENL) ** 0.6 * (HP * (Q + QLK) / (DENL * ZI * VOI)) ** (-0.4) * (DENL / DENG) ** 0.2
            DBMAX = DB / 0.4
            REB = DENL * VSR * DBMAX / VISL
            SR = DBMAX * OMEGA / VSR

            if REB < 50.0:
                CD = 24.0 / REB * (1.0 + 0.15 * REB ** 0.687) * (1.0 + 0.3 * SR ** 2.5)
            else:
                CD = 24.0 / REB * (1.0 + 0.15 * REB ** 0.687) * (1.0 + 0.55 * SR ** 2.0)

            VSRN = np.sqrt(4.0 * DBMAX * (DENL - DENG) * RI / (3.0 * CD * DENL)) * OMEGA
            ABV = np.abs((VSRN - VSR) / VSR)
            VSR = VSRN
            RS = VSR * (2.0 * pi * RI - ZI * TB) * YI / (Q + QLK)

            if GV < 0.0:
                GV = 0.0
            else:
                GV = (RS - 1.0 + np.sqrt((1.0 - RS) ** 2 + 4.0 * RS * GF)) / (2.0 * RS)

        return GV, CD/DBMAX

    @staticmethod
    def gv_sun(B1, B2, DENL, DENG, N, Q, QG, QLK, R1, R2, TB, VISL, VISG, YI1, YI2, YI, ZI):
        """
        Use Sun (2003) drag force coefficient
        :param B1: impeller inlet blade angle in radian
        :param B2: impeller outlet blade angle in radian
        :param DENL: liuqd density in kg/m3
        :param DENG: gas density in kg/m3
        :param GF: gas volumetric fraction at pump inlet
        :param N: rotational speed in rpm
        :param Q: liquid flow rate in m3/s
        :param QLK: leakage flow rate in m3/s
        :param QG: gas flow rate in m3/s
        :param R2: impeller outlet radius in m
        :param RI: impeller representative radius in m
        :param TB: impeller blade thickness in m
        :param VISL: liquid viscosity in Pas
        :param VISG: gas viscosity in Pas
        :param YI1: impeller inlet height in m
        :param YI2: impeller outlet height in m
        :param YI: impeller representative height in m
        :param ZI: impeller balde number
        :return: GV (in-situ gas void fraction), CD_to_rb(CD over rb)
        """
        GF = QG / (Q + QLK + QG)
        GV1, CD_to_rb = 1-GF, 0.3
        GV2 = GV1
        GV = (GV1 + GV2) / 2
        ABV1, ABV2 = 1., 1.
        VSR1, VSR2 = 0.1, 0.1
        CD_to_rb = 0.1
        counter = 0
        OMEGA = 2.0 * pi * N / 60.0

        C1M_L = (Q + QLK) / ((2.0 * pi * R1 - ZI * TB) * YI1)
        C2M_L = (Q + QLK) / ((2.0 * pi * R2 - ZI * TB) * YI2)
        W1_L, W2_L = C1M_L / np.sin(B1), C2M_L / np.sin(B2)

        C1M_G = QG / ((2.0 * pi * R1 - ZI * TB) * YI1)
        C2M_G = QG / ((2.0 * pi * R2 - ZI * TB) * YI2)
        W1_G, W2_G = C1M_G / np.sin(B1), C2M_G / np.sin(B2)

        while (ABV1 > E1) or (ABV2 > E1):
            counter += 1
            if counter > 500:
                return GV, CD_to_rb

            W1L = W1_L / (1 - GV1)
            W2L = W2_L / (1 - GV2)
            W1G = W1_G / GV1
            W2G = W2_G / GV2
            mium1 = (1 - GV1) * VISL + GV1 * VISG
            DENM1 = (1 - GV1) * DENL + GV1 * DENG
            mium2 = (1 - GV2) * VISL + GV2 * VISG
            DENM2 = (1 - GV2) * DENL + GV2 * DENG

            if GV1 <= 0.15:
                CD_to_rb1 = 12. * mium1 * (4.564e7 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)
            elif (GV1 > 0.15) and (GV1 < 0.5):
                CD_to_rb1 = 12. * mium1 * (6.39e7 * (1 - GV1) ** 1.5 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)
            else:
                CD_to_rb1 = 12. * mium1 * (9.13e7 * GV1 ** 2 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)

            VSRN1 = np.sqrt(4.0 * (DENL - DENG) * R1 / (3.0 * CD_to_rb1 * DENL)) * OMEGA
            ABV1 = np.abs((VSRN1 - VSR1) / VSR1)
            VSR1 = VSRN1
            RS1 = VSR1 * (2.0 * pi * R1 - ZI * TB) * YI / (Q + QLK)

            if GV2 <= 0.15:
                CD_to_rb2 = 12. * mium2 * (4.564e7 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)
            elif (GV2 > 0.15) and (GV2 < 0.5):
                CD_to_rb2 = 12. * mium2 * (6.39e7 * (1 - GV2) ** 1.5 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)
            else:
                CD_to_rb2 = 12. * mium2 * (9.13e7 * GV2 ** 2 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)

            VSRN2 = np.sqrt(4.0 * (DENL - DENG) * R2 / (3.0 * CD_to_rb2 * DENL)) * OMEGA
            ABV2 = np.abs((VSRN2 - VSR2) / VSR2)
            VSR2 = VSRN2
            RS2 = VSR2 * (2.0 * pi * R2 - ZI * TB) * YI / (Q + QLK)

            if GV1 <= 0.0:
                GV1 = GF
            else:
                GV1 = (RS1 - 1.0 + np.sqrt((1.0 - RS1) ** 2 + 4.0 * RS1 * GF)) / (2.0 * RS1)

            if GV2 <= 0.0:
                GV2 = GF
            else:
                GV2 = (RS2 - 1.0 + np.sqrt((1.0 - RS2) ** 2 + 4.0 * RS2 * GF)) / (2.0 * RS2)

            GV = (GV1 + GV2) / 2
            CD_to_rb = (CD_to_rb1 + CD_to_rb2) / 2

        return GV, CD_to_rb

    @staticmethod
    def gv_sun_gamboa(B1, B2, DENL, DENG, N, Q, QG, QLK, R1, R2, ST, TB, VISL, VISG, YI1, YI2, YI, ZI):
        """
        Use Sun (2003) drag force coefficient and Gamoba (2008) flow pattern transition
        """

        GF = QG / (Q + QLK + QG)
        GV1, CD_to_rb = 1 - GF, 0.3
        GV2 = GV1
        GV = (GV1 + GV2) / 2
        ABV1, ABV2 = 1., 1.
        VSR1, VSR2 = 0.1, 0.1
        CD_to_rb1, CD_to_rb2 = 0.1, 0.1
        CD_to_rb = 0.1
        counter = 0
        OMEGA = 2.0 * pi * N / 60.0

        qmax = 2.91492 * N * bbl_to_m3 / 3600 / 24  # for GC-6100 only
        qld = (Q + QLK) / qmax
        lambdaC1 = (2.12e-6 * (2.4 * ST / (OMEGA ** 2 * (DENL - DENG))) ** (1./3) * (ST / DENL) ** (-0.6) *
                    ((OMEGA * 2 * R1) ** 4 / (VISL / DENL)) ** 0.4 * (DENL / DENG) ** (-0.2) - 0.0052) ** (1 / 0.1997)

        C1M_L = (Q + QLK) / ((2.0 * pi * R1 - ZI * TB) * YI1)
        C2M_L = (Q + QLK) / ((2.0 * pi * R2 - ZI * TB) * YI2)
        W1_L, W2_L = C1M_L / np.sin(B1), C2M_L / np.sin(B2)

        C1M_G = QG / ((2.0 * pi * R1 - ZI * TB) * YI1)
        C2M_G = QG / ((2.0 * pi * R2 - ZI * TB) * YI2)
        W1_G, W2_G = C1M_G / np.sin(B1), C2M_G / np.sin(B2)

        while (ABV1 > E1) or (ABV2 > E1):
            counter += 1
            if counter > 1000:
                return GV, CD_to_rb

            W1L = W1_L / (1 - GV1)
            W2L = W2_L / (1 - GV2)
            W1G = W1_G / GV1
            W2G = W2_G / GV2
            mium1 = (1 - GV1) * VISL + GV1 * VISG
            DENM1 = (1 - GV1) * DENL + GV1 * DENG
            mium2 = (1 - GV2) * VISL + GV2 * VISG
            DENM2 = (1 - GV2) * DENL + GV2 * DENG

            if GV1 <= 0.15:
                CD_to_rb1 = 12. * mium1 * (4.564e7 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)
            elif (GV1 > 0.15) and (GV1 <= 0.5):
                CD_to_rb1 = 12. * mium1 * (6.39e7 * (1 - GV1) ** 1.5 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)
            elif GV1 > 0.5:
                CD_to_rb1 = 12. * mium1 * (9.13e7 * GV1 ** 2 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)

            VSRN1 = np.sqrt(4.0 * (DENL - DENG) * R1 / (3.0 * CD_to_rb1 * DENL)) * OMEGA
            ABV1 = np.abs((VSRN1 - VSR1) / VSR1)
            VSR1 = VSRN1
            RS1 = VSR1 * (2.0 * pi * R1 - ZI * TB) * YI / (Q + QLK)

            if GV2 <= 0.15:
                CD_to_rb2 = 12. * mium2 * (4.564e7 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)
            elif (GV2 > 0.15) and (GV2 <= 0.5):
                CD_to_rb2 = 12. * mium2 * (6.39e7 * (1 - GV2) ** 1.5 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)
            elif GV2 > 0.5:
                CD_to_rb2 = 12. * mium2 * (9.13e7 * GV2 ** 2 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)

            VSRN2 = np.sqrt(4.0 * (DENL - DENG) * R2 / (3.0 * CD_to_rb2 * DENL)) * OMEGA
            ABV2 = np.abs((VSRN2 - VSR2) / VSR2)
            VSR2 = VSRN2
            RS2 = VSR2 * (2.0 * pi * R2 - ZI * TB) * YI / (Q + QLK)

            if GV1 <= 0.0:
                GV1 = GF
            else:
                GV1 = (RS1 - 1.0 + np.sqrt((1.0 - RS1) ** 2 + 4.0 * RS1 * GF)) / (2.0 * RS1)

            if GV2 <= 0.0:
                GV2 = GF
            else:
                GV2 = (RS2 - 1.0 + np.sqrt((1.0 - RS2) ** 2 + 4.0 * RS2 * GF)) / (2.0 * RS2)

            GV = (GV1 + GV2) / 2
            CD_to_rb = (CD_to_rb1 + CD_to_rb2) / 2

        return GV, CD_to_rb

    @staticmethod
    def gv_sun_zhu(B1, B2, DENL, DENG, N, Q, QG, QLK, R1, R2, TB, VISL, VISG, YI1, YI2, YI, ZI, lambdaC2, lambdaC3):
        """
        Modified Sun (2003) drag force coefficient
        """
        if lambdaC2 > lambdaC3:
            lambdaC2, lambdaC3 = lambdaC3, lambdaC2

        GF = QG / (Q + QLK + QG)
        GV1, CD_to_rb = 1-GF, 0.3
        GV2 = GV1
        GV = (GV1 + GV2) / 2
        ABV1, ABV2 = 1., 1.
        VSR1, VSR2 = 0.1, 0.1
        CD_to_rb = 0.1
        counter = 0
        OMEGA = 2.0 * pi * N / 60.0

        C1M_L = (Q + QLK) / ((2.0 * pi * R1 - ZI * TB) * YI1)
        C2M_L = (Q + QLK) / ((2.0 * pi * R2 - ZI * TB) * YI2)
        W1_L, W2_L = C1M_L / np.sin(B1), C2M_L / np.sin(B2)

        C1M_G = QG / ((2.0 * pi * R1 - ZI * TB) * YI1)
        C2M_G = QG / ((2.0 * pi * R2 - ZI * TB) * YI2)
        W1_G, W2_G = C1M_G / np.sin(B1), C2M_G / np.sin(B2)

        while (ABV1 > E1) or (ABV2 > E1):
            counter += 1
            if counter > 1000:
                return GV, CD_to_rb

            W1L = W1_L / (1 - GV1 + E1)
            W2L = W2_L / (1 - GV2 + E1)
            W1G = W1_G / GV1
            W2G = W2_G / GV2
            mium1 = (1 - GV1) * VISL + GV1 * VISG
            DENM1 = (1 - GV1) * DENL + GV1 * DENG
            mium2 = (1 - GV2) * VISL + GV2 * VISG
            DENM2 = (1 - GV2) * DENL + GV2 * DENG

            if GF <= lambdaC2:
                CD_to_rb1 = 12. * mium1 * (4.564e7 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)
            elif (GF > lambdaC2) and (GF <= lambdaC3):
                CD_to_rb1 = 12. * mium1 * (6.39e7 * (1 - GV1) ** 1.5 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)
            else:
                CD_to_rb1 = 12. * mium1 * (9.13e7 * GV1 ** 2 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)

            VSRN1 = np.sqrt(4.0 * (DENL - DENG) * R1 / (3.0 * CD_to_rb1 * DENL)) * OMEGA
            ABV1 = np.abs((VSRN1 - VSR1) / VSR1)
            VSR1 = VSRN1
            RS1 = VSR1 * (2.0 * pi * R1 - ZI * TB) * YI / (Q + QLK)

            if GF <= lambdaC2:
                CD_to_rb2 = 12. * mium2 * (4.564e7 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)
            elif (GV2 > lambdaC2) and (GF <= lambdaC3):
                CD_to_rb2 = 12. * mium2 * (6.39e7 * (1 - GV2) ** 1.5 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)
            else:
                CD_to_rb2 = 12. * mium2 * (9.13e7 * GV2 ** 2 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)

            VSRN2 = np.sqrt(4.0 * (DENL - DENG) * R2 / (3.0 * CD_to_rb2 * DENL)) * OMEGA
            ABV2 = np.abs((VSRN2 - VSR2) / VSR2)
            VSR2 = VSRN2
            RS2 = VSR2 * (2.0 * pi * R2 - ZI * TB) * YI / (Q + QLK)

            if GV1 <= 0.0:
                GV1 = GF
            else:
                GV1 = (RS1 - 1.0 + np.sqrt((1.0 - RS1) ** 2 + 4.0 * RS1 * GF)) / (2.0 * RS1)

            if GV2 <= 0.0:
                GV2 = GF
            else:
                GV2 = (RS2 - 1.0 + np.sqrt((1.0 - RS2) ** 2 + 4.0 * RS2 * GF)) / (2.0 * RS2)

            GV = (GV1 + GV2) / 2
            CD_to_rb = (CD_to_rb1 + CD_to_rb2) / 2

        return GV, CD_to_rb

    @staticmethod
    def gv_barrios(DENL, DENG, GF, N, Q, QLK, R1, RI, ST, TB, VISL, YI, ZI):
        """
        Use Barrios (2007) drag force coefficient, which is the original form Dr Zhang's development
        :param DENL: liquid density in m3/s
        :param DENG: gas density in m3/s
        :param GF: gas volumetric fraction at pump inlet
        :param N: rotatonal speed in rpm
        :param Q: liquid flow rate in m3/s
        :param QLK: leakage flow rate in m3/s
        :param R1: radius of impeller at inlet in m
        :param RI: radius of impeller in m
        :param ST: interfacial tension in N/m
        :param TB: impeller blade thickness in m
        :param VISL: liquid viscosity in Pas
        :param YI: impeller outlet height
        :param ZI: impeller blade number
        :return: GV (in-situ gas void fraction), CD_to_rb(CD over rb)
        """
        GV, DB= GF, 0.3
        CD = 0.
        ABV = 1.
        VSR = 0.1
        counter = 0
        OMEGA = 2.0 * pi * N / 60.0

        while ABV > E1:
            counter += 1
            if counter > 10000:
                return GV, CD/DB

            # Bubble diameter - Barrios 2007 bubble diameter at surging
            # DB = 0.0348 * N ** 0.8809 * GF ** 0.25 * (self.ST / DENL) ** 0.6 / (N ** 3 * R1 ** 2) ** 0.4

            # A modification on bubble size is to use GV instead of GF
            DB = 3.0 * 0.0348 * N ** 0.8809 * GV ** 0.25 * (ST / DENL) ** 0.6 / (N ** 3 * R1 ** 2) ** 0.4

            # Gamboa bubble size model
            # DB = 14.2667*(ST/DENL)**0.6*((OMEGA*2*RI)**4/(VISL/DENL))**(-0.4)*(DENL/DENG)**0.2*(1+191.7*GV**0.1997)
            REB = DENL * VSR * DB / VISL

            if REB <= 0.1:
                REB = 0.1

            # Drag coefficient of gas bubble in radial direction
            YY = 0.00983 + 389.9 * REB / N ** 2
            CD = (24.0 + 5.48 * (REB * YY) ** 0.427 + 0.36 * REB * YY) / (REB * YY)

            VSRN = np.sqrt(4.0 * DB * (DENL - DENG) * RI / (3.0 * CD * DENL)) * OMEGA
            ABV = np.abs((VSRN - VSR) / VSR)
            VSR = VSRN
            RS = VSR * (2.0 * pi * RI - ZI * TB) * YI / (Q + QLK)

            if GV < 0.0:
                GV = 0.0
            else:
                GV = (RS - 1.0 + np.sqrt((1.0 - RS) ** 2 + 4.0 * RS * GF)) / (2.0 * RS)

        return GV, CD/DB

    def gl_calculate_old(self, Q, QG, flg='Z'):
        """
        Calcualte air-water flow performance of ESP
        :param QG: gas flow rate (bpd)
        :param Q: liquid flow rate (bpd)
        :param flg: 'Z': Zhu model; 'S': Sun model; 'B': Barrios model
        :return: PP, PE, PF, PT, PD, PRE, PLK, QLK, GV
        """
        # run single-phase calculation to initialize
        HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_old(Q)

        # convert field units
        PP = HP * psi_to_pa
        Q = Q * bbl_to_m3 / 24.0 / 3600.0
        QG = QG * bbl_to_m3 / 24.0 / 3600.0
        GF = Q / (Q + QG)

        icon = 0
        ABP = 1.
        PE, PEE = 0., 0
        PFI, PFD = 0., 0
        PTI, PTD = 0, 0
        GV = 0.
        PLK = 0.

        AIW = self.AB + self.ASF + self.ASB
        ADW = self.AV + self.ADF + self.ADB
        AI = self.VOI / self.LI
        AD = self.VOD / self.LD
        DI = 4.0 * self.VOI / AIW
        DD = 4.0 * self.VOD / ADW
        EDI = self.EA / DI
        EDD = self.EA / DD
        DEND = self.DENL * (1.0 - GF) + self.DENG * GF

        while ABP > E1:
            VI = (Q + QLK) / self.ZI / AI
            VD = Q / self.ZD / AD
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

            HE = (U2 ** 2 - U2 * C2M / np.tan(self.B2)) / G

            if GF <= E1:
                PP = HP * psi_to_pa
                break
            elif GF >= (1.-E1):
                PP = HP * psi_to_pa * self.DENG / DENW
                break

            if flg == 'Z':
                GV = self.bubble_flow_zhu(HP * psi_to_pa, GF, Q, QLK)
            elif flg == 'S':
                GV = self.bubble_flow_sun(GF, Q, QLK)
            elif flg == 'B':
                GV = self.bubble_flow_barrios(GF, Q, QLK)
            else:
                GV = self.bubble_flow_zhu(HP * psi_to_pa, GF, Q, QLK)

            DENI = self.DENL * (1.0 - GV) + self.DENG * GV
            PE = HE * DENI * G
            VLR = (Q + QLK) / ((2.0 * pi * self.RI - self.ZI * self.TB) * self.YI * (1.0 - GV))
            VGR = QG / ((2.0 * pi * self.RI - self.ZI * self.TB) * self.YI * GV)

            if Q < self.QBEM:
                VSH = U2 * (self.QBEM - Q) / self.QBEM
                C2F = C2B * Q / self.QBEM
                DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
                REC1 = DENI * VSH * DC / self.VISL
                # SGM = (self.VISW / self.VISL) ** 0.5 / (1.0 + 0.02 * REC1 ** 0.2)
                SGM = 0.0
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F)
                PEE = PE + DENI * (C2E ** 2 - C2 ** 2) / 2.0
            else:
                VSH = U2 * (Q - self.QBEM) / self.QBEM
                C2F = C2B * Q / self.QBEM
                C2E = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                PEE = PE + DENI * (C2E ** 2 - C2 ** 2) / 2.0

            REI = DENI * (W1 + W2) * DI / self.VISL / 2.0
            RED = DEND * VD * DD / self.VISL
            FFI = self.get_fff(REI, EDI)
            FFD = self.get_fff(RED, EDD)
            PFI = 1.5 * 4.0 * FFI * DENI * VI**2 * self.LI / (2.0*DI)
            PFD = 1.5 * 4.0 * FFD * DEND * VD ** 2 * self.LD / (2.0 * DD)

            if (self.DENL*4*self.R2**2*self.OMEGA**2/self.VISL) > 1e7:
                FTI = 6.5
                FTD = 2.5
            else:
                FTI = 12.
                FTD = 3.

            PTI = FTI * DENI * VI ** 2 / 2.0
            PTD = FTD * DEND * VD ** 2 / 2.0
            PPN = PEE - PFI - PFD - PTI - PTD

            UL = self.RLK * self.OMEGA
            PIO = PEE - PFI - PTI
            PLK = PIO - DENI * (U2 ** 2 - UL ** 2) / 8.0

            if PLK >= 0:
                VL = np.abs(QLK) / (2.0 * pi * self.RLK * self.SL)
                REL = np.abs(DEND * VL * self.SL / self.VISL)
                EDL = 0.0
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * PLK / (1.5 + 4.0 * FFL * self.LG / self.SL) / DEND)
                QLKN = 2.0 * pi * self.RLK * self.SL * VL
            else:
                VL = np.abs(QLK / (2.0 * pi * self.RLK * self.SL))
                REL = self.DENL * VL * self.SL / self.VISL
                EDL = 0.0
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * np.abs(PLK) / (1.5 + 4.0 * FFL * self.LG / self.SL) / DEND)
                QLKN = -2.0 * pi * self.RLK * self.SL * VL

            QLK = QLKN
            ABP = np.abs((PPN-PP)/PPN)
            PP = PPN

            if icon > 500:
                break
            else:
                icon += 1

        # return pressure in psi, flow rate in bpd
        PP = PP / psi_to_pa
        PE = PE / psi_to_pa
        PEE = PEE / psi_to_pa
        PF = (PFI + PFD) / psi_to_pa
        PT = (PTI + PTD) / psi_to_pa
        PD = (PFD + PTD) / psi_to_pa
        PRE = np.abs(PE - PEE)
        PLK = PLK / psi_to_pa
        QLK = QLK * 24.0 * 3600.0 / bbl_to_m3
        return PP, PE, PF, PT, PD, PRE, PLK, QLK, GV

    def gl_calculate_new(self, QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW, WC, flg):
        """
        Calcualte gas-liquid performance of ESP
        :param QL:  liquid flow rate in m3/s
        :param QG:  gas flow rate m3/s
        :param QEM: best match flow rate in m3/s
        :param flg: 'Z': Zhu model; 'S': Sun model; 'B': Barrios model
        :param DENG: gas density in kg/m3
        :param DENL: liquid density in kg/m3
        :param DENW: water density in kg/m3
        :param N:   rotational speed in rpm
        :param NS:  specific speed based on field units
        :param SGM: tuning factor
        :param SN:  stage number
        :param ST:  surface tension Nm
        :param VISG: gas viscosity in Pas
        :param VISL: liquid viscosity in Pas
        :param VISM: water viscosity in Pas
        :param WC:  water cut in %
        :return: PP, PE, PF, PT, PD, PRE, PLK, QLK, GV in field units
        """
        # # a tricky adjusting of rotational speed in Sun's model
        # if (flg == 'S') and N < 3500:
        #     N += 600

        # run single-phase calculation to initialize
        HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_new(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW,
                                                                   WC)

        # convert filed units
        PP = HP * psi_to_pa
        GF = QG / (QL + QG)

        icon = 0
        ABP = 1.
        PE, PEE = 0., 0
        PFI, PFD = 0., 0
        PTI, PTD = 0, 0
        GV = 0.
        PLK = 0.
        QLK = 0.02 * QL
        OMEGA = 2.0 * pi * N / 60.0

        AIW = self.AB + self.ASF + self.ASB
        ADW = self.AV + self.ADF + self.ADB
        AI = self.VOI / self.LI
        AD = self.VOD / self.LD
        DI = 4.0 * self.VOI / AIW
        DD = 4.0 * self.VOD / ADW
        EDI = self.EA / DI
        EDD = self.EA / DD
        DEND = DENL * (1.0 - GF) + DENG * GF

        # new QBEM due to liquid viscosity

        QBEM = QBEM * (N / 3600) * (VISL / VISW) ** (0.01 * (3448 / NS) ** 4)

        # check if emulsion occurs
        if WC > 0.:
            VISL = self.emulsion(self.VOI, self.R2, VISL, VISW, DENL, DENW, WC, ST, N, QL, SN)

        while ABP > E1:
            VI = (QL + QLK) / self.ZI / AI
            VD = QL / self.ZD / AD
            C1M = (QL + QLK) / ((2.0 * pi * self.R1 - self.ZI * self.TB) * self.YI1)
            C2M = (QL + QLK) / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            U1 = self.R1 * OMEGA
            U2 = self.R2 * OMEGA
            W1 = C1M / np.sin(self.B1)
            W2 = C2M / np.sin(self.B2)
            C1 = np.sqrt(C1M ** 2 + (U1 - C1M / np.tan(self.B1)) ** 2)
            C2 = np.sqrt(C2M ** 2 + (U2 - C2M / np.tan(self.B2)) ** 2)
            CMB = QBEM / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
            C2B = np.sqrt(CMB ** 2 + (U2 - CMB / np.tan(self.B2)) ** 2)

            HE = (U2 ** 2 - U2 * C2M / np.tan(self.B2)) / G

            if GF <= E1:
                PP = HP * psi_to_pa
                break
            elif GF >= (1-E1):
                PP = HP * psi_to_pa * DENG / DENW
                break

            if flg == 'Z':
                qld = (QL + QLK) / (2.91492 * N * bbl_to_m3 / 3600 / 24)        # for GC-6100 only
                qgd = (DENG / DENL) ** 0.2 * (OMEGA * (2 * self.R1) ** 2 / (VISL / DENL)) ** 0.4 * \
                      (0.102 * np.exp(qld)) ** 4.4682
                lambdaC2 = qgd / (qgd + qld)
                lambdaC3 = self.get_lambda_c3(QL + QLK)
                GV, _ = self.gv_sun_zhu(self.B1, self.B2, DENL, DENG, N, QL, QG, QLK, self.R1, self.R2, self.TB,
                                        VISL, VISG, self.YI1, self.YI2, self.YI, self.ZI, lambdaC2, lambdaC3)
            elif flg == 'S':
                GV, _ = self.gv_sun(self.B1, self.B2, DENL, DENG, N, QL, QG, QLK, self.R1, self.R2, self.TB,
                                    VISL, VISG, self.YI1, self.YI2, self.YI, self.ZI)
            elif flg == 'G':
                GV, _ = self.gv_sun_gamboa(self.B1, self.B2, DENL, DENG, N, QL, QG, QLK, self.R1, self.R2, ST, self.TB,
                                           VISL, VISG, self.YI1, self.YI2, self.YI, self.ZI)
            else:
                GV, _ = self.gv_sun(self.B1, self.B2, DENL, DENG, N, QL, QG, QLK, self.R1, self.R2, self.TB,
                                    VISL, VISG, self.YI1, self.YI2, self.YI, self.ZI)

            DENI = DENL * (1.0 - GV) + DENG * GV
            PE = HE * DENI * G
            VLR = (QL + QLK) / ((2.0 * pi * self.RI - self.ZI * self.TB) * self.YI * (1.0 - GV))
            VGR = QG / ((2.0 * pi * self.RI - self.ZI * self.TB) * self.YI * GV)

            if (QL + QLK) <= QBEM:
                VSH = U2 * (QBEM - (QL + QLK)) / QBEM
                C2F = C2B * (QL + QLK) / QBEM
                DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
                REC1 = DENI * VSH * DC / VISL
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F
                PEE = PE + DENI * (C2E ** 2 - C2 ** 2) / 2.0
            else:
                VSH = U2 * (QL + QLK - QBEM) / QBEM
                C2F = C2B * (QL + QLK) / QBEM
                C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
                C2E = C2F + SGM * (C2P - C2F) * (QL + QLK - QBEM) / QBEM
                PEE = PE + DENI * (C2E ** 2 - C2 ** 2) / 2.0

            REI = DENI * (W1 + W2) * DI / VISL / 2.0 / 2.0
            RED = DEND * VD * DD / VISL / 2.0
            FFI = self.get_fff(REI, EDI)
            FFD = self.get_fff(RED, EDD)
            PFI = 2.5 * 4.0 * FFI * DENI * VI ** 2 * self.LI / (2.0 * DI)
            PFD = 2.5 * 4.0 * FFD * DEND * VD ** 2 * self.LD / (2.0 * DD)

            FTI = 3.0
            FTD = 3.0
            PTI = FTI * DENI * VI ** 2 / 2.0
            PTD = FTD * DEND * VD ** 2 / 2.0
            PPN = PEE - PFI - PFD - PTI - PTD

            UL = self.RLK * OMEGA
            PIO = PEE - PFI - PTI
            PLK = PIO - DENI * (U2 ** 2 - UL ** 2) / 8.0

            if PLK >= 0.:
                VL = np.abs(QLK) / (2.0 * pi * self.RLK * self.SL)
                REL = np.abs(DEND * VL * self.SL / VISL)
                EDL = 0.0
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * PLK / (1.5 + 4.0 * FFL * self.LG / self.SL) / DEND)
                QLKN = 2.0 * pi * self.RLK * self.SL * VL
            else:
                VL = np.abs(QLK / (2.0 * pi * self.RLK * self.SL))
                REL = DENL * VL * self.SL / VISL
                EDL = 0.0
                FFL = self.get_fff(REL, EDL)
                VL = np.sqrt(2.0 * np.abs(PLK) / (1.5 + 4.0 * FFL * self.LG / self.SL) / DEND)
                QLKN = -2.0 * pi * self.RLK * self.SL * VL

            QLK = QLKN
            ABP = np.abs((PPN - PP) / PPN)
            PP = PPN

            if icon > 500:
                break
            else:
                icon += 1

        # return pressure in psi, flow rate in bpd
        PP = PP / psi_to_pa
        PE = PE / psi_to_pa
        PEE = PEE / psi_to_pa
        PF = (PFI + PFD) / psi_to_pa
        PT = (PTI + PTD) / psi_to_pa
        PD = (PFD + PTD) / psi_to_pa
        PRE = np.abs(PE - PEE)
        PLK = PLK / psi_to_pa
        QLK = QLK * 24.0 * 3600.0 / bbl_to_m3
        return PP, PE, PF, PT, PD, PRE, PLK, QLK, GV

    def flow_pattern(self):
        QL = np.arange(500)[1:] * 50.0
        QSG1 = []
        QSG2 = []
        QSG3 = []
        for ql in QL:
            # ql in bpd, HP in psi
            HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_old(ql)
            if HP > 0:
                Q = (ql + QLK) * bbl_to_m3 / 24 / 3600      # convert bpd to m3/s
                LambdaC1 = self.get_lambda_c1(Q, HP * psi_to_pa)
                LambdaC2 = self.get_lambda_c2(Q, HP * psi_to_pa)
                LambdaC3 = self.get_lambda_c3(Q)

                if LambdaC1 > LambdaC2: LambdaC2 = LambdaC1

                QSG1.append(LambdaC1 / (1.0 - LambdaC1) * ql)
                QSG1.append(LambdaC2 / (1.0 - LambdaC2) * ql)
                QSG1.append(LambdaC3 / (1.0 - LambdaC3) * ql)
            else:
                break
        return QL[:len(QSG1)], QSG1, QSG2, QSG3


#################################
class SinglePhaseCompare(object):
    def __init__(self, pump, conn):
        # QBEM = {'TE2700': 5600, 'DN1750': 3300, 'GC6100': 8800, 'P100': 12000}      # sgl_new
        QBEM = {'TE2700': 4500, 'DN1750': 3000, 'GC6100': 7800, 'P100': 11000}    # sgl_2018
        self.pump = pump
        self.ESP = ESP[pump]
        self.conn = conn
        self.df_catalog = pd.read_sql_query("SELECT * FROM Catalog_All;", self.conn)
        self.df_pump = self.df_catalog[self.df_catalog['Pump'] == pump]
        # print(self.df_pump)
        self.QBEM = QBEM[pump]

    def single_phase_water(self):
        sgl = SinglePhaseModel(self.ESP, self.QBEM)
        QL, hpsgl, _, _, _, _, _, _, _ = sgl.performance_curve()
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(111)
        ax.plot(QL, hpsgl, 'b-', label='model')
        ax.plot(self.df_pump['Flow_bpd'], self.df_pump['DP_psi'], 'ro', label='catalog')
        ax.set_xlabel(r'$Q_L$ (bpd)')
        ax.set_ylabel(r'$P$ (psi)')
        ax.set_title(r'{} ESP, $Q_B$={} bpd'.format(self.pump, self.QBEM))
        ax.legend(frameon=False)
        plt.show()

    def single_phase_viscous(self, QL, N, DENL, VISL, WC):
        """
        :param QL: the liquid flow rate in bpd
        :param DENL: liquid density in kg/m3
        :param VISL: liquid viscosity in cP
        :return: HP, a single-point calculation based on the input flow conditions
        """
        sgl = SinglePhaseModel(self.ESP, self.QBEM)
        NS = self.ESP['NS']
        SGM = self.ESP['SGM']
        SN = self.ESP['SN']
        ST = self.ESP['ST']
        QBEM = self.QBEM * bbl_to_m3 / 24.0 / 3600.0
        QL = QL * bbl_to_m3 / 24.0 / 3600.0
        VISL = VISL / 1000.
        # HP, _, _, _, _, _, _, _ = sgl.sgl_calculate_new(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)
        HP, _, _, _, _, _, _, _ = sgl.sgl_calculate_2018(QL, QBEM, DENL, DENW, N, NS, SGM, SN, 0.03, VISL, VISW, WC)
        return HP

    def viscous_fluid_performance(self, N, DENL, VISL, WC):
        """
        calculate single-phase viscous fluid flow H-Q performance with plotting
        :param vis: viscosity in Pas
        :param den: density in kg/m3
        :param WC: water cut
        :return: df_model
        """
        sgl = SinglePhaseModel(self.ESP, self.QBEM)
        NS = self.ESP['NS']
        SGM = self.ESP['SGM']
        SN = self.ESP['SN']
        ST = self.ESP['ST']
        QBEM = self.QBEM * bbl_to_m3 / 24.0 / 3600.0

        QL = np.arange(5000)[1:] * 50.0
        hpsgl = []
        VISL = VISL /1000.

        # q in bpd
        for ql in QL:
            # q in bpd
            ql = ql * bbl_to_m3 / 24.0 / 3600.0
            # HP, _, _, _, _, _, _, _ = sgl.sgl_calculate_new(ql, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)
            HP, _, _, _, _, _, _, _ = sgl.sgl_calculate_2018(ql, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)
            if HP > 0:
                hpsgl.append(HP)
            else:
                break
        df_model = pd.DataFrame({'QL': QL[:len(hpsgl)], 'DP': hpsgl})
        return df_model


#################################
class TwoPhaseCompare(object):
    def __init__(self, pump, conn):
        QBEM = {'TE2700': 8000, 'DN1750': 3300, 'GC6100': 8800, 'P100': 12000}  # gl_new
        QBEP = {'TE2700': 2700, 'DN1750': 1750, 'GC6100': 6100, 'P100': 9000}  # bep flow rate
        self.pump = pump
        self.ESP = ESP[pump]
        self.conn = conn
        self.QBEM = QBEM[pump]
        self.QBEP = QBEP[pump]

    def surging_cd_to_db(self):
        """
        :return: two dataframes for GV and CD_over_dB, the column names: zhu, barrios, sun
        """
        sgl = SinglePhaseModel(self.ESP, self.QBEM)
        gl = GasLiquidModel(self.ESP, self.QBEM)
        sgl_cal = np.vectorize(sgl.sgl_calculate_new)
        zhu_cal = np.vectorize(gl.gv_zhu)
        sun_cal = np.vectorize(gl.gv_sun)
        barrios_cal = np.vectorize(gl.gv_barrios)

        # use the best efficiency point to compare
        GF = np.arange(0.01, 1, 0.02)
        QL = gl.QL * bbl_to_m3 / 24.0 / 3600.0 * np.ones(GF.shape)
        QG = GF / (1 - GF) * QL
        QBEM = self.QBEM * bbl_to_m3 / 24.0 / 3600.0 * np.ones(GF.shape)
        QLK = 0 * np.ones(GF.shape)
        DENL = gl.DENL * np.ones(GF.shape)
        DENG = gl.DENG * np.ones(GF.shape)
        N = gl.N * np.ones(GF.shape)
        NS = gl.NS * np.ones(GF.shape)
        SGM = gl.SGM * np.ones(GF.shape)
        SN = gl.SN * np.ones(GF.shape)
        ST = gl.ST * np.ones(GF.shape)
        VISL = gl.VISL * np.ones(GF.shape)
        VISG = gl.VISG * np.ones(GF.shape)
        WC = gl.WC * np.ones(GF.shape)

        R1 = gl.R1 * np.ones(GF.shape)
        R2 = gl.R2 * np.ones(GF.shape)
        YI1 = gl.YI1 * np.ones(GF.shape)
        YI2 = gl.YI2 * np.ones(GF.shape)
        RI = gl.RI * np.ones(GF.shape)
        VOI = gl.VOI * np.ones(GF.shape)
        YI = gl.YI * np.ones(GF.shape)
        ZI = gl.ZI * np.ones(GF.shape)
        TB = gl.TB * np.ones(GF.shape)
        B1 = gl.B1 * np.ones(GF.shape)
        B2 = gl.B2 * np.ones(GF.shape)

        HP, _, _, _, _, _, _, _ = sgl_cal(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)
        gv_zhu, cd_over_db_zhu = zhu_cal(DENL, DENG, GF, HP * psi_to_pa, N, QL, QLK, RI, ST, TB, VISL, VOI, YI, ZI)
        gv_sun, cd_over_db_sun = sun_cal(B1, B2, DENL, DENG, N, QL, QG, QLK, R1, R2, RI, TB, VISL, VISG,
                                         YI1, YI2, YI, ZI)
        gv_barrios, cd_over_db_barrios = barrios_cal(DENL, DENG, GF, N, QL, QLK, R1, RI, ST, TB, VISL, YI, ZI)

        df_gv = pd.DataFrame({'zhu': gv_zhu, 'sun': gv_sun, 'barrios': gv_barrios, 'gf': GF})
        df_cd_over_db = pd.DataFrame({'zhu': cd_over_db_zhu, 'sun': cd_over_db_sun,
                                      'barrios': cd_over_db_barrios, 'gf': GF})
        return df_gv, df_cd_over_db

    def surging_performance(self, QL, maxGF, N, p, t):
        """
        :param QL: liquid flow rate in bpd
        :param maxGF: maximum GF for calculation
        :param N: array for rotational speed rpm
        :param p: array for gas pressure psi
        :param t: array for temperature F
        :return: dataframe of predicted pump heads under surging flow, the column names: zhu, barrios, sun
        """
        sgl = SinglePhaseModel(self.ESP, self.QBEM)
        gl = GasLiquidModel(self.ESP, self.QBEM)
        sgl_cal = np.vectorize(sgl.sgl_calculate_new)
        gl_cal = np.vectorize(gl.gl_calculate_new)

        GF = np.arange(0.0, maxGF + 0.02, 0.01)
        QL = QL * bbl_to_m3 / 24.0 / 3600.0 * np.ones(GF.shape)
        QBEM = self.QBEM * bbl_to_m3 / 24.0 / 3600.0 * np.ones(GF.shape)
        QG = GF / (1 - GF) * QL
        DENL = gl.DENL * np.ones(GF.shape)
        DENG = gasdensity(p, t, 0)
        DENG = DENG * np.ones(GF.shape)
        N = N * np.ones(GF.shape)
        NS = gl.NS * np.ones(GF.shape)
        SGM = gl.SGM * np.ones(GF.shape)
        SN = gl.SN * np.ones(GF.shape)
        ST = gl.ST * np.ones(GF.shape)
        VISL = gl.VISL * np.ones(GF.shape)
        VISG = gl.VISG * np.ones(GF.shape)
        WC = gl.WC * np.ones(GF.shape)

        # flg for different drag coefficient model
        flgz = np.empty(GF.shape, dtype='str')
        flgs = np.empty(GF.shape, dtype='str')
        flgb = np.empty(GF.shape, dtype='str')
        flgz[:] = 'Z'
        flgb[:] = 'B'
        flgs[:] = 'S'

        HP, _, _, _, _, _, _, _ = sgl_cal(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)

        PPS, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgs)

        df = pd.DataFrame({'gf': GF, 'sun': PPS / HP})

        return df

    def error_analysis(self, QL, GVF, N, p, t):
        """
        :param QL: array for liquid flow rate bpd
        :param QG: array for gas friction factor %
        :param N: array for rotational speed rpm
        :param p: array for gas pressure psi
        :param t: array for temperature F
        :return: a data frame for pressure increment
        """
        gl = GasLiquidModel(self.ESP, self.QBEM)
        gl_cal = np.vectorize(gl.gl_calculate_new)

        QL = QL * bbl_to_m3 / 24.0 / 3600.0
        QBEM = self.QBEM * bbl_to_m3 / 24.0 / 3600.0 * np.ones(QL.shape)
        QG = GVF / (100 - GVF) * QL
        DENL = gl.DENL * np.ones(QL.shape)
        h = np.zeros(QL.shape)
        DENG = gasdensity(p, t, h)
        NS = gl.NS * np.ones(QL.shape)
        SGM = gl.SGM * np.ones(QL.shape)
        SN = gl.SN * np.ones(QL.shape)
        ST = gl.ST * np.ones(QL.shape)
        VISL = gl.VISL * np.ones(QL.shape)
        VISG = gl.VISG * np.ones(QL.shape)
        WC = gl.WC * np.ones(QL.shape)

        # flg for different drag coefficient model
        flgs = np.empty(QL.shape, dtype='str')
        flgs[:] = 'S'

        PPS, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgs)
        df = pd.DataFrame({'sun': PPS})
        return df

    def mapping_performance(self, QG, maxQL, N, p, t):
        """
        :param QG: constant gas flow rate bpd
        :param maxQL: maximum liquid flow rate bpd
        :param N: rotational speed rpm
        :param p: array for gas pressure psi
        :param t: array for temperature F
        :return: dataframe of predicted pump heads under mapping flow, the column names: zhu, barrios, sun
        """
        gl = GasLiquidModel(self.ESP, self.QBEM)
        gl_cal = np.vectorize(gl.gl_calculate_new)

        QL = np.arange(0.01, 1, 0.02) * maxQL * bbl_to_m3 / 24.0 / 3600.0
        QG = QG * bbl_to_m3 / 24.0 / 3600.0 * np.ones(QL.shape)

        QBEM = self.QBEM * bbl_to_m3 / 24.0 / 3600.0 * np.ones(QL.shape)
        DENL = DENW * np.ones(QL.shape)
        DENG = gasdensity(p, t, 0)
        DENG = DENG * np.ones(QL.shape)
        NS = gl.NS * np.ones(QL.shape)
        N = N * np.ones(QL.shape)
        SGM = gl.SGM * np.ones(QL.shape)
        SN = gl.SN * np.ones(QL.shape)
        ST = gl.ST * np.ones(QL.shape)
        VISL = gl.VISL * np.ones(QL.shape)
        VISG = gl.VISG * np.ones(QL.shape)
        WC = gl.WC * np.ones(QL.shape)

        # flg for different drag coefficient model
        flgz = np.empty(QL.shape, dtype='str')
        flgs = np.empty(QL.shape, dtype='str')
        flgb = np.empty(QL.shape, dtype='str')
        flgz[:] = 'Z'
        flgb[:] = 'B'
        flgs[:] = 'S'

        # PPZ, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
        #                                      WC, flgz)
        PPS, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgs)
        # PPB, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM,
        #                                      SN, ST, VISG, VISL, VISW, WC, flgb)

        # df = pd.DataFrame({'ql': QL/bbl_to_m3 * 3600 * 24, 'zhu': PPZ, 'sun': PPS, 'barrios': PPB})
        df = pd.DataFrame({'ql': QL / bbl_to_m3 * 3600 * 24, 'sun': PPS})

        return df


# plot TE2700 ESP flow case
def vis_te2700_plot():
    conn, c = connect_db('ESP.db')
    te2700_case = SinglePhaseCompare('TE2700', conn)

    df_jc = pd.read_sql_query("SELECT * FROM TE2700_JiechengZhang;", conn)
    df_3500 = df_jc[(df_jc.RPM == 3500)]
    df_2400 = df_jc[(df_jc.RPM == 2400)]
    df_all = [df_3500, df_2400]

    for index, df in enumerate(df_all):
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(111)
        id = 0

        for visl in df.viscosity_cP.unique():
            denl = df[df.viscosity_cP == visl]['Density_kg/m3'].mean()  # average density
            df_model = te2700_case.viscous_fluid_performance(3500 - 1100 * index, denl, visl, 0)
            df_turzo = turzo_2000('TE2700', visl, denl, 3500 - 1100 * index)
            ax.scatter(df[df.viscosity_cP == visl]['flow_bpd'], df[df.viscosity_cP == visl]['DP_psi'],
                        marker=symbols[id], label=r'{} cp_exp'.format(visl), facecolors='none',
                        edgecolor='C{}'.format(id + 1), linewidths=0.75)
            ax.plot(df_model.QL, df_model.DP, c='C{}'.format(id + 1), label=r'{} cp_model'.format(visl), linewidth=0.75)
            ax.plot(df_turzo.Qvis, df_turzo.DPvis, c='C{}'.format(id + 1), label=r'{} cp_Turzo (2000)'.format(visl),
                     linestyle='--', linewidth=0.75)
            id += 1

        handles, labels = ax.get_legend_handles_labels()
        handles = [handles[10], handles[11], handles[12], handles[13], handles[14],
                   handles[0], handles[2], handles[4], handles[6], handles[8],
                   handles[1], handles[3], handles[5], handles[7], handles[9]]

        labels = [labels[10], labels[11], labels[12], labels[13], labels[14],
                  labels[0], labels[2], labels[4], labels[6], labels[8],
                  labels[1], labels[3], labels[5], labels[7], labels[9]]

        ax.set_xlabel(r'$Q_L$ (bpd)')
        ax.set_ylabel(r'$P$ (psi)')
        ax.legend(handles, labels, frameon=False, fontsize=5)
        # fig.show()

    # plot best match line curve
    fig3 = plt.figure(dpi=300)
    ax3 = fig3.add_subplot(111)
    model_predict = []

    for i in range(df_jc.shape[0]):
        ql = df_jc.iloc[i].flow_bpd
        N = df_jc.iloc[i].RPM
        DENL = df_jc.iloc[i]['Density_kg/m3']
        VISL = df_jc.iloc[i].viscosity_cP
        model_predict.append(te2700_case.single_phase_viscous(ql, N, DENL, VISL, 0))

    df_jc['model'] = model_predict
    x_max = np.array(range(int(max(df_jc.model) + 1)))

    df_2400 = df_jc[df_jc.RPM == 2400]
    df_3500 = df_jc[df_jc.RPM == 3500]
    ax3.scatter(df_2400.DP_psi, df_2400.model, facecolor='none', edgecolor='b',
                marker=symbols[0], label=r'$N=2400$ rpm', linewidth=0.75)
    ax3.scatter(df_3500.DP_psi, df_3500.model, facecolor='none', edgecolor='r',
                marker=symbols[1], label=r'$N=3500$ rpm', linewidth=0.75)
    ax3.plot(x_max, x_max, 'k--', label='perfect match')
    ax3.plot(x_max, 1.25 * x_max, ls=':', c='gray')
    ax3.plot(x_max, 0.85 * x_max, ls=':', c='gray')
    ax3.set_xlabel(r'$P_{exp}$ (psi)')
    ax3.set_ylabel(r'$P_{model}$ (psi)')
    ax3.legend(frameon=False, fontsize=8)
    fig3.show()

    # epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_jc.model, df_jc.DP_psi)
    # print(epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6)

    disconnect_db(conn)


def gl_te2700_plot():
    conn, c = connect_db('ESP.db')
    te2700_case = TwoPhaseCompare('TE2700', conn)
    dfgv, dfcd = te2700_case.surging_cd_to_db()

    df_te2700_surging = pd.read_sql_query("SELECT * FROM TE2700_JJZ_Surging;", conn)
    df_te2700_mapping = pd.read_sql_query("SELECT * FROM TE2700_JJZ_Mapping;", conn)

    fig1 = plt.figure(dpi=300)
    ax1 = fig1.add_subplot(111)
    ax1.plot(dfgv.gf * 100, dfgv.sun * 100, 'r--', label='Sun (2003)')
    ax1.plot(dfgv.gf * 100, dfgv.barrios * 100, 'k-.', label='Barrios (2007)')
    ax1.plot(dfgv.gf * 100, dfgv.zhu * 100, 'b-', label='Zhu & Zhang (2015)')

    ax1.set_xlabel(r'$\lambda_{G}$ (%)')
    ax1.set_ylabel(r'$\alpha_{G}$ (%)')
    ax1.legend(frameon=False)

    fig2 = plt.figure(dpi=300)
    ax2 = fig2.add_subplot(111)
    ax2.plot(dfgv.gf * 100, dfcd.sun, 'r--', label='Sun (2003)')
    ax2.plot(dfgv.gf * 100, dfcd.barrios, 'k-.', label='Barrios (2007)')
    ax2.plot(dfgv.gf * 100, dfcd.zhu, 'b-', label='Zhu & Zhang (2015)')

    ax2.set_xlabel(r'$\lambda_{G}$ (%)')
    ax2.set_ylabel(r'$\frac{C_{d}}{r_{b}}$ $(m^{-1})$')
    ax2.legend(frameon=False)
    ax2.set_yscale('log')

    # compare surging test for example, N = 3500, Psep = 150 psi, QBEP
    df_3500_QBEP_150 = df_te2700_surging[(df_te2700_surging.Psep == 150) & (df_te2700_surging.RPM == 3500) &
                                         (df_te2700_surging.QBEP == 0.75)]
    df_model = te2700_case.surging_performance(df_3500_QBEP_150.QL.mean(), df_3500_QBEP_150.GVF.max()/100., 3500, 150,
                                               (df_3500_QBEP_150.T1F.mean() + df_3500_QBEP_150.T10F.mean())/2)
    fig3 = plt.figure(dpi=300)
    ax3 = fig3.add_subplot(111)
    ax3.scatter(df_3500_QBEP_150.GVF, df_3500_QBEP_150['DP2-3']/ df_3500_QBEP_150['DP2-3'].max(),
                label='Experiment', facecolor='none', edgecolor='b')
    ax3.plot(df_model.gf * 100, df_model.sun, 'm--', label='Sun (2003)')
    ax3.plot(df_model.gf * 100, df_model.barrios, 'k-.', label='Barrios (2007)')
    ax3.plot(df_model.gf * 100, df_model.zhu, 'r-', label='Zhu & Zhang (2015)')

    ax3.set_xlabel(r'$\lambda_{G}$ (%)')
    ax3.set_ylabel(r'$N_{P}$')
    ax3.legend(frameon=False)

    # compare mapping test for example, N = 3500, 1800, Psep = 150 psi, Qgd = 0.03, 0.01
    df_mapping = df_te2700_mapping[(df_te2700_mapping.RPM == 3500) & (df_te2700_mapping.Psep == 150) &
                                    (df_te2700_mapping.Qgd == 0.01)]

    df_model = te2700_case.mapping_performance(df_mapping.QG.mean(), df_mapping.QL_vol.max() + 100, 3500,
                                               150, (df_mapping.T1F.mean() + df_mapping.T10F.mean())/2)

    fig4 = plt.figure(dpi=300)
    ax4 = fig4.add_subplot(111)
    ax4.scatter(df_mapping.QL_vol, df_mapping['DP2-3'], label='Experiment', facecolor='none', edgecolor='b')
    ax4.plot(df_model.ql, df_model.sun, 'm--', label='Sun (2003)')
    ax4.plot(df_model.ql, df_model.barrios, 'k-.', label='Barrios (2007)')
    ax4.plot(df_model.ql, df_model.zhu, 'r-', label='Zhu & Zhang (2015)')

    ax4.set_xlabel(r'$Q_{L}$ (bpd)')
    ax4.set_ylabel(r'$P$ (psi)')
    ax4.set_xlim(0, df_mapping.QL_vol.max() + 300)
    ax4.legend(frameon=False, fontsize=6)

    # plot surging test error analysis at stage 2-3
    fig5 = plt.figure(dpi=300)
    ax5 = fig5.add_subplot(111)
    x_max = np.arange(0, df_te2700_surging['DP2-3'].max() + 1, 0.1)
    df_model = te2700_case.error_analysis(df_te2700_surging.QL, df_te2700_surging.GVF, df_te2700_surging.RPM,
                                          df_te2700_surging.Psep, (df_te2700_surging.T1F + df_te2700_surging.T10F)/2)
    ax5.scatter(df_te2700_surging['DP2-3'], df_model.zhu, facecolor='none', edgecolor='b', label='')
    ax5.plot(x_max, x_max, 'k--', label='perfect match')
    ax5.plot(x_max, 1.25 * x_max, 'r:')
    ax5.plot(x_max, 0.75 * x_max, 'r:')

    ax5.set_xlabel(r'$P_{exp}$ (psi)')
    ax5.set_ylabel(r'$P_{model}$ (psi)')
    ax5.set_xlim(0, df_te2700_surging['DP2-3'].max() + 1)
    ax5.set_ylim(0, df_te2700_surging['DP2-3'].max() + 1)
    ax5.legend(frameon=False)

    # plot mapping test error analysis at stage 2-3
    fig6 = plt.figure(dpi=300)
    ax6 = fig6.add_subplot(111)
    x_max = np.arange(0, df_te2700_mapping['DP2-3'].max() + 1, 0.1)
    df_model = te2700_case.error_analysis(df_te2700_mapping.QL_vol, df_te2700_mapping.GVF, df_te2700_mapping.RPM,
                                          df_te2700_mapping.Psep, (df_te2700_mapping.T1F + df_te2700_mapping.T10F) / 2)
    ax6.scatter(df_te2700_mapping['DP2-3'], df_model.zhu, facecolor='none', edgecolor='b', label='')
    ax6.plot(x_max, x_max, 'k--', label='perfect match')
    ax6.plot(x_max, 1.25 * x_max, 'r:')
    ax6.plot(x_max, 0.75 * x_max, 'r:')

    ax6.set_xlabel(r'$P_{exp}$ (psi)')
    ax6.set_ylabel(r'$P_{model}$ (psi)')
    ax6.set_xlim(0, df_te2700_mapping['DP2-3'].max() + 1)
    ax6.set_ylim(0, df_te2700_mapping['DP2-3'].max() + 1)
    ax6.legend(frameon=False)

    # fig1.show()
    # fig2.show()
    fig3.show()
    fig4.show()
    # fig5.show()
    # fig6.show()
    disconnect_db(conn)


# plot P100 ESP
def vis_p100_plot():
    conn, c = connect_db('ESP.db')
    p100_case = SinglePhaseCompare('P100', conn)

    df_p100 = pd.read_sql_query("SELECT * FROM BHI_P100;", conn)
    df_p100 = df_p100[df_p100.Viscosity_cP != 0]
    df_vis = df_p100.Viscosity_cP.astype(int)
    df_p100.Viscosity_cP = df_vis
    df_3600 = df_p100[(df_p100.RPM == 3600)]
    df_3000 = df_p100[(df_p100.RPM == 3000)]
    df_2400 = df_p100[(df_p100.RPM == 2400)]
    df_all = [df_3600, df_3000, df_2400]

    for index, df in enumerate(df_all):
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(111)
        id = 0

        for visl in df.Viscosity_cP.unique():
            df_model = p100_case.viscous_fluid_performance(3600 - 600 * index, DENW, visl, 0)
            df_exp = df[df.Viscosity_cP == visl]
            df_turzo = turzo_2000('P100', visl, DENW, 3600 - 600 * index)

            # convert psi to ft
            ax.plot(df_model.QL, df_model.DP * psi_to_ft, c='C{}'.format(id + 1), label=r'{} cp_model'.format(visl),
                    linewidth=0.75)
            ax.plot(df_turzo.Qvis, df_turzo.DPvis * psi_to_ft, c='C{}'.format(id + 1),
                    label=r'{} cp_Turzo (2000)'.format(visl), linestyle='--', linewidth=0.75)
            ax.scatter(df_exp['Flow_bpd'], df_exp['Head_ft'], marker=symbols[id], label=r'{} cp_exp'.format(visl),
                        facecolors='none', edgecolor='C{}'.format(id + 1), s=11, linewidths=0.75)

            # plot affinity law
            rpm = df_exp.RPM.unique()
            if visl <= 30:
                dfn = df_3600[df_3600.Viscosity_cP == visl]
                if rpm == 3000 or rpm == 2400:
                    ax.plot(dfn['Flow_bpd'] * rpm / 3600., dfn['Head_ft'] * (rpm / 3600.)**2,
                            c='C{}'.format(id + 1), linestyle='-.', linewidth=0.75, label='')
            id += 1

        handles, labels = ax.get_legend_handles_labels()
        handles = [handles[16], handles[17], handles[18], handles[19], handles[20], handles[21], handles[22], handles[23],
                   handles[0], handles[2], handles[4], handles[6], handles[8], handles[10], handles[12], handles[14],
                   handles[1], handles[3], handles[5], handles[7], handles[9], handles[11], handles[13], handles[15]]

        labels = [labels[16], labels[17], labels[18], labels[19], labels[20], labels[21], labels[22], labels[23],
                   labels[0], labels[2], labels[4], labels[6], labels[8], labels[10], labels[12], labels[14],
                   labels[1], labels[3], labels[5], labels[7], labels[9], labels[11], labels[13], labels[15]]

        ax.set_xlabel(r'$Q_L$ (bpd)')
        ax.set_ylabel(r'$H$ (ft)')
        ax.legend(handles, labels, frameon=False, fontsize=4)
        fig.show()

    # plot best match line curve
    fig2 = plt.figure(dpi=300)
    ax2 = fig2.add_subplot(111)
    model_predict = []
    for i in range(df_p100.shape[0]):
        ql = df_p100.iloc[i].Flow_bpd
        N = df_p100.iloc[i].RPM
        VISL = df_p100.iloc[i].Viscosity_cP
        model_predict.append(p100_case.single_phase_viscous(ql, N, DENW, VISL, 0)*psi_to_ft)

    df_p100['model'] = model_predict
    # plot error less than 25%
    df_p100_selected = df_p100[abs(df_p100.model - df_p100.Head_ft) / df_p100.Head_ft < 0.3]
    x_max = np.array(range(int(max(df_p100.model) + 1)))

    for index, rpm in enumerate(sorted(df_p100_selected.RPM.unique())):
        df_sub = df_p100_selected[df_p100_selected.RPM == rpm]
        ax2.scatter(df_sub.Head_ft, df_sub.model, facecolor='none', edgecolor=colors[index],
                    linewidth=0.75, marker=symbols[index], label=r'$N$={} rpm'.format(rpm))

    ax2.plot(x_max, x_max, 'k--', label='perfect match')
    ax2.plot(x_max, 1.25 * x_max, ls=':', c='gray')
    ax2.plot(x_max, 0.75 * x_max, ls=':', c='gray')
    ax2.set_xlabel(r'$H_{exp}$ (ft)')
    ax2.set_ylabel(r'$H_{model}$ (ft)')
    ax2.legend(frameon=False, fontsize=8)
    fig2.show()

    # epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = \
    #     stats_analysis(df_p100_selected.model, df_p100_selected.Head_ft)
    # print(epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6)

    disconnect_db(conn)


# plot DN1750 ESP
def vis_dn1750_solano():
    conn, c = connect_db('ESP.db')
    dn1750_case = SinglePhaseCompare('DN1750', conn)

    df_dn1750 = pd.read_sql_query("SELECT * FROM DN1750_Solano;", conn)
    df_dn1750 = df_dn1750[df_dn1750.viscosity_cP != 0]

    # extract high and low flow rates > 50 bpd
    dfl = df_dn1750[df_dn1750.Flowrate_bpd < 50]
    dfh = df_dn1750[df_dn1750.Flowrate_bpd > 50]

    # plot best match line curve
    fig1 = plt.figure(dpi=300)
    ax1 = fig1.add_subplot(111)
    model_predict = []

    # random select samples
    dfls = dfl.sample(frac=0.00018, random_state=42)
    dfhs = dfh.sample(frac=0.003, random_state=42)
    df = pd.concat([dfls, dfhs])

    for i in range(df.shape[0]):
        ql = df.iloc[i].Flowrate_bpd
        N = df.iloc[i].N
        DENL = df.iloc[i]['Density_kg/m3']
        VISL = df.iloc[i].viscosity_cP
        model_predict.append(dn1750_case.single_phase_viscous(ql, N, DENL, VISL, 0))

    df['model'] = model_predict
    x_max = np.array(range(21))

    for index, rpm in enumerate(sorted(df['N'].unique())):
        df_sub = df[df.N == rpm]
        ax1.scatter(df_sub.Stage3_DeltaP_psi, df_sub.model, facecolor='none', edgecolor=colors[index],
                    linewidth=0.5, marker=symbols[index], label=r'$N$={} rpm'.format(rpm))

    ax1.plot(x_max, x_max, 'k--', label='perfect match')
    ax1.plot(x_max, 1.25 * x_max, ls=':', c='gray')
    ax1.plot(x_max, 0.7 * x_max, ls=':', c='gray')
    ax1.set_xlabel(r'$P_{exp}$ (psi)')
    ax1.set_ylabel(r'$P_{model}$ (psi)')
    ax1.legend(frameon=False, fontsize=7)
    fig1.show()

    # epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df.model, df.Stage3_DeltaP_psi)
    # print(epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6)

    disconnect_db(conn)


def vis_dn1750_banjar():
    conn, c = connect_db('ESP.db')
    dn1750_case = SinglePhaseCompare('DN1750', conn)

    df_dn1750 = pd.read_sql_query("SELECT * FROM DN1750_Banjar;", conn)
    df_oil_only = df_dn1750[df_dn1750.water_cut == 0]
    df_3500 = df_oil_only[df_oil_only.RPM == 3500]
    df_3000 = df_oil_only[df_oil_only.RPM == 3000]
    df_2500 = df_oil_only[df_oil_only.RPM == 2500]
    df_2000 = df_oil_only[df_oil_only.RPM == 2000]
    df_all = [df_3500, df_3000, df_2500, df_2000]

    for index, df in enumerate(df_all):
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(111)
        id = 0
        af = 0.9        # adjust factor for Banjar data to account for scales

        for visl in np.sort(df.oil_viscosity_cP.unique()):
            df_exp = df[df.oil_viscosity_cP == visl]
            denl = df_exp['Density_Kg/m3'].mean()
            df_model = dn1750_case.viscous_fluid_performance(3500 - index * 500, denl, visl, 0)
            df_turzo = turzo_2000('DN1750', visl, denl, 3500 - index * 500)
            ax.scatter(df_exp['QL_bpd'], df_exp['DPStage3_psi'], marker=symbols[id], label=r'{} cp_exp'.format(visl),
                        facecolors='none', edgecolor='C{}'.format(id + 1), linewidths=0.75)
            ax.plot(df_model.QL, af * df_model.DP, c='C{}'.format(id + 1), label=r'{} cp_model'.format(visl),
                     linewidth=0.75)
            # ax.plot(df_turzo.Qvis, df_turzo.DPvis, c='C{}'.format(id + 1), label=r'{} cp_Turzo (2000)'.format(visl),
            #          linestyle='--', linewidth=0.75)
            id += 1

        handles, labels = ax.get_legend_handles_labels()
        ax.set_xlabel(r'$Q_L$ (bpd)')
        ax.set_ylabel(r'$P$ (psi)')
        ax.legend(handles, labels, frameon=False, fontsize=6)
        # fig.show()

    # plot best match line curve
    fig3 = plt.figure(dpi=300)
    ax3 = fig3.add_subplot(111)
    model_predict = []

    for i in range(df_oil_only.shape[0]):
        ql = df_oil_only.iloc[i].QL_bpd
        N = df_oil_only.iloc[i].RPM
        DENL = df_oil_only.iloc[i]['Density_Kg/m3']
        VISL = df_oil_only.iloc[i].oil_viscosity_cP
        model_predict.append(dn1750_case.single_phase_viscous(ql, N, DENL, VISL, 0))

    df_oil_only['model'] = model_predict
    x_max = np.array(range(int(max(df_oil_only.model) + 1)))

    for index, rpm in enumerate(sorted(df_oil_only['RPM'].unique())):
        df_sub = df_oil_only[df_oil_only.RPM == rpm]
        ax3.scatter(df_sub.DPStage3_psi, af * df_sub.model, facecolor='none', edgecolor=colors[index],
                    linewidth=0.5, marker=symbols[index], label=r'$N$={} rpm'.format(rpm))
    ax3.plot(x_max, x_max, 'k--', label='perfect match')
    ax3.plot(x_max, 1.25 * x_max, ls=':', c='gray')
    ax3.plot(x_max, 0.85 * x_max, ls=':', c='gray')
    ax3.set_xlabel(r'$P_{exp}$ (psi)')
    ax3.set_ylabel(r'$P_{model}$ (psi)')
    ax3.legend(frameon=False, fontsize=8)
    # fig3.show()
    # epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = \
    #     stats_analysis(af * df_oil_only.model, df_oil_only.DPStage3_psi)
    # print(epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6)

    # plot emulsion error line
    df_emulsion = df_dn1750[(df_dn1750.oil_viscosity_cP != 0) & (df_dn1750.water_cut != 0)]
    df_emulsion = df_emulsion.copy()
    fig4 = plt.figure(dpi=300)
    ax4 = fig4.add_subplot(111)
    emulsion_predict = []

    for i in range(df_emulsion.shape[0]):
        ql = df_emulsion.iloc[i].QL_bpd
        N = df_emulsion.iloc[i].RPM
        DENL = df_emulsion.iloc[i]['Density_Kg/m3']
        VISL = df_emulsion.iloc[i].oil_viscosity_cP
        WC = df_emulsion.iloc[i].water_cut * 100
        emulsion_predict.append(dn1750_case.single_phase_viscous(ql, N, DENL, VISL, WC))

    df_emulsion['model'] = emulsion_predict
    x_max = np.array(range(int(max(df_emulsion.DPStage3_psi) + 2)))

    for index, waterCut in enumerate(sorted(df_emulsion['water_cut'].unique())):
        df_sub = df_emulsion[df_emulsion.water_cut == waterCut]
        ax4.scatter(df_sub.DPStage3_psi, af * df_sub.model, facecolor='none', edgecolor=colors[index],
                    linewidth=0.5, marker=symbols[index], label=r'$WC$={}%'.format(int(100 * waterCut)))

    ax4.plot(x_max, x_max, 'k--', label='perfect match')
    ax4.plot(x_max, 1.25 * x_max, 'r:')
    ax4.plot(x_max, 1.50 * x_max, 'r:')
    ax4.plot(x_max, 0.85 * x_max, 'r:')
    ax4.set_xlabel(r'$P_{exp}$ (psi)')
    ax4.set_ylabel(r'$P_{sim}$ (psi)')
    ax4.set_ylim(0, 14)
    ax4.legend(frameon=False, fontsize=7)
    fig4.show()

    # epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = \
    #     stats_analysis(af * df_emulsion.model, df_emulsion.DPStage3_psi)
    # print(epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6)

    disconnect_db(conn)


# plot GC6100 ESP
def gl_GC6100_plot():
    conn, c = connect_db('ESP.db')
    gc6100_case = TwoPhaseCompare('GC6100', conn)
    df_gc6100 = pd.read_sql_query("SELECT * FROM GC6100_Gamboa_All;", conn)

    # only look at positive head
    df_gc6100 = df_gc6100[df_gc6100.DPStage10_psi > 0]

    """
    # compare surging test for example, N = 3000 rpm, Pin = 150 psi, Ql = 3497bpd
    df_surging = df_gc6100[(df_gc6100.test_type == 'Surging') & (df_gc6100.N == 3000) &
                           (df_gc6100.Qstl == 3497) & (df_gc6100.Pin_psi == 150)]
    df_model = gc6100_case.surging_performance(df_surging.QL_bpd.mean(), df_surging['GVF_Stage10_%'].max() / 100.,
                                               3000, 150, df_surging.Temperature_F.mean())
    fig1 = plt.figure(dpi=300)
    ax1 = fig1.add_subplot(111)
    ax1.scatter(df_surging['GVF_Stage10_%'], df_surging.DPStage10_psi / df_surging.DPStage10_psi.max(),
                label='Experiment', facecolor='none', edgecolor='b')
    ax1.plot(df_model.gf * 100, df_model.sun, 'm--', label='Sun (2003)')

    ax1.set_xlabel(r'$\lambda_{G}$ (%)')
    ax1.set_ylabel(r'$N_{P}$')
    ax1.legend(frameon=False, fontsize=6)
    fig1.show()

    # compare mapping test for example
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(111)
    for qgd in [1, 2, 3, 4]:
        try:
            df_mapping = df_gc6100[(df_gc6100.test_type == 'Mapping') & (df_gc6100.N == 2400) &
                                   (df_gc6100.Qgd == qgd) & (df_gc6100.Pin_psi == 250)]

            df_qg = df_mapping['GVF_Stage10_%']/(100-df_mapping['GVF_Stage10_%'])*df_mapping.QL_bpd
            df_model = gc6100_case.mapping_performance(df_qg.mean(), df_mapping.QL_bpd.max() + 100, 2400 + 450,
                                                       250, df_mapping.Temperature_F.mean())

            ax.scatter(df_mapping.QL_bpd, df_mapping.DPStage10_psi, label=r'Experiment $q_{gd}$'+'={}'.format(qgd/100.),
                       marker=symbols[qgd], facecolor='none', edgecolor='C{}'.format(qgd-1), linewidths=0.75, s=8)
            ax.plot(df_model.ql, df_model.sun, linestyle='--', label=r'Model $q_{gd}$'+'={}'.format(qgd/100.),
                    c='C{}'.format(qgd-1), linewidth=0.75)
                    
        except ValueError:
            continue
            
    ax.set_xlabel(r'$Q_{L}$ (bpd)')
    ax.set_ylabel(r'$P$ (psi)')
    ax.set_xlim(0, df_mapping.QL_bpd.max() + 300)
    ax.set_title(r'$N=2400$ rpm, $P=250$ psig', fontsize=10)
    ax.legend(frameon=False, fontsize=4)
    fig.show()
    """

    # plot error analysis
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(111)
    df = df_gc6100[(df_gc6100.test_type == 'Mapping') & (df_gc6100.N == 3000) & (df_gc6100.Qgd <= 4)]
    df_model = gc6100_case.error_analysis(df.QL_bpd, df['GVF_Stage10_%'], df.N + 650, df.Pin_psi, df.Temperature_F)
    ax.scatter(df.DPStage10_psi, df_model.sun, edgecolor='b', s=10, facecolor='none', linewidths=0.5, label='')

    df_stats = pd.DataFrame({'pre': df_model.sun.values.tolist(), 'exp': df.DPStage10_psi.values.tolist()})
    epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
    print(epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6)

    x_max = np.arange(0, plt.gca().get_xlim()[1], 0.2)
    ax.plot(x_max, x_max, 'r--', label='perfect match', linewidth=0.75)

    ax.set_xlabel(r'$P_{exp}$ (psi)')
    ax.set_ylabel(r'$P_{model}$ (psi)')
    ax.set_title(r'$N=3000$ rpm', fontsize=10)
    ax.legend(frameon=False)

    fig.show()

    disconnect_db(conn)


# define Turzo et al. (2000) method for viscous fluid flow
def turzo_2000(pump, vis, den, N):
    """
    :param pump: pump name in string
    :param vis: viscosity in cP
    :param den: density in kg/m3
    :param N: rotational speed in rpm
    :return: boosting pressure at four different flow rates, 0.6, 0.8, 1.0. 1.2 Qbep
    """
    QBEM = {'TE2700': 4500, 'DN1750': 3000, 'GC6100': 7800, 'P100': 11000}
    QBEP = {'TE2700': 2700, 'DN1750': 1750, 'GC6100': 6100, 'P100': 9000}

    qbep = QBEP[pump] * bbl_to_m3 / 3600 / 24
    qbem = QBEM[pump] * bbl_to_m3 / 3600 / 24
    esp = ESP[pump]
    ns = esp['NS']
    sgm = esp['SGM']
    sn = esp['SN']
    st = esp['ST']
    denl = esp['DENL']
    visl = esp['VISL']
    visw = esp['VISW']
    wc = esp['WC']

    sgl = SinglePhaseModel(esp, qbem)
    DPbep, _, _, _, _, _, _, _ = sgl.sgl_calculate_2018(qbep, qbem, denl, DENW, N, ns, sgm, sn, st, visl, visw, wc)
    DPbep06, _, _, _, _, _, _, _ = sgl.sgl_calculate_2018(0.6 * qbep, qbem, denl, DENW, N, ns, sgm, sn, st, visl,
                                                          visw, wc)
    DPbep08, _, _, _, _, _, _, _ = sgl.sgl_calculate_2018(0.8 * qbep, qbem, denl, DENW, N, ns, sgm, sn, st, visl,
                                                          visw, wc)
    DPbep12, _, _, _, _, _, _, _ = sgl.sgl_calculate_2018(1.2 * qbep, qbem, denl, DENW, N, ns, sgm, sn, st, visl,
                                                          visw, wc)
    # convert units
    qbep = QBEP[pump] * 42 / 1440               # to 100 gpm
    vis = vis / (den / DENW)                    # to cSt
    Hbep = DPbep * psi_to_ft

    y = -7.5946 + 6.6504 * np.log(Hbep) + 12.8429 * np.log(qbep)
    Qstar = np.exp((39.5276 + 26.5605 * np.log(vis) - y)/51.6565)

    CQ = 1.0 - 4.0327e-3 * Qstar - 1.724e-4 * Qstar**2

    CH06 = 1.0 - 3.68e-3 * Qstar - 4.36e-5 * Qstar**2
    CH08 = 1.0 - 4.4723e-3 * Qstar - 4.18e-5 * Qstar**2
    CH10 = 1.0 - 7.00763e-3 * Qstar - 1.41e-5 * Qstar**2
    CH12 = 1.0 - 9.01e-3 * Qstar + 1.31e-5 * Qstar**2

    Qvis = CQ * np.array([0.6 * qbep, 0.8 * qbep, qbep, 1.2 * qbep]) * 1440 / 42        # to bpd
    DPvis = np.array([CH06, CH08, CH10, CH12]) * np.array([DPbep06, DPbep08, DPbep, DPbep12])

    df = pd.DataFrame({'Qvis': Qvis, 'DPvis': DPvis, 'N': [N] * 4})

    return df


# define a statistical function
def stats_analysis(df_pre, df_exp):
    df_relative = (df_pre - df_exp) / df_exp * 100
    df_actual = df_pre - df_exp

    epsilon1 = df_relative.mean()
    epsilon2 = df_relative.abs().mean()
    epsilon3 = df_relative.std()
    epsilon4 = df_actual.mean()
    epsilon5 = df_actual.abs().mean()
    epsilon6 = df_actual.std()
    return epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6


# def a function to plot inversion point curves
def plot_inversion(pump):
    QBEM_all = {'TE2700': 4500, 'DN1750': 3000, 'GC6100': 7800, 'P100': 11000}
    QBEM = QBEM_all[pump]

    # esp model
    esp = ESP[pump]
    esp_sgl = SinglePhaseModel(esp, QBEM)
    emulsion = np.vectorize(esp_sgl.emulsion)

    # input array
    WC = np.arange(0, 100, 1)
    VOI = esp['VOI'] * np.ones(WC.shape)
    R2 = esp['R2'] * np.ones(WC.shape)
    VISO = 0.01 * np.ones(WC.shape)
    VISW = 0.001 * np.ones(WC.shape)
    DENO = 850 * np.ones(WC.shape)
    DENW = 1000 * np.ones(WC.shape)
    ST = 0.02 * np.ones(WC.shape)
    SN = 1 * np.ones(WC.shape)
    N = 3600 * np.ones(WC.shape)
    Q = 1500 * bbl_to_m3 / 24.0 / 3600.0 * np.ones(WC.shape)
    mod = np.empty(WC.shape, dtype='str')
    mod[:] = 'zhu'

    fig1 = plt.figure(dpi=300)
    ax1 = fig1.add_subplot(111)
    ax1.plot(WC, emulsion(VOI, R2, VISO, VISW, DENO, DENW, WC, ST, N, Q, SN, mod) * 1000, 'r-', label=r'$Q_{L} = 1500$ bpd')
    ax1.plot(WC, emulsion(VOI, R2, VISO, VISW, DENO, DENW, WC, ST, N, Q/3., SN, mod) * 1000, 'b-', label=r'$Q_{L} = 500$ bpd')
    ax1.plot(WC, emulsion(VOI, R2, VISO, VISW, DENO, DENW, WC, ST, N, Q/10., SN, mod) * 1000, 'k-', label=r'$Q_{L} = 150$ bpd')
    ax1.set_xlabel(r'Water cut (%)')
    ax1.set_ylabel(r'Effective viscosity (cP)')
    ax1.legend(frameon=False)
    # fig1.show()

    fig2 = plt.figure(dpi=300)
    ax2 = fig2.add_subplot(111)
    ax2.plot(WC, emulsion(VOI, R2, VISO, VISW, DENO, DENW, WC, ST, 1.5*N, Q, SN, mod) * 1000, 'r-', label=r'$N = 4800$ rpm')
    ax2.plot(WC, emulsion(VOI, R2, VISO, VISW, DENO, DENW, WC, ST, N, Q, SN, mod) * 1000, 'b-', label=r'$N = 3600$ rpm')
    ax2.plot(WC, emulsion(VOI, R2, VISO, VISW, DENO, DENW, WC, ST, N/1.5, Q, SN, mod) * 1000, 'k-', label=r'$N = 2400$ rpm')
    ax2.set_xlabel(r'Water cut (%)')
    ax2.set_ylabel(r'Effective viscosity (cP)')
    ax2.legend(frameon=False)
    # fig2.show()

    fig3 = plt.figure(dpi=300)
    ax3 = fig3.add_subplot(111)
    ax3.plot(WC, emulsion(VOI, R2, VISO, VISW, DENO, DENW, WC, ST, N, Q, SN, mod) / (VISO[1]), 'r-',
             label=r'$\mu_{oil} = 10$ cP')
    ax3.plot(WC, emulsion(VOI, R2, 5*VISO, VISW, DENO, DENW, WC, ST, N, Q, SN, mod) / (5 * VISO[1]), 'b-',
             label=r'$\mu_{oil} = 50$ cP')
    ax3.plot(WC, emulsion(VOI, R2, 15*VISO, VISW, DENO, DENW, WC, ST, N, Q, SN, mod) / (15 * VISO[1]), 'k-',
             label=r'$\mu_{oil} = 150$ cP')
    ax3.set_xlabel(r'Water cut (%)')
    ax3.set_ylabel(r'Relative viscosity ($\frac{\mu_{eff}}{\mu_{oil}}$)')
    ax3.legend(frameon=False)
    fig3.show()


##########################
if __name__ == "__main__":

    # connect database
    # conn, c = connect_db('ESP.db')

    # te2700_case = SinglePhaseCompare('TE2700', conn)
    # te2700_case.single_phase_water()
    # dn1750_case = SinglePhaseCompare('DN1750', conn)
    # dn1750_case.single_phase_water()
    # gc6100_case = SinglePhaseCompare('GC6100', conn)
    # gc6100_case.single_phase_water()
    # p100_case = SinglePhaseCompare('P100', conn)
    # p100_case.single_phase_water()

    # vis_te2700_plot()
    # vis_p100_plot()
    vis_dn1750_banjar()
    # vis_dn1750_solano()

    # te2700_gl = TwoPhaseCompare('TE2700', conn)
    # dfgv, dfcd = te2700_gl.surging_cd_to_db()

    # gl_te2700_plot()
    # gl_GC6100_plot()

    # disconnect_db(conn)

    # plot_inversion('DN1750')




