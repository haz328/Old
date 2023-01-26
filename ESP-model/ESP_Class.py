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
from ESP_inputs import ESP, QBEM_default

from datetime import datetime
import time
from mpl_toolkits.mplot3d import Axes3D
from itertools import groupby


# global variables
factor = 1.2 # factor used in gas bubble size
DB_GVF_EFF = 1.02    ## factor used in gas bubble size
DB_NS_EFF = 2    ## factor used in gas bubble size

FTI = 3       # Original turning factor
FTD = 3       # Original turning factor
CD_Gas_EFF = 0.4    # Original drag coefficiet for gas velocity in modified Sun for INT flow
CD_Liquid_EFF = 0.4    # Original drag coefficiet for gas velocity in modified Sun for INT flow
CD_INT_EFF = 2
CD_GV_EFF = 2
error_control_high = 100  # reduce data noise for error > error_control_high
error_control_low = -50  # reduce data noise for error < error_control_low
ABSerror_control_high = 2  # reduce data noise for error > error_control_high
ABSerror_control_low = 2  # reduce data noise for error > error_control_high
transition_zone = 0.3   # transition zone between slug and bubble, transition_zone*QBEM
G = 9.81
pi = np.pi
E1 = 1e-3
DENW = 997.                 # water density
VISW = 1e-3                 # water viscosity
psi_to_pa = 1.013e5 / 14.7
psi_to_ft = 2.3066587368787
bbl_to_m3 = 0.15897
colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
symbols = ['o', 's', '^', '*', 'd', 'p', 'v', 'D', 'x', '+']
sgl_model = 'zhang_2016'
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
mpl.rcParams['text.color'] = 'black'
mpl.rcParams['axes.labelcolor'] = 'black'
mpl.rcParams['lines.linewidth'] = 0.75

plt.style.use('seaborn-ticks')
mpl.rcParams['figure.figsize'] = (4, 3)
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.color'] = 'black'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.color'] = 'black'
mpl.rcParams['markers.fillstyle'] = 'none'
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams["font.weight"] = "bold"
mpl.rcParams["axes.labelweight"] = "bold"
mpl.rcParams["axes.labelcolor"] = "black"
mpl.rcParams["axes.titleweight"] = "bold"
# mpl.rcParams["axes.titlecolor"] = "black"    
mpl.rcParams["axes.linewidth"] = 0.75

#xtick.labelsize:     medium  # fontsize of the tick labels


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
    def get_fff_leakage(re, ed, N, RLK, VLK, LG, SL):
        '''
        N: rotational speed RPM
        RLK: leakage radius (rotational effect, assume equal to RI=R1+R2)
        VLK: axial velocity in the leakage area
        LG: leakage length
        SL: leakage width



        '''
        # friction factor based on Childs 1983 and Zhu et al. 2019 10.4043/29480-MS
        REM = re
        OMEGA = 2.0 * pi * N / 60.0     # rotational speed in rads/s
        fff = LG/SL*0.066*REM**-0.25*(1+OMEGA**2*RLK**2/4/VLK**2)**0.375
        return fff
        
    def C_Drag_Sphere(self, ed, DI, VISL, DENL, VL, VSR, DB):
        # DI:               Flow passage hydraulic diameter
        # DENL:             Liquid density
        # VISL:             Liquid viscosity
        # VL:               Liquid velocity
        # VSR:              Slip velocity between gas and
        # n:                FLuid behavior index
        # K:                Non newtonian fluid viscosity, n and K are for non-newtonian fluid
        # CD:               Drag coefficient
        # CDU:              Coefficient to get CD by Graham & Jones
        # DB:               Bubble size

        # Get eta
        
        # Assume newtonian fluid
        n=1
        K=VISL
        Re_bed = DI ** n * DENL * VL**2/K     # Fluid velocity above sand bed (assume no sand bed in this equation, detail please refer the class notes)
        f_bed = self.get_fff(Re_bed, ed)
        if Re_bed == 0:
            f_bed = 0
        Tau_bed = f_bed * DENL * VL**2 / 2      # Average bed shear stess
        dUdy = (Tau_bed / K)**(1 / n) # Velocity gradient
        eta = DB/VL * dUdy

        # Reynold number around particle
        Rep = (DENL * VSR**(2-n) * DB**n)/K

        # CDU
        if Rep < 0.2 * 2**n:
            CDU = 24 * 2**(1-n) / Rep * (2-n)
        elif Rep < 24 * 2**n:
            CDU = 35.2 / (Rep / 2**n)**1.03 + n * (1 - 20.9 / (Rep / 2**n)**1.11)
        elif Rep < 100 * 2**n:
            CDU = 37 / (Rep / 2**n)**1.11 + 0.25 + 0.36 * n
        else:
            pass
        
        # corrected CD
        if Rep < 100 * 2**n:
            CD = 0.8 * CDU * (1 + (5e-4 * Rep + 0.0179)*eta)   # not sure
        else:
            Rep = 100 * 2**n
            CDU = 37 / (Rep / 2**n)**1.11 + 0.25 + 0.36 * n
            CD = 0.8 * CDU * (1 + (5e-4 * Rep + 0.0179)*eta)        # From drilling fluid ppt notes, might be improved in the future
        
        return CD

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

    def sgl_calculate_old(self, Q, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC):
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
            # FTI = 2.5
            # FTD = 2.5
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
        HLKloss = 0.
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
        if WC==100:
            VISL=VISW
        elif WC > 0.:
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
            #HFI = 2.5 * 4.0 * FFI * (W1 + W2) ** 2 * self.LI / (8.0 * G * DI)
            HFI = 2.5 * 4.0 * FFI * (VI) ** 2 * self.LI / (2.0 * G * DI)    # modified by Haiwen zhu (it is used in gl model)
            HFD = 2.5 * 4.0 * FFD * VD ** 2 * self.LD / (2.0 * G * DD)

            # turn loss
            # FTI = 3.0
            # FTD = 3.0
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
                FFL = self.get_fff_leakage(REL,EDL,N, (self.R1+self.R2)/2, VL, self.LG, self.SL)        # by Haiwen Zhu
                VL = np.sqrt(2.0 * G * HLK / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLK_new = 2.0 * pi * self.RLK * self.SL * VL
            else:
                VL = np.abs(QLK / (2.0 * pi * self.RLK * self.SL))
                REL = DENL * VL * self.SL / VISL
                EDL = 0.
                FFL = self.get_fff(REL, EDL)
                FFL = self.get_fff_leakage(REL,EDL,N, (self.R1+self.R2)/2, VL, self.LG, self.SL)        # by Haiwen Zhu
                VL = np.sqrt(2.0 * G * np.abs(HLK) / (1.5 + 4.0 * FFL * self.LG / self.SL))
                QLK_new = -2.0 * pi * self.RLK * self.SL * VL

            ABQ = np.abs((QLK_new - QLK) / QLK_new)
            QLK = QLK_new
            HLKloss = 20/2/G*(QLK/self.ZD/AI)**2        # by Haiwen Zhu
            HP_new -= HLKloss
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

    def sgl_calculate_jc(self, Q, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC):
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
        if WC==100:
            VISL=VISW
        elif WC > 0.:
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
            # FTI = 3.0
            # FTD = 3.0
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

    def performance_curve(self,sgl_model):
        # sgl_model, select single phase model
                #   'zhang_2015': sgl_calculate_new
                #   'zhang_2016': sgl_calculate_old
                #   'jiecheng_2017': sgl_calculate_jc
                #   'zhu_2018': sgl_calculate_2018
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
            if sgl_model=='zhang_2015':
                HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_old(ql, self.QBEM, self.DENL, DENW, self.N, self.NS,
                                                                        self.SGM, self.SN, self.ST, self.VISL, VISW,
                                                                        self.WC)              
            elif sgl_model=='zhang_2016':
                HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_new(ql, self.QBEM, self.DENL, DENW, self.N, self.NS,
                                                                        self.SGM, self.SN, self.ST, self.VISL, VISW,
                                                                        self.WC)              
            elif sgl_model=='jiecheng_2017':
                HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_jc(ql, self.QBEM, self.DENL, DENW, self.N, self.NS,
                                                                        self.SGM, self.SN, self.ST, self.VISL, VISW,
                                                                        self.WC)              
            elif sgl_model=='zhu_2018':
                HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_2018(ql, self.QBEM, self.DENL, DENW, self.N, self.NS,
                                                                        self.SGM, self.SN, self.ST, self.VISL, VISW,
                                                                        self.WC)                #new, 2018, jc are available. old need to change the inputs

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
    
    @staticmethod
    def CD_Cal(VSR, DENL, VISL, N, DB):
        # CD-Legendre & Magnaudet (1998), Clift et al. (1978), Rastello et al. (2011)
        REB = DENL * VSR * DB / VISL
        SR = DB * (2.0 * np.pi * N / 60.0) / VSR


        if REB <= 50:
            # if REB < 0.001:
            #     REB = 0.001
            # if SR < 0.001:
            #     SR = 0.001
            CD = 24.0 / REB * (1.0 + 0.15 * REB**0.687) * (1.0 + 0.3 * SR**2.5)
        else:
            CD = 24.0 / REB * (1.0 + 0.15 * REB**0.687) * (1.0 + 0.55 * SR**2.0)
    
        # #def CD_Barrios(VSR, DENL, VISL, N, DB):
        # REB = DENL * VSR * DB / VISL
        # # Barrios
        # YY = 0.00983 + 389.9 * REB / N ** 2
        # CD = (24.0 + 5.48 * (REB * YY) ** 0.427 + 0.36 * REB * YY) / (REB * YY)


        # YY = 0.00983 + 389.9 * REB / N ** 2
        # CD = (24.0 + 5.48 * (REB * YY) ** 0.427 + 0.36 * REB * YY) / (REB * YY)

        # # From drilling fluid, and Zhu's ABM report 2020 April
        # VL = Q/ self.ZI/AI/LambdaC2
        # CD = self.C_Drag_Sphere(EDI, DI, self.VISL, self.DENL, VL, VSR, DB)

        return CD

    def get_DB(self, Q, HP, N, GV):
        # use DB or DBMax to calculate bubble flow?
        # factor = 0.8/0.43       # 0.8 good for Flex31, TE2700
        
        # factor = 1.5       # 0.8 good for GC6100, TE2700

        # factor = 0.7      # 
        # DB = factor*6.034 * GV**DB_GVF_EFF * (self.ST / self.DENL)**0.6 * (HP * Q * (N / self.NS * 0.05/self.R2)**3 / (self.DENL * self.ZI * self.VOI))**(-0.4) * (self.DENL / self.DENG)**0.2  #original
        DB = factor*6.034 * GV**DB_GVF_EFF * (self.ST / self.DENL)**0.6 * (HP * Q * (N / 3500 * 0.05/self.R2)**DB_NS_EFF / (self.DENL * self.ZI * self.VOI))**(-0.4) * (self.DENL / self.DENG)**0.2  #original
        
        # DB = factor**6.034 * GV**1.2 * (self.ST / self.DENL)**0.6 * (HP * Q  / (self.DENL * self.ZI * self.VOI))**(-0.4) * (self.DENL / self.DENG)**0.2  #original
        # DB = DB*((self.NS/N*self.R2/0.05))**0.4       # 0.05 reference R2 of GC6100
        DBMAX = DB/0.43
        # DB=factor*1.4*(GV)**0.25*(self.ST/self.DENL)**0.6*(HP*Q/(self.DENL*self.ZI*self.VOI))**(-0.4)* \
        #                 (self.DENL/self.DENG)**0.2     # Bubble size-Zhu & Zhang, 2015
        # DBMAX = DB/0.43
        # DBMAX = DB/0.6
        return DB, DBMAX
    
    def get_lambda_c1(self, Q, HP, N):
        # critical bubble size in turbulent flow (Barnea 1982)
        DBCrit = 2.0 * (0.4 * 0.073 / (self.DENL - self.DENG) / ((2.0 * pi * N / 60.0)**2 * self.RI))**0.5
        LambdaC1 = DBCrit / (6.034 / 0.6 * (self.ST / self.DENL)**0.6 * (HP * Q / (self.DENL * self.ZI * self.VOI)) **  \
                             (-0.4) * (self.DENL / self.DENG)**0.2)
        return LambdaC1

    def get_lambda_c2(self, Q, HP, N):
        # rotary speed effect
        # alphaG_crit = pi / 6.0 - (pi / 6.0 - 0.25) * np.exp(-(N / 3600.0)**4)
        alphaG_crit = 0.5-0.25*np.exp(-(N/3500.0)**4.0)   #Fotran code selection
        alphaG = alphaG_crit
        VSR = 0.1
        LambdaC2 = 0.9
        ABV = 1.0
        icon = 0
        relax = 0.1

        AIW = self.AB + self.ASF + self.ASB
        DI = 4.0 * self.VOI / AIW
        EDI = self.EA / DI
        AI = self.VOI / self.LI

        while ABV > 0.00001:
            if icon > 1000:
                break
            else:
                icon += 1

            # DB = 6.034 * LambdaC2 * (self.ST / self.DENL)**0.6 * (HP * Q / (self.DENL * self.ZI * self.VOI))**(-0.4) * \
            #      (self.DENL / self.DENG)**0.2  #original
            # DB=1.4*(LambdaC2)**0.25*(self.ST/self.DENL)**0.6*(HP*Q/(self.DENL*self.ZI*self.VOI))**(-0.4)* \
                        # (self.DENL/self.DENG)**0.2     # Bubble size-Zhu & Zhang, 2015
            # DBMAX = DB / 0.6        #original
            # DBMAX = DB / 0.42      # Bubble size-Zhu & Zhang, 2015
            DB, DBMAX = self.get_DB(Q,HP,N,LambdaC2)
            # DB, DBMAX = self.get_DB(Q,HP,N,alphaG)
            # DBMAX *=1
            REB = self.DENL * VSR * DB / self.VISL
            SR = DB * (2.0 * pi * N / 60.0) / VSR


            # Original CD and bubble size is best for LambdaC2 prediction of GC6100

            # original CD
            if REB <= 50:
                if REB < 0.001:
                    REB = 0.001
                if SR < 0.001:
                    SR = 0.001
                CD = 24.0 / REB * (1.0 + 0.15 * REB**0.687) * (1.0 + 0.3 * SR**2.5)
            else:
                CD = 24.0 / REB * (1.0 + 0.15 * REB**0.687) * (1.0 + 0.55 * SR**2.0)

            CD = self.CD_Cal(VSR, self.DENL, self.VISL, N, DB)

            # # Barrios
            # YY = 0.00983 + 389.9 * REB / N ** 2
            # CD = (24.0 + 5.48 * (REB * YY) ** 0.427 + 0.36 * REB * YY) / (REB * YY)
            # YY = 0.00983 + 389.9 * REB / N ** 2
            # CD = (24.0 + 5.48 * (REB * YY) ** 0.427 + 0.36 * REB * YY) / (REB * YY)

            # # From drilling fluid, and Zhu's ABM report 2020 April
            # VL = Q/ self.ZI/AI/LambdaC2
            # CD = self.C_Drag_Sphere(EDI, DI, self.VISL, self.DENL, VL, VSR, DB)

            VSR = np.sqrt(4.0 * DB * (self.DENL - self.DENG) * self.RI / (3.0 * CD * self.DENL)) * (2.0 * pi * N / 60.0)
            RS = 2*VSR * (2.0 * pi * self.RI - self.ZI * self.TB) * self.YI / Q
            alphaG = (RS - 1.0 + np.sqrt((1.0 - RS)**2 + 4.0 * RS * LambdaC2)) / (2.0 * RS)

            # VSR1 = np.sqrt(4.0 * DB * (self.DENL - self.DENG) * self.R1 / (3.0 * CD * self.DENL)) * (2.0 * pi * N / 60.0)
            # RS1 = 2*VSR1 * (2.0 * pi * self.R1 - self.ZI * self.TB) * self.YI / Q
            # alphaG = (RS1 - 1.0 + np.sqrt((1.0 - RS1)**2 + 4.0 * RS1 * LambdaC2)) / (2.0 * RS1)

            # VSR2 = np.sqrt(4.0 * DB * (self.DENL - self.DENG) * self.R2 / (3.0 * CD * self.DENL)) * (2.0 * pi * N / 60.0)
            # RS2 = 2*VSR2 * (2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI / Q
            # alphaG = (RS2 - 1.0 + np.sqrt((1.0 - RS2)**2 + 4.0 * RS2 * LambdaC2)) / (2.0 * RS2)

            ABV = np.abs((alphaG - alphaG_crit) / alphaG_crit)
            if ABV < E1:
                break
            if alphaG > alphaG_crit:
                # LambdaC2 -= relax
                LambdaC2 *=0.9
            else:
                # LambdaC2 += relax
                # relax /= 5
                LambdaC2 *=1.1

            if LambdaC2 < 0:
                # endtime=time.time()
                # print('LambdaC2 time: ', starttime-endtime, 'icon: ', icon, 'ABV: ', ABV ,'\n')
                LambdaC2 = -LambdaC2
                return LambdaC2

            if icon > 1e5:
                # endtime=time.time()
                # print('LambdaC2 time: ', starttime-endtime, 'icon: ', icon, 'ABV: ', ABV ,'\n')
                return LambdaC2
            else:
                icon += 1
        # endtime=time.time()
        # print('LambdaC2 time: ', starttime-endtime, 'icon: ', icon, 'ABV: ', ABV ,'\n')

        return LambdaC2

    def get_lambda_c3(self, Q, N):
        LambdaC3 = 0.1
        ANG = -90.0 * pi / 180.0
        FIC = 0.0142
        AI = self.VOI / self.LI
        BI = (self.B1 + self.B2) / 2.0
        AIW = self.AB + self.ASF + self.ASB
        DI = 4.0 * self.VOI / AIW
        EDI = self.EA / DI
        GC = (2.0 * pi * N / 60.0)**2 * self.RI * np.sin(BI)
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
        # starttime=time.time()
        # while ABU > E1 or ABHLF > E1:
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
            #elif REMX > 4000.0: #Fotran code
            #    VAV = 1.2*VM    #Fotran code
                VAV = 1.3 * VM
            else:
                VAV = (2.0 - 0.7 * (REMX - 2000.0) / 2000.0) * VM
            #    VAV = (2.0 - 0.8 * (REMX - 2000.0) / 2000.0) * VM  #Fotran code

            VT = VAV + (0.54 * np.cos(ANG) + 0.35 * np.sin(ANG)) * np.sqrt(GC * DI * np.abs(self.DENL - self.DENG) / self.DENL)

            if VT < 0:
                VT = np.abs(VT)
                # VT=0.1     #Fotran code

            HLFN = ((HLS * (VT - VM) + VSL) * (VSG + VSL * FE) - VT * VSL * FE) / (VT * VSG)

            if HLFN < 0.0:
                HLFN = np.abs(HLFN)

            if HLFN > 1.0:
                HLFN = 1.0 / HLFN

            ABHLF = np.abs(HLF-HLFN)/HLFN
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
            THFO = THF * np.sqrt(np.abs(SHI / self.DENG)) / self.VISG
            FRA = FS * (1.0 + 13.8 * (THFO - 200.0 * np.sqrt(self.DENG / self.DENL)) * WED**0.2 / REG**0.6)
            FRA1 = self.get_fff(REC1, THF / DI)

            if FRA > FRA1:
                FRA = FRA1
            # elif FRA < 0:
            #     FRA = FS
            FIN = FRA

            if FIN > 0.5:       #Commentted in Fotran code
                FIN = 0.5

            FI = (FIN + 9.0 * FI) / 10.0
            ABCDN = (-(self.DENL - DENC) * GC + SF * FF * self.DENL * VF * VF / (2.0 * AF)) * 2 / (SI * FI * self.DENG * (1.0 / AF + 1.0 / AC))   # Original
            # ABCDN = - (-(self.DENL - DENC) * GC + SF * FF * self.DENL * VF * VF / (2.0 * AF)) * 2 / (SI * FI * self.DENG * (1.0 / AF + 1.0 / AC))   # by Haiwen
            ABCD = (ABCDN + 9.0 * ABCD) / 10.0

            if ABCD < 0.:
                VCN = VC * 0.9
            else:
                VCN = np.sqrt(ABCD) + VF
                #VCN = np.sqrt(ABCD)     #Fotran code
            if VCN < V24:
                VCN = V24

            ABU = np.abs((VCN - VC) / VC)
            VC = 0.5*VCN+0.5*VC
            VSG = VC * (1.0 - HLF) - VSL * FE
            VSG = VC * (1.0 - HLF) * (1 - FE)

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
        
        # endtime=time.time()
        # print('LambdaC3: ', LambdaC3, 'icon: ',icon, 'ABU: ', ABU, 'VC', VC, 'V24', V24, 'HLS', HLS, 'VSG', VSG, 'HLF', HLF, 'FE', FE,'\n')
        return LambdaC3

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
        factor = 1
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
                CD_to_rb1 = 12. * mium1 * (factor*4.564e7 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)
                # CD_to_rb1 = 12. * mium1 * (4.564e7 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENL)
            elif (GV1 > 0.15) and (GV1 < 0.5):
                CD_to_rb1 = 12. * mium1 * (factor*6.39e7 * (1 - GV1) ** 1.5 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)
                # CD_to_rb1 = 12. * mium1 * (6.39e7 * (1 - GV1) ** 1.5 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENL)
            else:
                CD_to_rb1 = 12. * mium1 * (factor*9.13e7 * GV1 ** 2 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)
                # CD_to_rb1 = 12. * mium1 * (9.13e7 * GV1 ** 2 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENL)

            VSRN1 = np.sqrt(4.0 * (DENL - DENG) * R1 / (3.0 * CD_to_rb1 * DENL)) * OMEGA
            ABV1 = np.abs((VSRN1 - VSR1) / VSR1)
            VSR1 = VSRN1
            RS1 = VSR1 * (2.0 * pi * R1 - ZI * TB) * YI / (Q + QLK)

            if GV2 <= 0.15:
                CD_to_rb2 = 12. * mium2 * (factor*4.564e7 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)
                # CD_to_rb2 = 12. * mium2 * (4.564e7 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENL)
            elif (GV2 > 0.15) and (GV2 < 0.5):
                CD_to_rb2 = 12. * mium2 * (factor*6.39e7 * (1 - GV2) ** 1.5 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)
                # CD_to_rb2 = 12. * mium2 * (6.39e7 * (1 - GV2) ** 1.5 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENL)
            else:
                CD_to_rb2 = 12. * mium2 * (factor*9.13e7 * GV2 ** 2 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)
                # CD_to_rb2 = 12. * mium2 * (9.13e7 * GV2 ** 2 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENL)

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
    def gv_sun_ANN(B1, B2, DENL, DENG, N, Q, QG, QLK, R1, R2, TB, VISL, VISG, YI1, YI2, YI, ZI):
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


            CD_to_rb1 = 12. * mium1 * (9.13e7 * GV1 ** 2 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)

            VSRN1 = np.sqrt(4.0 * (DENL - DENG) * R1 / (3.0 * CD_to_rb1 * DENL)) * OMEGA
            ABV1 = np.abs((VSRN1 - VSR1) / VSR1)
            VSR1 = VSRN1
            RS1 = VSR1 * (2.0 * pi * R1 - ZI * TB) * YI / (Q + QLK)

            CD_to_rb2 = 12. * mium2 * (9.13e7 * GV2 ** 2 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)

            VSRN2 = np.sqrt(4.0 * (DENL - DENG) * R2 / (3.0 * CD_to_rb2 * DENL)) * OMEGA
            ABV2 = np.abs((VSRN2 - VSR2) / VSR2)
            VSR2 = VSRN2
            RS2 = VSR2 * (2.0 * pi * R2 - ZI * TB) * YI / (Q + QLK)

            if GV1 <= GF:
                GV1 = GF
            else:
                GV1 = (RS1 - 1.0 + np.sqrt((1.0 - RS1) ** 2 + 4.0 * RS1 * GF)) / (2.0 * RS1)

            if GV2 <= GF:
                GV2 = GF
            else:
                GV2 = (RS2 - 1.0 + np.sqrt((1.0 - RS2) ** 2 + 4.0 * RS2 * GF)) / (2.0 * RS2)

            GV = (GV1 + GV2) / 2
            CD_to_rb = (CD_to_rb1 + CD_to_rb2) / 2

        return GV, CD_to_rb
    
    @staticmethod
    def gv_sun_INT(B1, B2, DENL, DENG, N, Q, QG, QLK, R1, R2, TB, VISL, VISG, YI1, YI2, YI, ZI):
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
        GV1, CD_to_rb = 0.8, 0.3 
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

            # if GV1 <= 0.15:
            #     CD_to_rb1 = 12. * mium1 * (4.564e7 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)
            #     # CD_to_rb1 = 12. * mium1 * (4.564e7 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENL)
            # elif (GV1 > 0.15) and (GV1 < 0.5):
            # CD_to_rb1 = 12. * mium1 * (6*6.39e7 * (1 - GV1) ** 1.5 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENM1)
            #     # CD_to_rb1 = 12. * mium1 * (6.39e7 * (1 - GV1) ** 1.5 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENL)
            # else:
            CD_to_rb1 = 12. * mium1 * (CD_INT_EFF*9.13e7 * GV1 ** 2 / W1_G ** CD_Gas_EFF / W1_L ** CD_Liquid_EFF) / (np.abs(W1G - W1L) * DENM1) * (N/3600)**2
                # CD_to_rb1 = 12. * mium1 * (9.13e7 * GV1 ** 2 / W1_G ** 0.4) / (np.abs(W1G - W1L) * DENL)

            VSRN1 = np.sqrt(4.0 * (DENL - DENG) * R1 / (3.0 * CD_to_rb1 * DENL)) * OMEGA
            ABV1 = np.abs((VSRN1 - VSR1) / VSR1)
            VSR1 = VSRN1
            RS1 = VSR1 * (2.0 * pi * R1 - ZI * TB) * YI / (Q + QLK)

            # if GV2 <= 0.15:
            #     CD_to_rb2 = 12. * mium2 * (4.564e7 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)
            #     # CD_to_rb2 = 12. * mium2 * (4.564e7 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENL)
            # elif (GV2 > 0.15) and (GV2 < 0.5):
            # CD_to_rb2 = 12. * mium2 * (6*6.39e7 * (1 - GV2) ** 1.5 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENM2)
            #     # CD_to_rb2 = 12. * mium2 * (6.39e7 * (1 - GV2) ** 1.5 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENL)
            # else:
            CD_to_rb2 = 12. * mium2 * (CD_INT_EFF*9.13e7 * GV2 ** 2 / W2_G ** CD_Gas_EFF) / (np.abs(W2G - W2L) * DENM2) * (N/3600)**2
                # CD_to_rb2 = 12. * mium2 * (9.13e7 * GV2 ** 2 / W2_G ** 0.4) / (np.abs(W2G - W2L) * DENL)
            VSRN2 = np.sqrt(4.0 * (DENL - DENG) * R2 / (3.0 * CD_to_rb2 * DENL)) * OMEGA
            ABV2 = np.abs((VSRN2 - VSR2) / VSR2)
            VSR2 = VSRN2
            RS2 = VSR2 * (2.0 * pi * R2 - ZI * TB) * YI / (Q + QLK)
            
            if GV1 <= GF:
                GV1 = GF
            else:
                GV1 = (RS1 - 1.0 + np.sqrt((1.0 - RS1) ** 2 + 4.0 * RS1 * GF)) / (2.0 * RS1)

            if GV2 <= GF:
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

    def gv_Barrios(self, DENL, DENG, GF, N, Q, QLK, R1, RI, ST, TB, VISL, YI, ZI):
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


        AIW = self.AB + self.ASF + self.ASB
        DI = 4.0 * self.VOI / AIW
        EDI = self.EA / DI
        AI = self.VOI / self.LI


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

    @staticmethod # with this statement, we can use self variable and function
    def DBFLOW(GF):
        GV = GF
        return GV

    def gv_zhu(self, DENL, DENG, GF, HP, N, Q, QLK, RI, ST, TB, VISL, VOI, YI, ZI):
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

            # DB = 6.034 * GF * (ST / DENL) ** 0.6 * (HP * (Q + QLK) / (DENL * ZI * VOI)) ** (-0.4) * (DENL / DENG) ** 0.2
            # DBMAX = DB / 0.4
            DB, DBMAX = self.get_DB(Q+QLK, HP,N, GF)
            # REB = DENL * VSR * DBMAX / VISL
            # SR = DBMAX * OMEGA / VSR
            REB = DENL * VSR * DB / VISL
            SR = DB * OMEGA / VSR

            if REB < 50.0:
                CD = 24.0 / REB * (1.0 + 0.15 * REB ** 0.687) * (1.0 + 0.3 * SR ** 2.5)
            else:
                CD = 24.0 / REB * (1.0 + 0.15 * REB ** 0.687) * (1.0 + 0.55 * SR ** 2.0)


            VSRN = np.sqrt(4.0 * DB * (DENL - DENG) * RI / (3.0 * CD * DENL)) * OMEGA
            ABV = np.abs((VSRN - VSR) / VSR)
            VSR = VSRN
            RS = VSR * (2.0 * pi * RI - ZI * TB) * YI / (Q + QLK)

            if GV < 0.0:
                GV = 0.0
            else:
                GV = (RS - 1.0 + np.sqrt((1.0 - RS) ** 2 + 4.0 * RS * GF)) / (2.0 * RS)
                # print(GV)
                
        return GV, CD/DBMAX

    def BUFLOW(self,DENL,DENG,GF,HP,N,Q,RI,ST,VISL,VISG,VOI,YI,ZI,QLK,TB):
        """
        # Bubble flow model based on Zhu (2017)
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
        VSR = 0.1
        ABV = 1.0
        GV  = 0.5
        REB = 0.0
        RS = 0.0
        SR = 0.0
        VSRN = 0.0
        O = 2.0 * np.pi * N / 60.0
        icon=0

        AIW = self.AB + self.ASF + self.ASB
        DI = 4.0 * self.VOI / AIW
        EDI = self.EA / DI
        AI = self.VOI / self.LI
        # starttime=time.time()
        while(ABV > E1):
            if icon > 1000:
                GV = GF
                break
            else:
                icon += 1
            # Bubble size-Zhu & Zhang, 2015
            # DB = 6.034 * GV * (ST / DENL)**0.6 * (HP * Q / (DENL * ZI * VOI))**(-0.4) * (DENL / DENG)**0.2       # original
            # DB=1.4*(GV)**0.25*(ST/DENL)**0.6*(HP*Q/(DENL*ZI*VOI))**(-0.4)*(DENL/DENG)**0.2     # Bubble size-Zhu & Zhang, 2015


            # Original bubble size + original and Barrios combine CD is best for LambdaC2 prediction of GC6100



            # if(DB <= 0.0003) DB = 0.0003 
            # DBMAX = DB/0.43
            DB, DBMAX = self.get_DB(Q, HP,N, GF)
            REB = DENL*VSR*DB/VISL
            SR  = DB*O/VSR
                
            # # CD-Legendre & Magnaudet (1998), Clift et al. (1978), Rastello et al. (2011)
            # if(REB <= 50.0):
            #     CD  = 24.0/REB*(1.0+0.15*REB**0.687)*(1.0+0.3*SR**2.5)
            # else:
            #     CD  = 24.0/REB*(1.0+0.15*REB**0.687)*(1.0+0.55*SR**2.0)


        # Original bubble size + Barrios combine CD is best for LambdaC2 prediction of GC6100


            # # Barrios
            # YY = 0.00983 + 389.9 * REB / N ** 2
            # CD = (24.0 + 5.48 * (REB * YY) ** 0.427 + 0.36 * REB * YY) / (REB * YY)


            CD = self.CD_Cal(VSR, DENL, VISL, N, DB)



          
            # # Combine Legendre and Barrios
            # if(REB <= 40.0):
            #     CD  = 24.0/REB*(1.0+0.15*REB**0.687)*(1.0+0.3*SR**2.5)
            # elif(REB> 40.0 and REB <= 60.0):
            #     # transision by HWZ
            #     CD1  = 24.0/REB*(1.0+0.15*REB**0.687)*(1.0+0.3*SR**2.5)
            #     YY = 0.00983 + 389.9 * REB / N ** 2
            #     CD2 = (24.0 + 5.48 * (REB * YY) ** 0.427 + 0.36 * REB * YY) / (REB * YY)
            #     CD = CD1 * (60-REB)/20 + CD2 * (REB-40)/20
            # else:
            #     YY = 0.00983 + 389.9 * REB / N ** 2
            #     CD = (24.0 + 5.48 * (REB * YY) ** 0.427 + 0.36 * REB * YY) / (REB * YY)

            
            # # From drilling fluid, and Zhu's ABM report 2020 April
            # if GV < GF:
            #     GV = GV
            # VL = Q/ self.ZI/AI/GV
            # CD = self.C_Drag_Sphere(EDI, DI, VISL, DENL, VL, VSR, DBMAX)


            VSRN    = np.sqrt(4.0*DB*(DENL-DENG)*RI/(3.0*CD*DENL))*O
            ABV     = np.abs((VSRN-VSR)/VSR)
            VSR     = (VSRN + 9.0*VSR)/10.0
            RS      = VSR*(2.0*np.pi*RI-ZI*TB)*YI/(Q+QLK)

            if(GF < 0.0): 
                GV  = 0.0
            else:
                GV  = (RS-1.0+np.sqrt((1.0-RS)**2+4.0*RS*GF))/(2.0*RS)
        
        alphaG_crit = 0.5-0.25*np.exp(-(N/3500.0)**4.0)   #Fotran code selection
        if GV > alphaG_crit:
            #  GV= alphaG_crit
            # GV, _ = self.gv_sun_INT(self.B1, self.B2, DENL, DENG, N, Q, Q*GF/(1-GF), QLK, self.R1, self.R2, self.TB,
            #                         VISL, VISG, self.YI1, self.YI2, self.YI, self.ZI)

            GV, _, _, _ = self.ITFLOW_HWZ(AI, self.B1, self.B2, DI, EDI, DENL, DENG, self.LI, N, self.RI, VISL, VISG, ST, 
                        self.ZI, Q,  Q*GF/(1-GF), QLK, self.YI, self.TB, HP)
        # else:
            # FGL = 'BUB'
        # endtime=time.time()
        # print('BUFLOW time: ', starttime-endtime, 'icon: ', icon, 'Residual: ', ABV, '\n')
        return GV, CD, DB, VSR

    def ITFLOW(self,AI,B1,B2,DI,EDI,DENL,DENG,LI,N,RI,VISL,VISG,ST,ZI,Q,QG,QLK):
        # real*8,   intent(in )     ::  AI               ! representative cross sectional area of impeller channel (m2)
        # real*8,   intent(in )     ::  B1               ! blade angle from tangential at impeller inlet (degree)
        # real*8,   intent(in )     ::  B2               ! blade angle from tangential at impeller outlet (degree)	
        # real*8,   intent(in )     ::  DI               ! representative impeller diameter (m)
        # real*8,   intent(in )     ::  EDI              ! relative impeller wall roughness 
        # real*8,   intent(in )     ::  Q                ! liquid flow rate (m3/s)
        # real*8,   intent(in )     ::  QG               ! gas flow rate (m3/s)   
        # real*8,   intent(in )     ::  QLK              ! leakage flow rate (m3/s)
        # real*8,   intent(in )     ::  DENL             ! liquid density (kg/m3)
        # real*8,   intent(in )     ::  DENG             ! gas density (kg/m3)
        # real*8,   intent(in )     ::  G                ! gravitational acceleration (m/s2)
        # real*8,   intent(in )     ::  O                ! impeller angular velocity (radians/s)
        # real*8,   intent(in )     ::  LI               ! impeller channel length (m)	
        # real*8,   intent(in )     ::  VISL             ! liquid viscosity (kg/m-s)	
        # real*8,   intent(in )     ::  VISG             ! gas viscosity (kg/m-s)	
        # real*8,   intent(in )     ::  RI               ! average impeller radius (m)	
        # real*8,   intent(in )     ::  ST               ! surface tension (N/m)
        # real*8,   intent(in )     ::  ZI               ! impeller blade number
        # real*8,   intent(out)     ::  HL               ! sulg unit liquid holdup
        # real*8,   intent(out)     ::  GV
        # real*8 :: ANG = -90       ! Inclination angle, assume downward flow (deg)
        # real*8 :: FIC = 0.0142
        # real*8 VSL                ! superfacial liquid velocity or relative velocity
        # real*8 VSG                ! surperfacial gas velocity on INT-ANN boundary
        # real*8 AB,ABV 
        # real*8 ABCD, ABCDN
        # real*8 ABU 
        # real*8 AC                 ! cross sectional area of gas core (m2)             
        # real*8 AD  
        # real*8 AF                 ! cross sectional area of film (m2)
        # real*8 BI                 ! representative blade angle in impeller(deg)
        # real*8 CE         
        # real*8 CS                 ! slug length (m)
        # real*8 CF,CFN             ! slug film/Taylor bubble length (m)
        # real*8 CU                 ! slug unit length (m)
        # real*8 DC
        # real*8 DeltaL             ! Taylor bubble film thickness (m) 
        # real*8 DENC
        # real*8 DF 
        # real*8 DPAL,DPEX,DPSL
        # real*8 FF,FFN             ! friction factor between film and impeller wall
        # real*8 FI,FIM,FRA,FRA1     ! interfacial friction factor between gas and film
        # real*8 FIN
        # real*8 FM,FS,FSL 
        # real*8 FRO                ! Froude number
        # real*8 GC                 ! centrifugal accelaration (in radial direction)            
        # real*8 HLF                ! liquid holdup in film                           
        # real*8 HLFN   
        # real*8 HLS                ! liquid holdup in slug body                              
        # real*8 HLSN   
        # integer ICON   
        # real*8 REC1,REF,REG
        # real*8 REMX,REM
        # real*8 RESG   
        # real*8 RESL
        # real*8 RSU                ! ratio of slug length to the unit length   
        # real*8 RTH
        # real*8 SF                 ! perimeter wetted by film (m)
        # real*8 SHI
        # real*8 SI                 ! perimeter of interface (m)
        # real*8 THF,THFO           ! film thickness (m)                                                                   
        # real*8 VAV    
        # real*8 VC,VCN
        # real*8 V24    
        # real*8 VF                 ! film velocity (m/s)                             
        # real*8 VFN,VFN1 
        # real*8 VM                 ! mixture or slug velocity (m/s)                                                   
        # real*8 VT                 ! slug traslational (tail and front) velocity (m/s)            
        # real*8 VMT
        # real*8 WEB,WED            ! Weber number
        ANG = -90
        FIC = 0.0142        
        BI=(B1+B2)/2.0
        O = 2.0 * np.pi * N / 60.0
        GC=O**2*RI*np.sin(BI)
        VSL=(Q+QLK)/ZI/AI
        VSG=QG/ZI/AI
        VM=VSL+VSG
        HLS=1.0/(1.0+np.abs(VM/8.66)**1.39)
        CF = 0.5
        
        if(HLS < 0.24): 
            HLS=0.24
        
        # Translational velocity according to Nicklin (1962), Bendiksen (1984) and Zhang et al. (2000)
        REMX=np.abs(DENL*DI*VM/VISL)
        if(REMX < 2000.0):
            VAV=2.0*VM
        elif(REMX > 4000.0):
            VAV=1.2*VM
        elif(REMX > 2000.0 and REMX < 4000.0):
            VAV=(2.0-0.8*(REMX-2000.0)/2000.0)*VM
        VT=VAV+(0.54*np.cos(ANG)+0.35*np.sin(ANG))*np.sqrt(G*DI*np.abs(DENL-DENG)/DENL)
        if(VT < 0):
            VT=0.1
        VMT=(VM+VT)/2.0
        
        # Slug length      
        CS = (32.0*np.cos(ANG)**2+16.0*np.sin(ANG)**2)*DI
        CF = CS
        CU = CF+CS
        # CS = 30*DI
        CE=1.25-0.5*np.abs(np.sin(ANG))
        HLF= 0.001
        VF=VSL/2.0
        VFN1=VSL/2.0
        VC=VM
        ABCD=0.0
        FF=0.1
        FI=FIC
        # starttime=time.time()
        for i in range (1, 1001):
            # Overall liquid holdup
            HL=(HLS*(VT-VM)+VSL)/VT
            # HL=(HLF*CF+HLS*CS)/CU
            if (HL>1.0): 
                HL = 1.0/HL
            if (HL<0.0): 
                HL = -HL
            HLF=HLS*(VT-VM)/(VT-VF)
            if(HLF < 0.0): 
                HLF=-HLF
            if(HLF > HL): 
                HLF=0.99*HL              
            VCN=(VM-HLF*VF)/(1.0-HLF) 
            if(VCN < 0.0): 
                VCN=-VCN
            VC=VCN
            CFN=(HLS-HL)*CS/(HL-HLF)
            if (CFN<0.0): 
                CFN=-CFN
            ABU=np.abs((CFN-CF)/CF)
            CF=(CFN+9.0*CF)/10.0
            CU = CS + CF
            RSU=CS/CU
            # Taylor bubble geometries         
            DeltaL = DI/2.0*(1-np.sqrt(1.0-HLF))
            AC = np.pi*(DI-2.0*DeltaL)**2/4.0
            AF = np.pi*DeltaL*(DI-1.0*DeltaL)
            SI = np.pi*(DI-2.0*DeltaL)
            SF = np.pi*DI
            DF = 4.0*DeltaL*(DI-DeltaL)/DI
            DC = DI-2.0*DeltaL
            THF= DeltaL
            # Slug liquid holdup   
            DPEX=(DENL*(VM-VF)*(VT-VF)*HLF+DENG*(VM-VC)*(VT-VC)*(1.0-HLF))*DI/CS/4.0
            REM=np.abs(DENL*VM*DI/VISL)
        
            FM=self.get_fff(REM,EDI)
        
            DPSL=FM*DENL*VM*VM/2.0
            DPAL=DPSL+DPEX
            if(REM < 5000.0):
                DPAL=DPAL*(REM/5000.0)
            AD=DPAL/(3.16*CE*np.sqrt(ST*np.abs(DENL-DENG)*GC))
            HLSN=1.0/(1.0+AD)
            if(HLSN < 0.24):
                HLSN=0.24
            HLS=HLSN
        
            if(REM < 1500.0 and HLS < 0.24):
                HLS=0.24
        
            # if(VSL/VM > HLS or np.abs(CF) < DI) return
            # Reynolds numbers 
            REF=np.abs(DENL*VF*DF/VISL) 
            REC1=np.abs(DENG*VC*DC/VISG)
            # Friction factors
            FFN=self.get_fff(REF,EDI)
            FIM=self.get_fff(REC1,THF/DI)
            FF=(FFN+9.0*FF)/10.0
            # Interfacial friction factor (annular) according to Ambrosini et al. (1991)
            REG=np.abs(VC*DENG*DI/VISG)
            WED=DENG*VC*VC*DI/ST
            FS=0.046/REG**0.2
            SHI=np.abs(FI*DENG*(VC-VF)**2/2.0) 
            THFO=THF*np.sqrt(np.abs(SHI*DENG))/VISG
            FRA=FS*(1.0+13.8*(THFO-200.0*np.sqrt(DENG/DENL))*WED**0.2/REG**0.6)
            FRA1=self.get_fff(REC1,THF/DI) 
            if(FRA > FRA1):
                FRA=FRA1
            FIN = FRA
            # if(FIN > 0.5) FIN=0.5 
            FI=(FIN+9.0*FI)/10.0      
            # Calculate film length CF using the combined momentum eqaution    
            FSL=(DENL*(VM-VF)*(VT-VF)-DENG*(VM-VC)*(VT-VC))/CF
            ABCDN=(FSL-SI*FI*DENG*(VC-VF)**2/2.0*(1.0/AF+1.0/AC)+(DENL-DENG)*GC)*2.0*AF/(SF*FF*DENL)
            ABCD=(ABCDN+19.0*ABCD)/20.0
            if(ABCD > 0.0): 
                VFN=np.sqrt(ABCD) 
                # if(VFN > VM):
                    # VFN=VM
            else: 
                VFN=np.sqrt(-ABCD) 
                # if(VFN < -VM):
                   # VFN=-VM
            ABV=np.abs((VFN-VF)/VF)
            VF=(VFN+9.0*VF)/10.0  
        
            if(ABU < 0.0001 and ABV < 0.0001):
                break

        if (LI<=CU):
            GV=1.0-HL*(LI/CU)
        else:
            GV=1.0-HL

        # endtime=time.time()
        # print('ITFLOW time: ', starttime-endtime, 'icon: ', i, 'Residual: ', ABU, '\n')
        return GV, HL

    def ITFLOW_HWZ(self,AI,B1,B2,DI,EDI,DENL,DENG,LI,N,RI,VISL,VISG,ST,ZI,Q,QG,QLK,YI,TB,HP):
        # Based on Sun drag coefficient
        VSR = 0.1
        ABV = 1.0
        GV  = 0.5
        REB = 0.0
        RS = 0.0
        SR = 0.0
        VSRN = 0.0
        O = 2.0 * np.pi * N / 60.0
        icon=0

        AIW = self.AB + self.ASF + self.ASB
        DI = 4.0 * self.VOI / AIW
        EDI = self.EA / DI
        AI = self.VOI / self.LI
        # starttime=time.time()

        GF = QG / (Q + QLK + QG)
        CD_to_rb = 0.1
        counter = 0
        OMEGA = 2.0 * pi * N / 60.0

        C1M_L = (Q + QLK) / ((2.0 * pi * RI - ZI * TB) * YI)
        W1_L = C1M_L / np.sin((B1+B2/2))
        C1M_G = QG / ((2.0 * pi * RI - ZI * TB) * YI)
        W1_G = C1M_G / np.sin((B1+B2/2))


        while(ABV > E1):
            if icon > 1000:
                GV = GF
                break
            else:
                icon += 1
            
            DB, DBMAX = self.get_DB(Q, HP,N, GV)
            REB = DENL*VSR*DB/VISL

            W1L = W1_L / (1 - GV)
            W1G = W1_G / GV
            mium1 = (1 - GV) * VISL + GV * VISG
            DENM1 = (1 - GV) * DENL + GV * DENG
            CD_to_rb1 = 12. * mium1 * (CD_INT_EFF*9.13e7 * GV ** CD_GV_EFF / W1_G ** CD_Gas_EFF * W1_L ** CD_Liquid_EFF ) / (np.abs(W1G - W1L) * DENM1) * (N/3600)**2

            SR  = DB*O/VSR
                
            # # CD  = 24.0/REB*(1.0+0.15*REB**0.687)*(1.0+0.55*SR**2.0)
            # CD  = 24.0/REB*(1.0+CD_Liquid_EFF*REB**0.687)*(1.0+CD_INT_EFF*SR**2.0)
            # CD_to_rb1 = CD/DB

            VSRN    = np.sqrt(4.0 * (DENL - DENG) * RI / (3.0 * CD_to_rb1 * DENL)) * OMEGA
            ABV     = np.abs((VSRN-VSR)/VSR)
            VSR     = (VSRN + 9.0*VSR)/10.0
            RS      = VSR*(2.0*np.pi*RI-ZI*TB)*YI/(Q+QLK)
            

            if(GF < 0.0): 
                GV  = 0.0
            else:
                GV  = (RS-1.0+np.sqrt((1.0-RS)**2+4.0*RS*GF))/(2.0*RS)
        
        return GV, CD_to_rb1, CD_to_rb1, VSR

    def ANFLOW(self,AI,B1,B2,DI,EDI,DENL,DENG,N,RI,VISL,VISG,ST,ZI,Q,QG,QLK):

        # real*8,   intent(in )     ::  AI            ! representative cross sectional area of impeller channel (m2)
        # real*8,   intent(in )     ::  B1            ! blade angle from tangential at impeller inlet (degree)
        # real*8,   intent(in )     ::  B2            ! blade angle from tangential at impeller outlet (degree)   
        # real*8,   intent(in )     ::  DI            ! representative impeller diameter (m)
        # real*8,   intent(in )     ::  EDI           ! relative impeller wall roughness 
        # real*8,   intent(in )     ::  Q             ! liquid flow rate (m3/s)
        # real*8,   intent(in )     ::  QG            ! gas flow rate (m3/s)   
        # real*8,   intent(in )     ::  QLK           ! leakage flow rate (m3/s)
        # real*8,   intent(in )     ::  DENL          ! liquid density (kg/m3)
        # real*8,   intent(in )     ::  DENG          ! gas density (kg/m3)
        # real*8,   intent(in )     ::  G             ! gravitational acceleration (m/s2)
        # real*8,   intent(in )     ::  O             ! impeller angular velocity (radians/s) 
        # real*8,   intent(in )     ::  VISL          ! liquid viscosity (kg/m-s) 
        # real*8,   intent(in )     ::  VISG          ! gas viscosity (kg/m-s)    
        # real*8,   intent(in )     ::  PI            ! 3.141592657
        # real*8,   intent(in )     ::  RI            ! average impeller radius (m)   
        # real*8,   intent(in )     ::  ST            ! surface tension (N/m)
        # real*8,   intent(in )     ::  ZI            ! impeller blade number
        # character*(*), intent(out)::  FGL            ! flow pattern
        # real*8,   intent(out)     ::  GV

        # real*8 :: ANG = -90       ! Inclination angle, assume downward flow (deg)
        # real*8 :: FIC = 0.0142
        # real*8 VSL                ! superfacial liquid velocity or relative velocity
        # real*8 VSG                ! surperfacial gas velocity on INT-ANN boundary
        # real*8 AB 
        # real*8 ALPHAC
        # real*8 ABCD
        # real*8 AC                 ! cross sectional area of gas core (m2)             
        # real*8 AD  
        # real*8 AF                 ! cross sectional area of film (m2)
        # real*8 BI                 ! representative blade angle in impeller(deg)
        # real*8 CCC         
        # real*8 CS                 ! slug length (m)
        # real*8 CF                 ! slug film/Taylor bubble length (m)
        # real*8 CU                 ! slug unit length (m)
        # real*8 DC
        # real*8 DeltaL             ! Taylor bubble film thickness (m) 
        # real*8 DENC
        # real*8 DF
        # real*8 FE,FEN 
        # real*8 FF                 ! friction factor between film and impeller wall
        # real*8 FI,FIM,FRA,FRA1     ! interfacial friction factor between gas and film
        # real*8 FIN
        # real*8 FM,FS 
        # real*8 FRO                ! Froude number
        # real*8 GC                 ! centrifugal accelaration (in radial direction)            
        # real*8 HLF                ! liquid holdup in film                           
        # integer ICON   
        # real*8 REC1,REF,REG
        # real*8 REMX
        # real*8 RESG   
        # real*8 RESL
        # real*8 SF                 ! perimeter wetted by film (m)
        # real*8 SHI
        # real*8 SI                 ! perimeter of interface (m)
        # real*8 THF,THFO           ! film thickness (m)                                                                       
        # real*8 VC,VCN   
        # real*8 VF,VFN             ! film velocity (m/s)                              
        # real*8 VM                 ! mixture or slug velocity (m/s)                                                   
        # real*8 VISC               ! gas core viscosity (kg/m-s) 
        # real*8 WEB,WED            ! Weber number
        ANG = -90
        FIC = 0.0142
        FI=FIC
        BI=(B1+B2)/2.0
        O = 2.0 * np.pi * N / 60.0
        GC=O**2*RI*np.sin(BI)
        VSL=(Q+QLK)/ZI/AI
        VSG=QG/ZI/AI
        VM=VSL+VSG
        VF=VSL/2.0
        VC=VM
        # Guess a HLF
        HLF=0.999
        # starttime=time.time()
        for i in range(1, 1001):
        # Update HLF
            HLF=HLF-0.001           
        # Entrainement fraction basd on Oliemans et al's (1986)
            WEB=DENG*VSG*VSG*DI/ST
            FRO=np.sqrt(G*DI)/VSG
            RESG=DENG*VSG*DI/VISG
            RESL=DENL*VSL*DI/VISL
            CCC=0.003*WEB**1.8*FRO**0.92*RESL**0.7*(DENL/DENG)**0.38*(VISL/VISG)**0.97/RESG**1.24
            FEN=CCC/(1.0+CCC)
            if(FEN > 0.9):
                FE=0.9
            FE=FEN
        # Film geometries
            DeltaL = DI/2.0*(1-np.sqrt(1.0-HLF))
            AC = np.pi*(DI-2.0*DeltaL)**2/4.0
            AF = np.pi*DeltaL*(DI-1.0*DeltaL)
            SI = np.pi*(DI-2.0*DeltaL)
            SF = np.pi*DI
            DF = 4.0*DeltaL*(DI-DeltaL)/DI
            DC = DI-2.0*DeltaL
            THF= DeltaL
            
            VFN= VSL*(1.0-FE)/HLF
            VF = (VFN+9.0*VF)/10.0
            VCN=(VSG+VSL*FE)/(1-HLF)
            VC = (VCN+9.0*VC)/10.0
                    
            ALPHAC = VSG/(VSG+VSL*FE)
            DENC=DENG*ALPHAC+DENL*(1.0-ALPHAC)
            VISC=VISG*ALPHAC+VISL*(1.0-ALPHAC)
            REF=np.abs(DENL*VF*DF/VISL) 
            REC1=np.abs(DENG*VC*DC/VISC)

        # Friction factors
            FF = self.get_fff(REF,EDI)
            FIM = self.get_fff(REC1,THF/DI)
            
        # Interfacial friction factor (annular) according to Ambrosini et al. (1991)
            REG=np.abs(VC*DENG*DI/VISG)
            WED=DENG*VC*VC*DI/ST
            FS=0.046/REG**0.2
            SHI=np.abs(FI*DENG*(VC-VF)**2/2.0) 
            THFO=THF*np.sqrt(np.abs(SHI*DENG))/VISG
            FRA=FS*(1.0+13.8*(THFO-200.0*np.sqrt(DENG/DENL))*WED**0.2/REG**0.6)
            FRA1 = self.get_fff(REC1,THF/DI) 
            if(FRA > FRA1):
                FRA=FRA1
            FIN=FRA
            # if(FIN > 0.5): FIN=0.5 
            FI=(FIN+9.0*FI)/10.0 
            
            ABCD=(DENL-DENC)*GC-SF*FF*DENL*VF*VF/(2.0*AF)+SI*FI*DENC*(VC-VF)*np.abs(VC-VF)*(1.0/2.0/AF+1.0/2.0/AC)
        
            if (ABCD < 0.001):
                break
            if (HLF < 0.0001):
                FGL='GAS'
                break
        # write(*,*) "ANFLOW",ICON, HLF, DI

        GV=ALPHAC*(1.0-2.0*DeltaL/DI)**2.0

        # endtime=time.time()
        # print('ANFLOW time: ', starttime-endtime, 'icon: ', i, 'Residual: ', ABCD, 'HLF: ', HLF,'\n')
        return GV
  
    def gv_zhu_flowpattern(self, GF, DENL, DENG, VISL, VISG, HP, N, QL, QG, QLK, ST, AI, DI, EDI, DENW):
        GVC1 = self.get_lambda_c1(QL+QLK,HP,N)
        GVC2 = self.get_lambda_c2(QL+QLK,HP,N)
        GVC3 = self.get_lambda_c3(QL+QLK,N)
        if GVC1 >= GVC2:
            GVC2 = GVC1
        if GVC2 >= GVC3:
            GVC3 = GVC2

        CF = 0
        if GF < GVC1:
            FGL = 'D-B'
            # GV = self.DBFLOW(GF)
            
            GV, _, _, _ = self.BUFLOW(DENL, DENG, GF, HP, N, QL, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
            GV1, _, _, _ = self.ITFLOW_HWZ(AI, self.B1, self.B2, DI, EDI, DENL, DENG, self.LI, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK, self.YI, self.TB,HP)
            GV2, _, _, _ = self.BUFLOW(DENL, DENG, GF, HP, N, QL, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
            if GV2<GV1:
                GV=GV2
            else:
                GVC2=GVC2*0.9   # safety factor, somnetime Pbubble < Pslug
                QLC = QG / (GVC2/(1-GVC2))
                Qrange = transition_zone*self.QBEM
                if QL > QLC-Qrange:
                    GV2, _, _, _ = self.BUFLOW(DENL, DENG, GF, HP, N, QLC, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
                    if GV2 < GV1:
                        GV = (QLC-QL)/Qrange*GV1+(1-(QLC-QL)/Qrange)*GV2
                    else:
                        GV = GV1
        elif GF < GVC2:
            FGL = 'BUB'
            # GV, _, _, FGL = self.BUFLOW(DENL, DENG, GF, HP, N, QL, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
            GV1, _, _, _ = self.ITFLOW_HWZ(AI, self.B1, self.B2, DI, EDI, DENL, DENG, self.LI, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK, self.YI, self.TB,HP)
            GV2, _, _, _ = self.BUFLOW(DENL, DENG, GF, HP, N, QL, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
            if GV2<GV1:
                GV=GV2
            else:
                GVC2=GVC2*0.9   # safety factor, somnetime Pbubble < Pslug
                QLC = QG / (GVC2/(1-GVC2))
                Qrange = transition_zone*self.QBEM
                if QL > QLC-Qrange:
                    GV2, _, _, _ = self.BUFLOW(DENL, DENG, GF, HP, N, QLC, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
                    if GV2 < GV1:
                        GV = (QLC-QL)/Qrange*GV1+(1-(QLC-QL)/Qrange)*GV2
                    else:
                        GV = GV1
        elif GF < GVC3:
            FGL = 'INT'
            # GV, _ = self.ITFLOW(AI, self.B1, self.B2, DI, EDI, DENL, DENG, self.LI, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK)
            # GV, _ = self.gv_sun_INT(self.B1, self.B2, DENL, DENG, N, QL, QG, QLK, self.R1, self.R2, self.TB,
            #                         VISL, VISG, self.YI1, self.YI2, self.YI, self.ZI)
            
            # GV, _, _, _ = self.ITFLOW_HWZ(AI, self.B1, self.B2, DI, EDI, DENL, DENG, self.LI, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK, self.YI, self.TB,HP)
            GVC2=GVC2*0.9   # safety factor, somnetime Pbubble < Pslug
            QLC = QG / (GVC2/(1-GVC2))
            Qrange = transition_zone*self.QBEM
            if QL > QLC-Qrange:
                GV1, _, _, _ = self.ITFLOW_HWZ(AI, self.B1, self.B2, DI, EDI, DENL, DENG, self.LI, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK, self.YI, self.TB,HP)
                GV2, _, _, _ = self.BUFLOW(DENL, DENG, GF, HP, N, QLC, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
                if GV2 < GV1:
                    GV = (QLC-QL)/Qrange*GV1+(1-(QLC-QL)/Qrange)*GV2
                else:
                    GV = GV1
            else:
                GV, _, _, _ = self.ITFLOW_HWZ(AI, self.B1, self.B2, DI, EDI, DENL, DENG, self.LI, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK, self.YI, self.TB,HP)
            
        else:
            FGL = 'ANN'
            # GV = self.ANFLOW(AI, self.B1, self.B2, DI, EDI, DENL, DENG, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK)
            GV, _, _, _ = self.ITFLOW_HWZ(AI, self.B1, self.B2, DI, EDI, DENL, DENG, self.LI, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK, self.YI, self.TB,HP)
            
        return GV, FGL, CF

    def gv_zhu_flowpattern_JJZ(self, GF, DENL, DENG, VISL, VISG, HP, N, QL, QG, QLK, ST, AI, DI, EDI, DENW):
        GVC1 = self.get_lambda_c1(QL+QLK,HP,N)
        GVC2 = self.get_lambda_c2(QL+QLK,HP,N)
        GVC3 = self.get_lambda_c3(QL+QLK,N)
        if GVC1 >= GVC2:
            GVC2 = GVC1
        if GVC2 >= GVC3:
            GVC3 = GVC2

        CF = 0
        if GF < GVC1:
            FGL = 'D-B'
            # GV = self.DBFLOW(GF)
            
            GV, _, _, _ = self.BUFLOW(DENL, DENG, GF, HP, N, QL, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
            GV1, _, _, _ = self.ITFLOW_HWZ(AI, self.B1, self.B2, DI, EDI, DENL, DENG, self.LI, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK, self.YI, self.TB,HP)
            GV2, _, _, _ = self.BUFLOW(DENL, DENG, GF, HP, N, QL, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
            if GV2<GV1:
                GV=GV2
            else:
                GVC2=GVC2*0.9   # safety factor, somnetime Pbubble < Pslug
                QLC = QG / (GVC2/(1-GVC2))
                Qrange = transition_zone*self.QBEM
                if QL > QLC-Qrange:
                    GV2, _, _, _ = self.BUFLOW(DENL, DENG, GF, HP, N, QLC, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
                    if GV2 < GV1:
                        GV = (QLC-QL)/Qrange*GV1+(1-(QLC-QL)/Qrange)*GV2
                    else:
                        GV = GV1
        elif GF < GVC2:
            FGL = 'BUB'
            # GV, _, _, FGL = self.BUFLOW(DENL, DENG, GF, HP, N, QL, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
            GV1, _, _, _ = self.ITFLOW_HWZ(AI, self.B1, self.B2, DI, EDI, DENL, DENG, self.LI, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK, self.YI, self.TB,HP)
            GV2, _, _, _ = self.BUFLOW(DENL, DENG, GF, HP, N, QL, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
            if GV2<GV1:
                GV=GV2
            else:
                GVC2=GVC2*0.9   # safety factor, somnetime Pbubble < Pslug
                QLC = QG / (GVC2/(1-GVC2))
                Qrange = transition_zone*self.QBEM
                if QL > QLC-Qrange:
                    GV2, _, _, _ = self.BUFLOW(DENL, DENG, GF, HP, N, QLC, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
                    if GV2 < GV1:
                        GV = (QLC-QL)/Qrange*GV1+(1-(QLC-QL)/Qrange)*GV2
                    else:
                        GV = GV1
        elif GF < GVC3:
            FGL = 'INT'
            GV, _ = self.ITFLOW(AI, self.B1, self.B2, DI, EDI, DENL, DENG, self.LI, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK)
            
        else:
            FGL = 'ANN'
            GV = self.ANFLOW(AI, self.B1, self.B2, DI, EDI, DENL, DENG, N, self.RI, VISL, VISG, ST, self.ZI, QL, QG, QLK)
            
        return GV, FGL, CF

    # def gl_calculate_old(self, Q, QG, flg='Z'):
    #     """
    #     Calcualte air-water flow performance of ESP
    #     :param QG: gas flow rate (bpd)
    #     :param Q: liquid flow rate (bpd)
    #     :param flg: 'Z': Zhu model; 'S': Sun model; 'B': Barrios model
    #     :return: PP, PE, PF, PT, PD, PRE, PLK, QLK, GV
    #     """
    #     # run single-phase calculation to initialize, options includes: old, new, 2018, jc... 
    #     HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_new(ql, self.QBEM, self.DENL, DENW, self.N, self.NS,
    #                                                                     self.SGM, self.SN, self.ST, self.VISL, VISW,
    #                                                                     self.WC)

    #     # convert field units
    #     PP = HP * psi_to_pa
    #     Q = Q * bbl_to_m3 / 24.0 / 3600.0
    #     QG = QG * bbl_to_m3 / 24.0 / 3600.0
    #     GF = QG / (Q + QG)

    #     icon = 0
    #     ABP = 1.
    #     PE, PEE = 0., 0.
    #     PFI, PFD = 0., 0.
    #     PTI, PTD = 0., 0.
    #     GV = 0.
    #     PLK = 0.

    #     AIW = self.AB + self.ASF + self.ASB
    #     ADW = self.AV + self.ADF + self.ADB
    #     AI = self.VOI / self.LI
    #     AD = self.VOD / self.LD
    #     DI = 4.0 * self.VOI / AIW
    #     DD = 4.0 * self.VOD / ADW
    #     EDI = self.EA / DI
    #     EDD = self.EA / DD
    #     QBEM=self.QBEM * (self.N / 3600) * (self.VISL / self.VISW) ** (0.01 * (3448 / self.NS) ** 4)
    #     DEND = self.DENL * (1.0 - GF) + self.DENG * GF

    #     while ABP > E1:
    #         VI = (Q + QLK) / self.ZI / AI
    #         VD = Q / self.ZD / AD
    #         C1M = (Q + QLK) / ((2.0 * pi * self.R1 - self.ZI * self.TB) * self.YI1)
    #         C2M = (Q + QLK) / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
    #         U1 = self.R1 * self.OMEGA
    #         U2 = self.R2 * self.OMEGA
    #         W1 = C1M / np.sin(self.B1)
    #         W2 = C2M / np.sin(self.B2)
    #         C1 = np.sqrt(C1M ** 2 + (U1 - C1M / np.tan(self.B1)) ** 2)
    #         C2 = np.sqrt(C2M ** 2 + (U2 - C2M / np.tan(self.B2)) ** 2)
    #         CMB = QBEM / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)
    #         C2B = np.sqrt(CMB ** 2 + (U2 - CMB / np.tan(self.B2)) ** 2)

    #         HE = (U2 ** 2 - U2 * C2M / np.tan(self.B2)) / G

    #         if GF <= E1:
    #             PP = HP * psi_to_pa
    #             break
    #         elif GF >= (1.-E1):
    #             PP = HP * psi_to_pa * self.DENG / DENW
    #             break
            
    #         #Flow pattern transision and GV calculation
    #         if flg == 'Z':
    #             GV = self.bubble_flow_zhu(HP * psi_to_pa, GF, Q, QLK)
    #         elif flg == 'S':
    #             GV = self.bubble_flow_sun(GF, Q, QLK)
    #         elif flg == 'B':
    #             GV = self.bubble_flow_Barrios(GF, Q, QLK)
    #         else:
    #             GV = self.bubble_flow_zhu(HP * psi_to_pa, GF, Q, QLK)
    #         #End flow pattern and GV

    #         DENI = self.DENL * (1.0 - GV) + self.DENG * GV
    #         PE = HE * DENI * G
    #         VLR = (Q + QLK) / ((2.0 * pi * self.RI - self.ZI * self.TB) * self.YI * (1.0 - GV))
    #         VGR = QG / ((2.0 * pi * self.RI - self.ZI * self.TB) * self.YI * GV)

    #         if Q+QLK < QBEM:
    #             VSH = U2 * (QBEM - Q) / QBEM
    #             C2F = C2B * Q / QBEM
    #             DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
    #             REC1 = DENI * VSH * DC / self.VISL
    #             # SGM = (self.VISW / self.VISL) ** 0.5 / (1.0 + 0.02 * REC1 ** 0.2)
    #             SGM = 0.0
    #             C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
    #             C2E = C2F + SGM * (C2P - C2F)
    #             PEE = PE + DENI * (C2E ** 2 - C2 ** 2) / 2.0
    #         else:
    #             VSH = U2 * (Q - QBEM) / QBEM
    #             C2F = C2B * Q / QBEM
    #             C2E = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
    #             PEE = PE + DENI * (C2E ** 2 - C2 ** 2) / 2.0

    #         REI = DENI * (W1 + W2) * DI / self.VISL / 2.0
    #         RED = DEND * VD * DD / self.VISL
    #         FFI = self.get_fff(REI, EDI)
    #         FFD = self.get_fff(RED, EDD)
    #         PFI = 1.5 * 4.0 * FFI * DENI * VI**2 * self.LI / (2.0*DI)
    #         PFD = 1.5 * 4.0 * FFD * DEND * VD ** 2 * self.LD / (2.0 * DD)

    #     #Turning loss coefficient
    #         # if (self.DENL*4*self.R2**2*self.OMEGA**2/self.VISL) > 1e7:
    #         #     FTI = 6.5
    #         #     FTD = 2.5
    #         # else:
    #         #     FTI = 12.
    #         #     FTD = 3.
    #         FTI=3.0
    #         FTD=3.0
    #     #Turning loss coefficient

    #         PTI = FTI * DENI * VI ** 2 / 2.0
    #         PTD = FTD * DEND * VD ** 2 / 2.0
    #         PPN = PEE - PFI - PFD - PTI - PTD

    #         UL = self.RLK * self.OMEGA
    #         PIO = PEE - PFI - PTI
    #         PLK = PIO - DENI * (U2 ** 2 - UL ** 2) / 8.0

    #         if PLK >= 0:
    #             VL = np.abs(QLK) / (2.0 * pi * self.RLK * self.SL)
    #             REL = np.abs(DEND * VL * self.SL / self.VISL)
    #             EDL = 0.0
    #             FFL = self.get_fff(REL, EDL)
    #             VL = np.sqrt(2.0 * PLK / (1.5 + 4.0 * FFL * self.LG / self.SL) / DEND)
    #             QLKN = 2.0 * pi * self.RLK * self.SL * VL
    #         else:
    #             VL = np.abs(QLK / (2.0 * pi * self.RLK * self.SL))
    #             REL = self.DENL * VL * self.SL / self.VISL
    #             EDL = 0.0
    #             FFL = self.get_fff(REL, EDL)
    #             VL = np.sqrt(2.0 * np.abs(PLK) / (1.5 + 4.0 * FFL * self.LG / self.SL) / DEND)
    #             QLKN = -2.0 * pi * self.RLK * self.SL * VL

    #         QLK = QLKN
    #         ABP = np.abs((PPN-PP)/PPN)
    #         PP = PPN

    #         if icon > 500:
    #             break
    #         else:
    #             icon += 1

    #     # return pressure in psi, flow rate in bpd
    #     PP = PP / psi_to_pa
    #     PE = PE / psi_to_pa
    #     PEE = PEE / psi_to_pa
    #     PF = (PFI + PFD) / psi_to_pa
    #     PT = (PTI + PTD) / psi_to_pa
    #     PD = (PFD + PTD) / psi_to_pa
    #     PRE = np.abs(PE - PEE)
    #     PLK = PLK / psi_to_pa
    #     QLK = QLK * 24.0 * 3600.0 / bbl_to_m3
    #     return PP, PE, PF, PT, PD, PRE, PLK, QLK, GV

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
        FGL = 'Other'   # flow pattern
        # # a tricky adjusting of rotational speed in Sun's model
        # if (flg == 'S') and N < 3500:
        #     N += 600

        # run single-phase calculation to initialize
        if sgl_model=='zhang_2015':
            HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_old(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW,
                                                                   WC)         
        elif sgl_model=='zhang_2016':
            HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_new(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW,
                                                                   WC)   
        elif sgl_model=='jiecheng_2017':
            HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_jc(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW,
                                                                   WC)   
        elif sgl_model=='zhu_2018':
            HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_2018(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW,
                                                                   WC)   

        if HP < 0:
            HP = HP
        # convert filed units
        PP = HP * psi_to_pa
        # QG = 
        GF = QG / (QL + QG)

        icon = 0
        ABP = 1.
        PE, PEE = 0., 0
        PFI, PFD = 0., 0
        PTI, PTD = 0, 0
        GV = 0.
        PLK = 0.
        HLKloss = 0.
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
        VISD = VISL * (1.0 - GF) + VISG * GF
        CF = 0
        # new QBEM due to liquid viscosity

        QBEM = QBEM * (N / 3600) * (VISL / VISW) ** (0.01 * (3448 / NS) ** 4)   #original
        # QBEM = QBEM * (N / 3600) * (VISL / VISW) ** (0.01 * (3448 / NS) ** 4)*DENL/(DENL * (1.0 - GF) + DENG * GF)   #add density effect by Haiwen

        # check if emulsion occurs
        if WC > 0.:
            VISL = self.emulsion(self.VOI, self.R2, VISL, VISW, DENL, DENW, WC, ST, N, QL, SN)

        # while ABP > E1 or ABQ > E1:
        while ABP > E1:
            VI = (QL + QG + QLK) / self.ZI / AI     # add QG by Haiwen Zhu
            VD = (QL + QG) / self.ZD / AD       # add QG by Haiwen Zhu
            C1M = (QL + QG + QLK) / ((2.0 * pi * self.R1 - self.ZI * self.TB) * self.YI1)        # add QG by Haiwen Zhu
            C2M = (QL + QG + QLK) / ((2.0 * pi * self.R2 - self.ZI * self.TB) * self.YI2)        # add QG by Haiwen Zhu
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
                PP = HP * psi_to_pa * (DENG * GF + DENL * (1 - GF)) / DENL
                GV = GF
                FGL = 'D-B'
                break
            elif GF >= (1-E1):
                PP = HP * psi_to_pa * (DENG * GF + DENL * (1 - GF)) / DENL
                GV = GF
                FGL = 'ANN'
                break

            if flg == 'Z':
                # qld = (QL + QLK) / (2.91492 * N * bbl_to_m3 / 3600 / 24)        # for GC-6100 only
                # qgd = (DENG / DENL) ** 0.2 * (OMEGA * (2 * self.R1) ** 2 / (VISL / DENL)) ** 0.4 * \
                #       (0.102 * np.exp(qld)) ** 4.4682
                # lambdaC2 = qgd / (qgd + qld)
                # lambdaC3 = self.get_lambda_c3(QL + QLK)
                GV, _ = self.gv_zhu(DENL, DENG, GF, HP*psi_to_pa, N, QL, QLK, self.RI, self.ST, self.TB, VISL, self.VOI, self.YI, self.ZI)
                # GV, _ = self.gv_zhu(DENL, DENG, GF, HP*psi_to_pa, N, QL, QLK, self.RI, ST, self.TB, VISL, self.VOI, self.YI, self.ZI)
                # GV, _, _, FGL = self.BUFLOW(DENL, DENG, GF, HP*psi_to_pa, N, QL, self.RI, ST, VISL, VISG, self.VOI, self.YI, self.ZI, QLK, self.TB)
           
            elif flg == 'S':
                # GV, _ = self.gv_sun(self.B1, self.B2, DENL, DENG, N, QL, QG, QLK, self.R1, self.R2, self.TB,
                #                     VISL, VISG, self.YI1, self.YI2, self.YI, self.ZI)
                GV, FGL, CF = self.gv_zhu_flowpattern_JJZ(GF, DENL, DENG, VISL, VISG, HP*psi_to_pa, N, QL, QG, QLK, ST, AI, DI, EDI, DENW)
            elif flg == 'G':
                GV, _ = self.gv_sun_gamboa(self.B1, self.B2, DENL, DENG, N, QL, QG, QLK, self.R1, self.R2, ST, self.TB,
                                           VISL, VISG, self.YI1, self.YI2, self.YI, self.ZI)
            elif flg == 'B':
                GV, _= self.gv_Barrios(self.DENL, DENG, GF, N, QL, QLK, self.R1, self.RI, ST, self.TB, VISL, self.YI, self.ZI)
            elif flg == 'F':
                
                GV, FGL, CF = self.gv_zhu_flowpattern(GF, DENL, DENG, VISL, VISG, HP*psi_to_pa, N, QL, QG, QLK, ST, AI, DI, EDI, DENW)

                # GV, _ = self.gv_zhu(DENL, DENG, GF, HP*psi_to_pa, N, QL, QLK, self.RI, ST, self.TB, VISL, self.VOI, self.YI, self.ZI)



                # if PP < 0:
                #     GV, FGL = self.gv_zhu_flowpattern(GF, DENL, DENG, VISL, VISG, HP*psi_to_pa*(1-GV), N, QL, QG, QLK, ST, AI, DI, EDI)
                # else:
                #     GV, FGL = self.gv_zhu_flowpattern(GF, DENL, DENG, VISL, VISG, PP, N, QL, QG, QLK, ST, AI, DI, EDI)
                # if FGL == 'INT':
                #     GV2, _ = self.gv_Barrios(self.DENL, self.DENG, GF, N, QL, QLK, self.R1, self.RI, ST, self.TB, VISL, self.YI, self.ZI)
                #     GV = (0.5*GV+0.5*GV2)
            elif flg == 'H':
                GV = self.DBFLOW(GF)
            else:
                GV, _ = self.gv_sun(self.B1, self.B2, DENL, DENG, N, QL, QG, QLK, self.R1, self.R2, self.TB,
                                    VISL, VISG, self.YI1, self.YI2, self.YI, self.ZI)

            # DENI = DENL * (1.0 - GV) + DENG * GV
            # VISI = VISL * (1.0 - GV) + VISG * GV        # revised by Haiwen zhu, should consider gas effect to viscosity
            # PE = HE * DENI * G
            # VLR = (QL + QLK) / ((2.0 * pi * self.RI - self.ZI * self.TB) * self.YI * (1.0 - GV))
            # VGR = QG / ((2.0 * pi * self.RI - self.ZI * self.TB) * self.YI * GV)

            # if (QL + QG + QLK) <= QBEM:  # add QG by Haiwen Zhu
            #     VSH = U2 * (QBEM - (QL + QG + QLK)) / QBEM   # add QG by Haiwen Zhu
            #     C2F = C2B * (QL + QG + QLK) / QBEM
            #     DC = 2.0 * pi * self.R2 * np.sin(self.B2) / self.ZI
            #     REC1 = DENI * VSH * DC / VISI       # changed VISL to VISI by Haiwen Zhu
            #     C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
            #     C2E = C2F
            #     PEE = PE + DENI * (C2E ** 2 - C2 ** 2) / 2.0
            # else:   
            #     VSH = U2 * (QL + QG + QLK - QBEM) / QBEM # add QG by Haiwen Zhu
            #     C2F = C2B * (QL + QG + QLK) / QBEM   # add QG by Haiwen Zhu
            #     C2P = (C2 ** 2 + C2F ** 2 - VSH ** 2) / (2.0 * C2F)
            #     C2E = C2F + SGM * (C2P - C2F) * (QL + QG + QLK - QBEM) / QBEM    # add QG by Haiwen Zhu
            #     PEE = PE + DENI * (C2E ** 2 - C2 ** 2) / 2.0

            # # REI = DENI * (W1 + W2) * DI / VISL / 2.0 / 2.0        #not sure why divided by 2
            # # RED = DEND * VD * DD / VISL / 2.0         # not sure why divided by 2
            # REI = DENI * (W1 + W2) * DI / VISI / 2.0        #changed VISL to VISI, removed /2 by Haiwen
            # RED = DEND * VD * DD / VISD       # changed VISD to VISD, removed /2 by Haiwen Zhu
            # FFI = self.get_fff(REI, EDI)
            # FFD = self.get_fff(RED, EDD)
            # PFI = 2.5 * 4.0 * FFI * DENI * VI ** 2 * self.LI / (2.0 * DI)
            # PFD = 2.5 * 4.0 * FFD * DEND * VD ** 2 * self.LD / (2.0 * DD)

            DENI = DENL * (1.0 - GV) + DENG * GV
            # DEND = DENL * (1.0 - (0.9*GF+0.1*GV)) + DENG * GF       # correction for GV in diffuser? try
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


            # FTI = 3.0
            # FTD = 3.0
            # FTI = 2.5
            # FTD = 2.5
            PTI = FTI * DENI * VI ** 2 / 2.0
            PTD = FTD * DEND * VD ** 2 / 2.0
            PPN = PEE - PFI - PFD - PTI - PTD

            UL = self.RLK * OMEGA
            PIO = PEE - PFI - PTI
            PLK = PIO - DENI * (U2 ** 2 - UL ** 2) / 8.0

            if PLK >= 0.:
                VL = np.abs(QLK) / (2.0 * pi * self.RLK * self.SL)
                REL = np.abs(DEND * VL * self.SL / VISD)    # changed VISL to VISD by Haiwen Zhu
                EDL = 0.0
                FFL = self.get_fff(REL, EDL)
                FFL = self.get_fff_leakage(REL,EDL,N, (self.R1+self.R2)/2, VL, self.LG, self.SL)        # by Haiwen Zhu
                VL = np.sqrt(2.0 * PLK / (1.5 + 4.0 * FFL * self.LG / self.SL) / DEND)
                QLKN = 2.0 * pi * self.RLK * self.SL * VL
            else:
                VL = np.abs(QLK / (2.0 * pi * self.RLK * self.SL))
                REL = DENL * VL * self.SL / VISL
                EDL = 0.0
                FFL = self.get_fff(REL, EDL)
                FFL = self.get_fff_leakage(REL,EDL,N, (self.R1+self.R2)/2, VL, self.LG, self.SL)        # by Haiwen Zhu
                VL = np.sqrt(2.0 * np.abs(PLK) / (1.5 + 4.0 * FFL * self.LG / self.SL) / DEND)
                QLKN = -2.0 * pi * self.RLK * self.SL * VL

            HLKloss = 0.25* DENI*(QLK/AI)**2        # by Haiwen Zhu
            PPN = PEE - PFI - PFD - PTI - PTD - HLKloss                              # by Haiwen Zhu
            ABQ = np.abs((QLKN - QLK) / QLK)
            QLK = QLKN
            ABP = np.abs((PPN - PP) / PPN)

            if icon > 200:
                break
                # PP = 0.9*PPN+0.1*PP
                # icon += 1
            else:
                PP = 0.5*PPN+0.5*PP
                icon += 1

        # return pressure in psi, flow rate in bpd
        PP = PP / psi_to_pa
        PE = PE / psi_to_pa
        PEE = PEE / psi_to_pa
        PF = (PFI + PFD) / psi_to_pa
        PT = (PTI + PTD) / psi_to_pa
        PD = (PFD + PTD) / psi_to_pa
        PRE = np.abs(PE - PEE)
        HLKloss = HLKloss / psi_to_pa
        QLK = QLK * 24.0 * 3600.0 / bbl_to_m3
        if FGL == 'BUB':
            FGL == 'BUB'
        # print('QL:', round(QL* 24.0 * 3600.0 / bbl_to_m3, 6), 'PP', round(PP,6), 'PRE', round(PRE,6), 'PF', round(PF,6), 'PT', round(PT,6),
        #              'PD', round(PD,6), 'HLKloss', round(HLKloss,6), 'GF', round(GF,6), 'GV', round(GV,6), 'FGL', FGL, 'icon', icon, 'CF', CF)
        return PP, PE, PF, PT, PD, PRE, PLK, QLK, GV

    def flow_pattern(self, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW, WC, flg):
        QL = np.arange(500)[1:] * 50.0
        QSG1 = []
        QSG2 = []
        QSG3 = []
        for ql in QL:
            
            ql = ql * bbl_to_m3 / 24 / 3600      # convert bpd to m3/s
            # ql in bpd, HP in psi
            if sgl_model=='zhang_2015':
                HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_old(ql, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW,
                                                                        WC)         
            elif sgl_model=='zhang_2016':
                HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_new(ql, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW,
                                                                        WC)   
            elif sgl_model=='jiecheng_2017':
                HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_jc(ql, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW,
                                                                        WC)   
            elif sgl_model=='zhu_2018':
                HP, HE, HF, HT, HD, HRE, HLK, QLK = self.sgl_calculate_2018(ql, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW,
                                                                   WC)   

            PP = HP * psi_to_pa

            QLK = QLK * bbl_to_m3 / 24 / 3600

            if HP < 0: break
            # LambdaC3
            LambdaC1 = self.get_lambda_c1(ql+QLK, PP , N)
            LambdaC2 = self.get_lambda_c2(ql+QLK, PP , N)
            LambdaC3 = self.get_lambda_c3(ql+QLK, N)
            # LambdaC1 = self.get_lambda_c1(ql+QLK, HP * psi_to_pa , N)
            # LambdaC2 = self.get_lambda_c2(ql+QLK, HP * psi_to_pa , N)
            # LambdaC3 = self.get_lambda_c3(ql+QLK, N)
            
            # LambdaC2
            ABP=1.
            icon = 0
            LambdaC2 = 0.01
            while ABP > E1:
                icon += 1
                if icon >200: break
                # PP = PP / psi_to_pa
                LambdaC2N = self.get_lambda_c2(ql+QLK, PP , N)
                QG2 = LambdaC2N / (1.0 - LambdaC2N) * ql
                PP,_,_,_,_,_,_,QLK,_ = self.gl_calculate_new(ql, QG2, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW, WC, 'B')
                PP = PP*psi_to_pa
                if PP < 0:
                    break
                QLK = QLK * bbl_to_m3 / 24 / 3600
                ABP = np.abs((LambdaC2-LambdaC2N)/LambdaC2)
                LambdaC2 = 0.5*LambdaC2N + 0.5*LambdaC2 
            if PP < 0:
                    break
            # LambdaC1
            ABP=1.
            icon = 0
            LambdaC1 = 0.01
            while ABP > E1:
                icon += 1
                if icon >200: break
                # PP = PP / psi_to_pa
                LambdaC1N = self.get_lambda_c1(ql+QLK, PP , N)
                QG1 = LambdaC1N / (1.0 - LambdaC1N) * ql
                PP,_,_,_,_,_,_,QLK,_ = self.gl_calculate_new(ql, QG1, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW, WC, 'B')
                PP = PP*psi_to_pa
                if PP < 0:
                    break
                QLK = QLK * bbl_to_m3 / 24 / 3600
                ABP = np.abs((LambdaC1-LambdaC1N)/LambdaC1) 
                LambdaC1 = 0.5*LambdaC1 + 0.5*LambdaC1N 
            if PP < 0:
                break

            if LambdaC1 > LambdaC2: LambdaC2 = LambdaC1
            if LambdaC2 > LambdaC3: LambdaC3 = LambdaC2

            QSG1.append(LambdaC1 / (1.0 - LambdaC1) * ql / (bbl_to_m3 / 24 / 3600))
            QSG2.append(LambdaC2 / (1.0 - LambdaC2) * ql / (bbl_to_m3 / 24 / 3600))
            QSG3.append(LambdaC3 / (1.0 - LambdaC3) * ql / (bbl_to_m3 / 24 / 3600))

        return QL, QSG1, QSG2, QSG3


#################################
class SinglePhaseCompare(object):
    def __init__(self, pump, conn):
        # QBEM = {'TE2700': 5600, 'DN1750': 3300, 'GC6100': 8800, 'P100': 12000}      # sgl_new
        # QBEM = {'TE2700': 4500, 'DN1750': 3000, 'GC6100': 7800, 'P100': 11000, 'Flex31': 5000}    # sgl_2018
        QBEM = QBEM_default
        self.pump = pump
        self.ESP = ESP[pump]
        self.conn = conn
        self.df_catalog = pd.read_sql_query("SELECT * FROM Catalog_All;", self.conn)
        self.df_pump = self.df_catalog[self.df_catalog['Pump'] == pump]
        # print(self.df_pump)
        self.QBEM = QBEM[pump]

    def single_phase_water(self, sgl_model):
        sgl = SinglePhaseModel(self.ESP, self.QBEM)
        QL, hpsgl, _, _, _, _, _, _, _ = sgl.performance_curve(sgl_model)       #four model options
        # fig = plt.figure(dpi=300)
        # ax = fig.add_subplot(111)
        # ax.plot(QL, hpsgl, 'b-', label='model')
        # ax.plot(self.df_pump['Flow_bpd'], self.df_pump['DP_psi'], 'ro', label='catalog')
        # ax.set_xlabel(r'$Q_L$ (bpd)')
        # ax.set_ylabel(r'$P$ (psi)')
        # ax.set_title(r'{} ESP, $Q_B$={} bpd'.format(self.pump, self.QBEM))
        # ax.legend(frameon=False)
        # 
        df = pd.DataFrame({'ql': QL, 'hp': hpsgl})
        return df

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
            HP, _, _, _, _, _, _, _ = sgl.sgl_calculate_new(ql, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)
            #HP, _, _, _, _, _, _, _ = sgl.sgl_calculate_2018(ql, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)
            if HP > 0:
                hpsgl.append(HP)
            else:
                break
        df_model = pd.DataFrame({'QL': QL[:len(hpsgl)], 'DP': hpsgl})
        return df_model

#################################
class TwoPhaseCompare(object):
    def __init__(self, pump, conn):
        # QBEM = {'TE2700': 6000, 'DN1750': 3300, 'GC6100': 8800, 'P100': 12000, 'Flex31': 6200}  # gl_new
        QBEM = QBEM_default
        QBEP = {'TE2700': 2700, 'DN1750': 1750, 'GC6100': 6100, 'P100': 9000, 'Flex31': 3100}  # bep flow rate
        self.pump = pump
        self.ESP = ESP[pump]
        self.conn = conn
        self.QBEM = QBEM[pump]
        self.QBEP = QBEP[pump]

    def surging_cd_to_db(self):
        """
        :return: two dataframes for GV and CD_over_dB, the column names: zhu, Barrios, sun
        """
        sgl = SinglePhaseModel(self.ESP, self.QBEM)
        gl = GasLiquidModel(self.ESP, self.QBEM)
        sgl_cal = np.vectorize(sgl.sgl_calculate_new)      #four options: old (zhu and zhang), new (Dr. Zhang update), jc (jiecheng), 2018 (jiecheng-jianjun)
        gl_cal = np.vectorize(gl.gl_calculate_new)  #old, new
        zhu_cal = np.vectorize(gl.gv_zhu)
        sun_cal = np.vectorize(gl.gv_sun)
        Barrios_cal = np.vectorize(gl.gv_Barrios)

        # use the best efficiency point to compare
        GF =  np.arange(0.01, 1, 0.05)
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
        gv_sun, cd_over_db_sun = sun_cal(B1, B2, DENL, DENG, N, QL, QG, QLK, R1, R2, TB, VISL, VISG, YI1, YI2, YI, ZI)
        gv_Barrios, cd_over_db_Barrios = Barrios_cal(DENL, DENG, GF, N, QL, QLK, R1, RI, ST, TB, VISL, YI, ZI)

        df_gv = pd.DataFrame({'zhu': gv_zhu, 'sun': gv_sun, 'Barrios': gv_Barrios, 'gf': GF})
        df_cd_over_db = pd.DataFrame({'zhu': cd_over_db_zhu, 'sun': cd_over_db_sun,
                                      'Barrios': cd_over_db_Barrios, 'gf': GF})
        
        # flg for different drag coefficient model
        flgz = np.empty(GF.shape, dtype='str')
        flgs = np.empty(GF.shape, dtype='str')
        flgb = np.empty(GF.shape, dtype='str')
        flgz[:] = 'F'
        flgb[:] = 'B'
        flgs[:] = 'S'

        HP, _, _, _, _, _, _, _ = sgl_cal(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)

        PPZ, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgz)
        # df = pd.DataFrame({'gf': GF, 'zhu': PPZ/HP})

        PPS, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgs)
        #df = pd.DataFrame({'gf': GF, 'sun': PPS/HP})
        
        PPB, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM,
                                             SN, ST, VISG, VISL, VISW, WC, flgb)
        # df = pd.DataFrame({'gf': GF, 'Barrios': PPB/HP})
        

        df = pd.DataFrame({'gf': GF, 'zhu': PPZ/HP, 'sun': PPS/HP, 'Barrios': PPB/HP})
        
        return df_gv, df_cd_over_db, df

    def surging_performance(self, QL, maxGF, N, p, t):
        """
        :param QL: liquid flow rate in bpd
        :param maxGF: maximum GF for calculation
        :param N: array for rotational speed rpm
        :param p: array for gas pressure psi
        :param t: array for temperature F
        :return: dataframe of predicted pump heads under surging flow, the column names: zhu, Barrios, sun
        """
        sgl = SinglePhaseModel(self.ESP, self.QBEM)
        gl = GasLiquidModel(self.ESP, self.QBEM)
        sgl_cal = np.vectorize(sgl.sgl_calculate_new) #four options: old (zhu and zhang), new (Dr. Zhang update), jc (jiecheng), 2018 (jiecheng-jianjun)
        gl_cal = np.vectorize(gl.gl_calculate_new)  #old, new

        GF = np.arange(0.0, maxGF + 0.02, 0.02)
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
        flgz[:] = 'F'
        flgb[:] = 'B'
        flgs[:] = 'S'

        HP, _, _, _, _, _, _, _ = sgl_cal(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)

        PPZ, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgz)
        # df = pd.DataFrame({'gf': GF, 'zhu': PPZ/HP})

        # PPS, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
        #                                      WC, flgs)
        #df = pd.DataFrame({'gf': GF, 'sun': PPS/HP})
        
        # PPB, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM,
        #                                      SN, ST, VISG, VISL, VISW, WC, flgb)
        # df = pd.DataFrame({'gf': GF, 'Barrios': PPB/HP})
        

        # df = pd.DataFrame({'gf': GF, 'zhu': PPZ/HP, 'sun': PPS/HP, 'Barrios': PPB/HP})
        # df = pd.DataFrame({'gf': GF, 'zhu': PPZ, 'sun': PPS, 'Barrios': PPB})
        df = pd.DataFrame({'gf': GF, 'zhu': PPZ})

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
        flgz = np.empty(QL.shape, dtype='str')
        flgs = np.empty(QL.shape, dtype='str')
        flgb = np.empty(QL.shape, dtype='str')
        flgz[:] = 'F'
        flgb[:] = 'B'
        flgs[:] = 'S'

        PPZ, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgz)
        # PPS, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
        #                                      WC, flgs)
        # PPB, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM,
        #                                      SN, ST, VISG, VISL, VISW, WC, flgb)

        # df = pd.DataFrame({'zhu': PPZ, 'sun': PPS, 'Barrios': PPB})
        df = pd.DataFrame({'zhu': PPZ})
        return df

    def mapping_performance(self, QG, maxQL, N, p, t):
        """
        :param QG: constant gas flow rate bpd
        :param maxQL: maximum liquid flow rate bpd
        :param N: rotational speed rpm
        :param p: array for gas pressure psi
        :param t: array for temperature F
        :return: dataframe of predicted pump heads under mapping flow, the column names: zhu, Barrios, sun
        """
        gl = GasLiquidModel(self.ESP, self.QBEM)
        gl_cal = np.vectorize(gl.gl_calculate_new)

        # if QG/maxQL<0.02:
        #     # one point wrong, no idea, avoid the point
        #     QL1 =  np.arange(0.01, 0.34, 0.02) * maxQL * bbl_to_m3 / 24.0 / 3600.0
        #     QL2 =  np.arange(0.34, 1.2, 0.08) * maxQL * bbl_to_m3 / 24.0 / 3600.0
        #     QL = np.hstack((QL1,QL2))   # horizontal combine
        #     # QL = np.vstack((QL1,QL2))   # vertical combine
        # else:
        #     QL =  np.arange(0.01, 1.1, 0.02) * maxQL * bbl_to_m3 / 24.0 / 3600.0
        QL =  np.arange(0.01, 1.1, 0.02) * maxQL * bbl_to_m3 / 24.0 / 3600.0
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
        flgzz = np.empty(QL.shape, dtype='str')
        flgh = np.empty(QL.shape, dtype='str')
        flgz[:] = 'F'   # new Zhu model
        flgb[:] = 'B'   # Barrios
        flgs[:] = 'S'   # Sun
        flgzz[:] = 'Z'  # old Zhu model
        flgh[:] = 'H'   # homogenous model

        PPZ, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgz)
        PPS, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgs)
        PPB, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM,
                                             SN, ST, VISG, VISL, VISW, WC, flgb)
        PPZZ, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgzz)
        PPH, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgh)
        # df = pd.DataFrame({'ql': QL/bbl_to_m3 * 3600 * 24, 'zhu': PPZ, 'sun': PPS, 'Barrios': PPB})
        # df = pd.DataFrame({'ql': QL/bbl_to_m3 * 3600 * 24, 'zhu': PPZ})
        df = pd.DataFrame({'ql': QL/bbl_to_m3 * 3600 * 24, 'zhu': PPZ, 'sun': PPS, 'Barrios': PPB, 'Old_zhu': PPZZ, 'Homo': PPH})

        return df

    def GVF_performance(self, GF, maxQL, minQL, N, p, t):
        """
        :param QG: constant gas flow rate bpd
        :param maxQL: maximum liquid flow rate bpd
        :param N: rotational speed rpm
        :param p: array for gas pressure psi
        :param t: array for temperature F
        :return: dataframe of predicted pump heads under mapping flow, the column names: zhu, Barrios, sun
        """
        
        sgl = SinglePhaseModel(self.ESP, self.QBEM)
        gl = GasLiquidModel(self.ESP, self.QBEM)
        sgl_cal = np.vectorize(sgl.sgl_calculate_new)
        gl_cal = np.vectorize(gl.gl_calculate_new)
        QLmin = minQL/maxQL
        if QLmin == 0:
            QLmin = 0.01
        QL = np.arange(0.3*QLmin, 1.2, 0.02) * maxQL * bbl_to_m3 / 24.0 / 3600.0
        QG = GF / (1 - GF) * QL

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
        flgz[:] = 'F'
        flgb[:] = 'B'    
        flgs[:] = 'S'
        
        HP, _, _, _, _, _, _, _ = sgl_cal(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)

        PPZ, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
                                             WC, flgz)
        # df = pd.DataFrame({'ql': QL/bbl_to_m3 * 3600 * 24, 'zhu': PPZ/HP})

        # PPS, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM, SN, ST, VISG, VISL, VISW,
        #                                      WC, flgs)
        # #df = pd.DataFrame({'ql': QL/bbl_to_m3 * 3600 * 24, 'sun': PPS, 'water': HP})
        
        # PPB, _, _, _, _, _, _, _, _ = gl_cal(QL, QG, QBEM, DENG, DENL, DENW, N, NS, SGM,
        #                                      SN, ST, VISG, VISL, VISW, WC, flgb)
        # df = pd.DataFrame({'ql': QL/bbl_to_m3 * 3600 * 24, 'Barrios': PPB/HP})
        

        df = pd.DataFrame({'ql': QL/bbl_to_m3 * 3600 * 24, 'zhu': PPZ})
        # df = pd.DataFrame({'ql': QL/bbl_to_m3 * 3600 * 24, 'zhu': PPZ, 'sun': PPS, 'Barrios': PPB, 'water': HP})
        #df = pd.DataFrame({'ql': QL / bbl_to_m3 * 3600 * 24, 'sun': PPS/HP})

        return df

# fitting water curve
def fitting_water(pump_name, QBEM, sgl_model):
    conn, c = connect_db('ESP.db')
    df_data = pd.read_sql_query("SELECT * FROM Catalog_All;", conn)
    df_data = df_data[df_data.Pump == pump_name]
    df_data = df_data[df_data.Flow_bpd != 0]
    df_data=df_data.reset_index(drop=True)
    
    # fig=plt.figure(dpi=128, figsize=(10,6))
    # ax = fig.add_subplot(111)
    # ax.scatter(df_data.Flow_bpd, df_data.DP_psi,label='test', color = 'red')

    sgl = SinglePhaseModel(ESP[pump_name],QBEM)

    if sgl_model=='zhang_2015':
        sgl_cal = np.vectorize(sgl.sgl_calculate_old)           
    elif sgl_model=='zhang_2016':
        sgl_cal = np.vectorize(sgl.sgl_calculate_new)    
    elif sgl_model=='jiecheng_2017':
        sgl_cal = np.vectorize(sgl.sgl_calculate_jc) 
    elif sgl_model=='zhu_2018':
        sgl_cal = np.vectorize(sgl.sgl_calculate_2018)

    ABV = 1
    error = 1
    icon = 0
    QL = df_data.Flow_bpd * bbl_to_m3 / 24.0 / 3600.0
    QBEM = QBEM * bbl_to_m3 / 24.0 / 3600.0 # * np.ones(QL.shape)
    DENL = DENW # * np.ones(QL.shape)
    N = ESP[pump_name]['N'] # * np.ones(QL.shape)
    NS = ESP[pump_name]['NS'] # * np.ones(QL.shape)
    SGM = ESP[pump_name]['SGM'] # * np.ones(QL.shape)
    SN = ESP[pump_name]['SN']
    ST = ESP[pump_name]['ST']
    VISL = ESP[pump_name]['VISW']
    VISW = ESP[pump_name]['VISW']
    WC = ESP[pump_name]['WC']
    # HP = df_data.DP_psi
    HP = np.ones(df_data.Flow_bpd.shape)
    for i in range(HP.shape[0]):
        HP[i] = df_data.DP_psi[i]


    while np.abs(ABV) > E1 or np.abs(error) >E1:
        error = 0
        ABV = 0
        icon += 1
        if icon > 1000:
            break
        HPN, _, _, _, _, _, _, _ = sgl_cal(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)

        for i in range(HP.shape[0]):
            error += (HPN[i]-df_data.DP_psi[i])/(df_data.DP_psi[i]) 
            ABV += (HP[i]-HPN[i])/HP[i]
            HP[i] = HPN[i]
        if error >0:
            if QBEM*0.1 < error* bbl_to_m3 / 24.0 / 3600.0:
                QBEM = QBEM * 1.1
            else:
                QBEM = QBEM + error* bbl_to_m3 / 24.0 / 3600.0
        else:
            if QBEM*0.1 > np.abs(error* bbl_to_m3 / 24.0 / 3600.0):
                QBEM = QBEM * 0.9
            else:
                QBEM = QBEM + error* bbl_to_m3 / 24.0 / 3600.0

    QL = QL / (bbl_to_m3 / 24.0 / 3600.0)
    QBEM = QBEM / (bbl_to_m3 / 24.0 / 3600.0)
    
    # df_model = ESP_case.single_phase_water(sgl_model)
    # ax.plot(QL, HP, label='model', color = 'blue')
    # ax.set_title(pump_name + ' catalog at ''QBEM= '+ str(QBEM))
    # ax.legend(frameon=False, fontsize=5)
    # 
    return QBEM, QL, HP

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
        #fig.show()
    fig.savefig('VisTe2700.jpg')

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
    ax3.set_ylabel(r'$P_{sim}$ (psi)')
    ax3.legend(frameon=False, fontsize=8)
    fig3.show()
    fig3.savefig('VisTe2700Error.jpg')
    
    # epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_jc.model, df_jc.DP_psi)
    # print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'##########')

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
    # print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'##########')

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
    print(df)
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
        ax1.scatter(df_sub.Flowrate_bpd, df_sub.Stage3_DeltaP_psi, facecolor='none', edgecolor=colors[index],
                    linewidth=0.5, marker=symbols[index], label=r'$N$={} rpm'.format(rpm))

        # ax1.scatter(df_sub.Stage3_DeltaP_psi, df_sub.model, facecolor='none', edgecolor=colors[index],
        #             linewidth=0.5, marker=symbols[index], label=r'$N$={} rpm'.format(rpm))

    # ax1.plot(x_max, x_max, 'k--', label='perfect match')
    # ax1.plot(x_max, 1.25 * x_max, ls=':', c='gray')
    # ax1.plot(x_max, 0.7 * x_max, ls=':', c='gray')
    # ax1.set_xlabel(r'$P_{exp}$ (psi)')
    # ax1.set_ylabel(r'$P_{sim}$ (psi)')

    ax1.legend(frameon=False, fontsize=7)
    fig1.show()
    

    # epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df.model, df.Stage3_DeltaP_psi)
    # print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'##########')

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
            ax.plot(df_turzo.Qvis, df_turzo.DPvis, c='C{}'.format(id + 1), label=r'{} cp_Turzo (2000)'.format(visl),
                     linestyle='--', linewidth=0.75)
            id += 1

        handles, labels = ax.get_legend_handles_labels()
        ax.set_xlabel(r'$Q_L$ (bpd)')
        ax.set_ylabel(r'$P$ (psi)')
        ax.legend(handles, labels, frameon=False, fontsize=4)
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
    ax3.set_ylabel(r'$P_{sim}$ (psi)')
    ax3.legend(frameon=False, fontsize=8)
    fig3.show()
    # epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = \
    #     stats_analysis(af * df_oil_only.model, df_oil_only.DPStage3_psi)
    # print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'##########')

    # plot emulsion error line
    df_emulsion = df_dn1750[(df_dn1750.oil_viscosity_cP != 0) & (df_dn1750.Oil != 'ND20') & (df_dn1750.water_cut != 0)]
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
    x_max = np.array(range(int(max(df_oil_only.model) + 1)))

    ax4.scatter(df_emulsion.DPStage3_psi, af * df_emulsion.model, facecolor='none', edgecolor='b', label='')
    ax4.set_title('Emulsion', fontsize=8)
    ax4.plot(x_max, x_max, 'k--', label='perfect match')
    ax4.plot(x_max, 1.25 * x_max, 'r:')
    ax4.plot(x_max, 1.50 * x_max, 'r:')
    ax4.plot(x_max, 0.85 * x_max, 'r:')
    ax4.set_xlabel(r'$P_{exp}$ (psi)')
    ax4.set_ylabel(r'$P_{sim}$ (psi)')
    ax4.legend(frameon=False)
    # fig4.show()
    
    # epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = \
    #     stats_analysis(af * df_emulsion.model, df_emulsion.DPStage3_psi)
    # print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'##########')

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
    # QBEM = {'TE2700': 4500, 'DN1750': 3000, 'GC6100': 7800, 'P100': 11000}
    QBEM = QBEM_default
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

def GVF_test():
    conn, c = connect_db('ESP.db')

    te2700_case = TwoPhaseCompare('TE2700', conn)
    df_1GVF = te2700_case.GVF_performance(0.001, 10000, 0, 3500, 100, 100)
    df_2GVF = te2700_case.GVF_performance(0.05, 10000, 0, 3500, 100, 100)
    df_5GVF = te2700_case.GVF_performance(0.10, 10000, 0, 3500, 100, 100)
    
    dfp = pd.DataFrame({'water': df_1GVF.water,
                        'dp_1zhu': df_1GVF.zhu,
                        'dp_2zhu': df_2GVF.zhu,
                        'dp_5zhu': df_5GVF.zhu,
                        'dp_1sun': df_1GVF.sun,
                        'dp_2sun': df_2GVF.sun,
                        'dp_5sun': df_5GVF.sun,
                        'dp_1Barrios': df_1GVF.Barrios,
                        'dp_2Barrios': df_2GVF.Barrios,
                        'dp_5Barrios': df_5GVF.Barrios,
                        'gf_water': df_1GVF.ql,
                        'gf_1GVF': df_1GVF.ql,
                        'gf_2GVF': df_2GVF.ql,
                        'gf_5GVF': df_5GVF.ql})

   
    fig1 = plt.figure(dpi=300)
    ax1 = fig1.add_subplot(111)
    ax1.plot(dfp['gf_water'] , dfp['dp_1zhu'], 'mo', label='1%')
    ax1.plot(dfp['gf_1GVF'] , dfp['dp_2zhu'], 'ks', label='5%')
    ax1.plot(dfp['gf_2GVF'] , dfp['dp_5zhu'], 'r^', label='10%')
    ax1.plot(dfp['gf_5GVF'] , dfp['water'], 'g*', label='water')

    ax1.set_xlabel(r'$Q_{L}$ (bpd)')
    ax1.set_ylabel(r'$P$ (psi)')
    ax1.set_xlim(0, 10000)
    ax1.set_ylim(0, 30)
    ax1.legend(frameon=False)
    fig1.show()
    fig1.savefig('GVFfig1zhu.jpg')

    fig2 = plt.figure(dpi=300)
    ax2 = fig2.add_subplot(111)
    ax2.plot(dfp['gf_water'] , dfp['dp_1sun'], 'mo', label='1%')
    ax2.plot(dfp['gf_1GVF'] , dfp['dp_2sun'], 'ks', label='5%')
    ax2.plot(dfp['gf_2GVF'] , dfp['dp_5sun'], 'r^', label='10%')
    ax2.plot(dfp['gf_5GVF'] , dfp['water'], 'g*', label='water')

    ax2.set_xlabel(r'$Q_{L}$ (bpd)')
    ax2.set_ylabel(r'$P$ (psi)')
    ax2.set_xlim(0, 10000)
    ax2.set_ylim(0, 30)
    ax2.legend(frameon=False)
    fig2.show()
    fig2.savefig('GVFfig2sun.jpg')

    fig3 = plt.figure(dpi=300)
    ax3 = fig3.add_subplot(111)
    ax3.plot(dfp['gf_water'] , dfp['dp_1Barrios'], 'mo', label='1%')
    ax3.plot(dfp['gf_1GVF'] , dfp['dp_2Barrios'], 'ks', label='5%')
    ax3.plot(dfp['gf_2GVF'] , dfp['dp_5Barrios'], 'r^', label='10%')
    ax3.plot(dfp['gf_5GVF'] , dfp['water'], 'g*', label='water')

    ax3.set_xlabel(r'$Q_{L}$ (bpd)')
    ax3.set_ylabel(r'$P$ (psi)')
    ax3.set_xlim(0, 10000)
    ax3.set_ylim(0, 30)
    ax3.legend(frameon=False)
    fig3.show()
    fig3.savefig('GVFfig3Barrios.jpg')
    

    disconnect_db(conn)

def plot_flow_pattern(pump_name, QBEM):
    ESP_case=GasLiquidModel(ESP[pump_name], QBEM)
    QBEM = QBEM * bbl_to_m3 / 24.0 / 3600.0
    fig = plt.figure(dpi =128, figsize = (6,4.5))
    ax = fig.add_subplot(1,1,1)
    fig2, ax2 = plt.subplots(dpi = 128, figsize = (6,4.5))

    QSL, QSG1, QSG2, QSG3 = ESP_case.flow_pattern(QBEM, ESP[pump_name]['DENG'], ESP[pump_name]['DENL'], ESP[pump_name]['DENW'], \
                    ESP[pump_name]['N'], ESP[pump_name]['NS'], ESP[pump_name]['SGM'], ESP[pump_name]['SN'], ESP[pump_name]['ST'], \
                        ESP[pump_name]['VISG'], ESP[pump_name]['VISL'], ESP[pump_name]['VISW'], ESP[pump_name]['WC'], 'F')
    GVF1 = QSG1/(QSG1+QSL[:len(QSG1)])
    GVF2 = QSG2/(QSG2+QSL[:len(QSG2)])
    GVF3 = QSG3/(QSG3+QSL[:len(QSG3)])
    # QSG vs. QSL
    ax.plot(QSG1, QSL[:len(QSG1)], label='DB-BUB', color='red')
    ax.plot(QSG2, QSL[:len(QSG2)], label='BUB-INT', color='yellow')
    ax.plot(QSG3, QSL[:len(QSG3)], label='INT-ANN', color='green')
    # GVF vs. QSL
    ax2.plot(GVF1, QSL[:len(QSG1)], label='DB-BUB', color='red')
    ax2.plot(GVF2, QSL[:len(QSG2)], label='BUB-INT', color='yellow')
    ax2.plot(GVF3, QSL[:len(QSG3)], label='INT-ANN', color='green')
    ax.legend(frameon=False)
    ax2.legend(frameon=False)
    # ax.set_xlim(0,0.2*QSL[:len(QSG1)].max())

def Flex31_surging_3100bpd(SQLname,tablename,casename, casename2, casename3, TargetQL, xlabel, ylabel, title, figurename, pump_name, test_type):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)

    df_data=df_data[((df_data.Case == casename) | (df_data.Case == casename2) | (df_data.Case == casename3)) & (df_data.TargetQL_bpd == TargetQL)]
    df_data=df_data.reset_index(drop=True)

    fig = plt.figure (dpi =128, figsize = (6,4.5))
    ax = fig.add_subplot(111)
    


    ESP_case = TwoPhaseCompare(pump_name, conn)

    # 1st surging 190 psi 3100bpd
    if casename != 'none':
        df_surging=df_data[df_data.Case == casename]

        ax.scatter(df_surging['HG_%'],df_surging['dpAve_psi'],label=casename,color='red')
        df_model1 = ESP_case.surging_performance(3100, 0.18, 3600, 190, 40)
        # df_model1 = ESP_case.surging_performance(df_surging.Qw_bpd.mean(), df_surging['HG_%'].max()/100., df_surging.RPM.mean(), df_surging.Ptank_psi.mean(),
        #                                         df_surging.Tin_F.mean())
        ax.plot(df_model1.gf * 100, df_model1.sun, 'r--', label='Sun'+casename)
        ax.plot(df_model1.gf * 100, df_model1.Barrios, 'r-.', label='Barrios'+casename)
        ax.plot(df_model1.gf * 100, df_model1.zhu, 'r-', label='Zhu'+casename)

    # 2st surging 150 psi 3100bpd
    if casename2 != 'none':
        df_surging=df_data[df_data.Case == casename2]

        ax.scatter(df_surging['HG_%'],df_surging['dpAve_psi'],label=casename2,color='b')
        df_model1 = ESP_case.surging_performance(3100, 0.18, 3600, 150, 40)
        # df_model1 = ESP_case.surging_performance(df_surging.Qw_bpd.mean(), df_surging['HG_%'].max()/100., df_surging.RPM.mean(), df_surging.Ptank_psi.mean(),
        #                                         df_surging.Tin_F.mean())
        ax.plot(df_model1.gf * 100, df_model1.sun, 'b--', label='Sun'+casename2)
        ax.plot(df_model1.gf * 100, df_model1.Barrios, 'b-.', label='Barrios'+casename2)
        ax.plot(df_model1.gf * 100, df_model1.zhu, 'b-', label='Zhu'+casename2)

    # 3st surging 35 psi 3100bpd
    if casename3 != 'none':
        df_surging=df_data[df_data.Case == casename3]
        ax.scatter(df_surging['HG_%'],df_surging['dpAve_psi'],label=casename3,color='c')
        df_model1 = ESP_case.surging_performance(3100, 0.18, 3600, 35, 40)
        # df_model1 = ESP_case.surging_performance(df_surging.Qw_bpd.mean(), df_surging['HG_%'].max()/100., df_surging.RPM.mean(), df_surging.Ptank_psi.mean(),
        #                                         df_surging.Tin_F.mean())
        ax.plot(df_model1.gf * 100, df_model1.sun, 'c--', label='Sun'+casename3)
        ax.plot(df_model1.gf * 100, df_model1.Barrios, 'c-.', label='Barrios'+casename3)
        ax.plot(df_model1.gf * 100, df_model1.zhu, 'c-', label='Zhu'+casename3)


    ax.set_xlabel(r'$\lambda_{G}$ (%)')
    ax.set_ylabel(r'$N_{P}$')
    ax.legend(frameon=False)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax.set_zlabel(zlabel)
    ax.set_title(title)
    
    fig.savefig(figurename)

    disconnect_db(conn)

def Flex31_surging_35psi(SQLname,tablename,casename, TargetQL1, TargetQL2, TargetQL3, xlabel, ylabel, title, figurename, pump_name, test_type):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)

    df_data=df_data[((df_data.TargetQL_bpd == TargetQL1) | (df_data.TargetQL_bpd == TargetQL2) | (df_data.TargetQL_bpd == TargetQL3)) & (df_data.Case == casename)]
    df_data=df_data.reset_index(drop=True)

    fig = plt.figure (dpi =128, figsize = (6,4.5))
    ax = fig.add_subplot(111)
    


    ESP_case = TwoPhaseCompare(pump_name, conn)

    # 1st 
    if TargetQL1 != 'none':
        df_surging=df_data[df_data.TargetQL_bpd == TargetQL1]

        ax.scatter(df_surging['HG_%'],df_surging['dpAve_psi'],label='Test'+str(TargetQL1),color='red')
        df_model1 = ESP_case.surging_performance(df_surging.Qw_bpd.mean(), df_surging['HG_%'].max()/100., df_surging.RPM.mean(), df_surging.Ptank_psi.mean(),
                                                df_surging.Tin_F.mean())
        # ax.plot(df_model1.gf * 100, df_model1.sun, 'r--', label='Sun'+str(TargetQL1))
        # ax.plot(df_model1.gf * 100, df_model1.Barrios, 'r-.', label='Barrios'+str(TargetQL1))
        ax.plot(df_model1.gf * 100, df_model1.zhu, 'r-', label='Zhu'+str(TargetQL1))

    # 2st 
    if TargetQL2 != 'none':
        df_surging=df_data[df_data.TargetQL_bpd == TargetQL2]

        ax.scatter(df_surging['HG_%'],df_surging['dpAve_psi'],label='Test'+str(TargetQL2),color='b')
        df_model1 = ESP_case.surging_performance(df_surging.Qw_bpd.mean(), df_surging['HG_%'].max()/100., df_surging.RPM.mean(), df_surging.Ptank_psi.mean(),
                                                df_surging.Tin_F.mean())
        # ax.plot(df_model1.gf * 100, df_model1.sun, 'b--', label='Sun'+str(TargetQL2))
        # ax.plot(df_model1.gf * 100, df_model1.Barrios, 'b-.', label='Barrios'+str(TargetQL2))
        ax.plot(df_model1.gf * 100, df_model1.zhu, 'b-', label='Zhu'+str(TargetQL2))

    # 3st 
    if TargetQL3 != 'none':
        df_surging=df_data[df_data.TargetQL_bpd == TargetQL3]
        ax.scatter(df_surging['HG_%'],df_surging['dpAve_psi'],label='Test'+str(TargetQL3),color='c')
        df_model1 = ESP_case.surging_performance(df_surging.Qw_bpd.mean(), df_surging['HG_%'].max()/100., df_surging.RPM.mean(), df_surging.Ptank_psi.mean(),
                                                df_surging.Tin_F.mean())
        # ax.plot(df_model1.gf * 100, df_model1.sun, 'c--', label='Sun'+str(TargetQL3))
        # ax.plot(df_model1.gf * 100, df_model1.Barrios, 'c-.', label='Barrios'+str(TargetQL3))
        ax.plot(df_model1.gf * 100, df_model1.zhu, 'c-', label='Zhu'+str(TargetQL3))


    ax.set_xlabel(r'$\lambda_{G}$ (%)')
    ax.set_ylabel(r'$N_{P}$')
    ax.legend(frameon=False)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax.set_zlabel(zlabel)
    ax.set_title(title)
    
    fig.savefig(figurename)

    disconnect_db(conn)

def plot_3d_surface(SQLname,tablename,casename, time1, time2, time3, xlabel, ylabel, zlabel, title, figurename):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)

    df_data=df_data[((df_data.Time == time1) | (df_data.Time == time2) | (df_data.Time == time3)) & (df_data.Case == casename)]
    df_data=df_data.sort_values(by=['Time'])
    df_data=df_data.reset_index(drop=True)

    fig = plt.figure (dpi =128, figsize = (6,4.5))
    ax = Axes3D(fig)
    
    HG = np.empty([3,15], dtype = float) 
    Qw = np.empty([3,15], dtype = float) 
    dpave = np.empty([3,15], dtype = float)
    dp3 = np.empty([3,15], dtype = float)
    dp6 = np.empty([3,15], dtype = float)
    dp9 = np.empty([3,15], dtype = float)
    dp12 = np.empty([3,15], dtype = float)
    Time = np.empty([3,15], dtype = float)
    i=0
    j=0
    HGlist=[]
    Qwlist=[]
    dpavelist=[]
    dp3list=[]
    dp6list=[]
    dp9list=[]
    dp12list=[]

    for k in range(df_data.shape[0]):
        HGlist.append(df_data['HG_%'][k])
        Qwlist.append(df_data.Qw_bpd[k])
        dpavelist.append(df_data.dpAve_psi[k])
        dp3list.append(df_data.dp3_psi[k])
        dp6list.append(df_data.dp6_psi[k])
        dp9list.append(df_data.dp9_psi[k])
        dp12list.append(df_data.dp12_psi[k])
        if k==df_data.shape[0]-1:
            if xlabel == 'GVF (%)':
                parameter1 = np.polyfit(HGlist,dpavelist,2)
                f1 = np.poly1d(parameter1)
                parameter2 = np.polyfit(HGlist,dp3list,2)
                f2 = np.poly1d(parameter2)
                parameter3 = np.polyfit(HGlist,dp6list,2)
                f3 = np.poly1d(parameter3)
                parameter4 = np.polyfit(HGlist,dp9list,2)
                f4 = np.poly1d(parameter4)
                parameter5 = np.polyfit(HGlist,dp12list,2)
                f5 = np.poly1d(parameter5)
                for j in range(15):
                    HG[i][j]=(j)/14*np.max(df_data['HG_%'])
                    Qw[i][j]=(j)/14*(np.max(df_data.Qw_bpd)-np.min(df_data.Qw_bpd))+np.min(df_data.Qw_bpd)
                    dpave[i][j]=f1(HG[i][j])
                    dp3[i][j]=f2(HG[i][j])
                    dp6[i][j]=f3(HG[i][j])
                    dp9[i][j]=f4(HG[i][j])
                    dp12[i][j]=f5(HG[i][j])
                    Time[i][j]=df_data.Time[k]
            else:
                parameter1 = np.polyfit(Qwlist,dpavelist,2)
                f1 = np.poly1d(parameter1)
                parameter2 = np.polyfit(Qwlist,dp3list,2)
                f2 = np.poly1d(parameter2)
                parameter3 = np.polyfit(Qwlist,dp6list,2)
                f3 = np.poly1d(parameter3)
                parameter4 = np.polyfit(Qwlist,dp9list,2)
                f4 = np.poly1d(parameter4)
                parameter5 = np.polyfit(Qwlist,dp12list,2)
                f5 = np.poly1d(parameter5)
                for j in range(15):
                    HG[i][j]=(j)/14*np.max(df_data['HG_%'])
                    Qw[i][j]=(j)/14*(np.max(df_data.Qw_bpd)-np.min(df_data.Qw_bpd))+np.min(df_data.Qw_bpd)
                    dpave[i][j]=f1(Qw[i][j])
                    dp3[i][j]=f2(Qw[i][j])
                    dp6[i][j]=f3(Qw[i][j])
                    dp9[i][j]=f4(Qw[i][j])
                    dp12[i][j]=f5(Qw[i][j])
                    Time[i][j]=df_data.Time[k]
            HGlist.clear()
            Qwlist.clear()
            dpavelist.clear()
            dp3list.clear()
            dp6list.clear()
            dp9list.clear()
            dp12list.clear()
        else:
            if df_data.Time[k+1] != df_data.Time[k]:
                if xlabel == 'GVF (%)':
                    parameter1 = np.polyfit(HGlist,dpavelist,2)
                    f1 = np.poly1d(parameter1)
                    parameter2 = np.polyfit(HGlist,dp3list,2)
                    f2 = np.poly1d(parameter2)
                    parameter3 = np.polyfit(HGlist,dp6list,2)
                    f3 = np.poly1d(parameter3)
                    parameter4 = np.polyfit(HGlist,dp9list,2)
                    f4 = np.poly1d(parameter4)
                    parameter5 = np.polyfit(HGlist,dp12list,2)
                    f5 = np.poly1d(parameter5)
                    for j in range(15):
                        HG[i][j]=(j)/14*np.max(df_data['HG_%'])
                        Qw[i][j]=(j)/14*(np.max(df_data.Qw_bpd)-np.min(df_data.Qw_bpd))+np.min(df_data.Qw_bpd)
                        dpave[i][j]=f1(HG[i][j])
                        dp3[i][j]=f2(HG[i][j])
                        dp6[i][j]=f3(HG[i][j])
                        dp9[i][j]=f4(HG[i][j])
                        dp12[i][j]=f5(HG[i][j])
                        Time[i][j]=df_data.Time[k]
                else:
                    parameter1 = np.polyfit(Qwlist,dpavelist,2)
                    f1 = np.poly1d(parameter1)
                    parameter2 = np.polyfit(Qwlist,dp3list,2)
                    f2 = np.poly1d(parameter2)
                    parameter3 = np.polyfit(Qwlist,dp6list,2)
                    f3 = np.poly1d(parameter3)
                    parameter4 = np.polyfit(Qwlist,dp9list,2)
                    f4 = np.poly1d(parameter4)
                    parameter5 = np.polyfit(Qwlist,dp12list,2)
                    f5 = np.poly1d(parameter5)
                    for j in range(15):
                        HG[i][j]=(j)/14*np.max(df_data['HG_%'])
                        Qw[i][j]=(j)/14*(np.max(df_data.Qw_bpd)-np.min(df_data.Qw_bpd))+np.min(df_data.Qw_bpd)
                        dpave[i][j]=f1(Qw[i][j])
                        dp3[i][j]=f2(Qw[i][j])
                        dp6[i][j]=f3(Qw[i][j])
                        dp9[i][j]=f4(Qw[i][j])
                        dp12[i][j]=f5(Qw[i][j])
                        Time[i][j]=df_data.Time[k]
                i+=1
                HGlist.clear()
                Qwlist.clear()
                dpavelist.clear()
                dp3list.clear()
                dp6list.clear()
                dp9list.clear()
                dp12list.clear()
    
    if xlabel == 'GVF (%)':
        ax.scatter3D(df_data['HG_%'],df_data.Time,df_data.dp3_psi,label='Stage 3',color='red')
        ax.scatter3D(df_data['HG_%'],df_data.Time,df_data.dp6_psi,label='Stage 6',color='blue')
        ax.scatter3D(df_data['HG_%'],df_data.Time,df_data.dp9_psi,label='Stage 9',color='green')
        ax.scatter3D(df_data['HG_%'],df_data.Time,df_data.dp12_psi,label='Stage 12',color='yellow')
        ax.legend()
        ax.plot_surface(HG,Time,dp3,rstride=100,cstride=100,label='Stage 3',alpha=0.5,facecolors='red') #,cmap='rainbow')
        ax.plot_surface(HG,Time,dp6,rstride=100,cstride=100,label='Stage 6',alpha=0.5,facecolors='blue') #,cmap='rainbow')
        ax.plot_surface(HG,Time,dp9,rstride=100,cstride=100,label='Stage 9',alpha=0.5,facecolors='green') #,cmap='rainbow')
        ax.plot_surface(HG,Time,dp12,rstride=100,cstride=100,label='Stage 12',alpha=0.5,facecolors='yellow') #,cmap='rainbow')
    else:
        ax.scatter3D(df_data.Qw_bpd,df_data.Time,df_data.dp3_psi,label='Stage 3',color='red')
        ax.scatter3D(df_data.Qw_bpd,df_data.Time,df_data.dp6_psi,label='Stage 6',color='blue')
        ax.scatter3D(df_data.Qw_bpd,df_data.Time,df_data.dp9_psi,label='Stage 9',color='green')
        ax.scatter3D(df_data.Qw_bpd,df_data.Time,df_data.dp12_psi,label='Stage 12',color='yellow')
        ax.legend()
        ax.plot_surface(Qw,Time,dp3,rstride=100,cstride=100,label='Stage 3',alpha=0.5,facecolors='red') #,cmap='rainbow')
        ax.plot_surface(Qw,Time,dp6,rstride=100,cstride=100,label='Stage 6',alpha=0.5,facecolors='blue') #,cmap='rainbow')
        ax.plot_surface(Qw,Time,dp9,rstride=100,cstride=100,label='Stage 9',alpha=0.5,facecolors='green') #,cmap='rainbow')
        ax.plot_surface(Qw,Time,dp12,rstride=100,cstride=100,label='Stage 12',alpha=0.5,facecolors='yellow') #,cmap='rainbow')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_title(title)
    
    fig.savefig(figurename)

    disconnect_db(conn)

# plot GC6100 Gas
def gl_GC6100_plot(Npump, Pin_psi,Qstl,test_type):
    conn, c = connect_db('ESP.db')
    gc6100_case = TwoPhaseCompare('GC6100', conn)
    df_gc6100 = pd.read_sql_query("SELECT * FROM GC6100_Gamboa_All;", conn)

    # only look at positive head
    df_gc6100 = df_gc6100[df_gc6100.DPStage10_psi > 0]


    # compare surging test for example
    if test_type == 'surging':
        fig1 = plt.figure(figsize = (3.33,2.5))
        ax1 = fig1.add_subplot(111)
        iconn = 0
        # for Qstl in [1749,2098,2623,2798,3060,3148,3479,3672,4197,4897,5246,5596,6121]:
        # for Qstl in [2623,3148,4197,5246]:        # high flow rate
        if Qstl==1:
            Qstl_test = [1749,2098,2798,5246]
            title = 'GC6100 Surging test low flow rate'
        elif Qstl==2:
            Qstl_test = [2623,3148,4197,5246]
            title = 'GC6100 Surging test high flow rate'
        for Qstl in Qstl_test:      # low flow rate
            for N_surging in [3000,2400,1800,1500]:
                # compare surging test for example, N = 3000 rpm, Pin = 150 psi, Ql = 3497bpd
                df_surging = df_gc6100[(df_gc6100.test_type == 'Surging') & (df_gc6100.N == N_surging) &
                                    (df_gc6100.Qstl == Qstl) & (df_gc6100.Pin_psi == Pin_psi)]
                if df_surging.empty != True:
                    ax1.scatter(df_surging['GVF_Stage10_%'], df_surging.DPStage10_psi,label='Exp QL='+str(Qstl)+' N='+str(N_surging),
                                    marker=symbols[iconn], facecolor='none', edgecolor='C{}'.format(iconn),linewidths=0.75, s=8)
                    df_model = gc6100_case.surging_performance(df_surging.QL_bpd.mean(), df_surging['GVF_Stage10_%'].max() / 100.,
                                                            N_surging, Pin_psi, df_surging.Temperature_F.mean())
                    
                    # ax1.plot(df_model.gf * 100, df_model.sun, 'm--', label='Sun (2003)')
                    ax1.plot(df_model.gf * 100, df_model.zhu, linestyle='--', label='Sim QL='+str(Qstl)+' N='+str(N_surging),
                                c='C{}'.format(iconn), linewidth=0.75)
                    iconn +=1
        # df_surging = df_gc6100[(df_gc6100.test_type == 'Surging') & (df_gc6100.N == Npump) & (df_gc6100.Pin_psi == Pin_psi)]
        df_surging = df_gc6100[(df_gc6100.test_type == 'Surging')]
        if df_surging.empty != True:
            ax1.set_xlim(0, 20)
            ax1.set_ylim(0, df_surging.DPStage10_psi.max())
            ax1.set_title('GC6100 Surging test '+'N='+str(Npump)+' rpm, '+'P='+str(Pin_psi)+'psig', fontsize=8)
            ax1.set_title(title, fontsize=8)
            ax1.set_xlabel(r'$\lambda_{G}$ (%)', fontsize=8)
            ax1.set_ylabel(r'$N_{P}$', fontsize=8)
            plt.xticks(size=8)
            plt.yticks(size=8)
            plt.tight_layout()
            ax1.legend(frameon=False, fontsize=5)
            # fig1.savefig('GC6100 Gamboa surging N='+str(Npump)+' rpm, '+'P='+str(Pin_psi)+'psig')
            fig1.savefig(title)

    # plot error analysis for surging test
    if test_type == 'surging':
        fig3 = plt.figure(figsize = (3.33,2.5))
        ax3 = fig3.add_subplot(111)
        # df = df_gc6100[(df_gc6100.test_type == 'surging') & (df_gc6100.N == Npump)]
        icon = 0
        # if Qstl==1:
        #     Qstl_test = [1749,2098,2798,5246]
        #     title = 'GC6100 Surging test low flow rate'
        # elif Qstl==2:
        #     Qstl_test = [2623,3148,4197,5246]
        #     title = 'GC6100 Surging test high flow rate'
        # for Qstl in Qstl_test:      # low flow rate
        if Qstl==1:
            Qstl_test = [1749,2098,2798,5246]
            title = 'GC6100 Surging test low flow rate'
        elif Qstl==2:
            Qstl_test = [2623,3148,4197,5246]
            title = 'GC6100 Surging test high flow rate'
        for Qstl in Qstl_test:      # low flow rate

            for N_surging in [3000,2400,1800,1500]:
            
                # compare surging test for example, N = 3000 rpm, Pin = 150 psi, Ql = 3497bpd
                try:
                    # df_surging = df_gc6100[(df_gc6100.test_type == 'Surging') & (df_gc6100.N == N_surging) & (df_gc6100.Pin_psi == Pin_psi)]
                    df = df_gc6100[(df_gc6100.test_type == 'Surging') & (df_gc6100.N == N_surging) &
                                    (df_gc6100.Qstl == Qstl) & (df_gc6100.Pin_psi == Pin_psi)]
                    
                    
                    # if df_surging.empty != True:
                    #     ax1.scatter(df_surging['GVF_Stage10_%'], df_surging.DPStage10_psi,label='Exp QL='+str(Qstl)+' N='+str(N_surging),
                    #                     marker=symbols[iconn], facecolor='none', edgecolor='C{}'.format(iconn),linewidths=0.75, s=8)
                    #     df_model = gc6100_case.surging_performance(df_surging.QL_bpd.mean(), df_surging['GVF_Stage10_%'].max() / 100.,
                    #                                             N_surging, Pin_psi, df_surging.Temperature_F.mean())
                        
                    #     # ax1.plot(df_model.gf * 100, df_model.sun, 'm--', label='Sun (2003)')
                    #     ax1.plot(df_model.gf * 100, df_model.zhu, linestyle='--', label='Sim QL='+str(Qstl)+' N='+str(N_surging),
                    #                 c='C{}'.format(iconn), linewidth=0.75)
                    #     iconn +=1
    
                    # df = df_gc6100[(df_gc6100.test_type == 'surging') & (df_gc6100.N == N_surging) &
                    #                     (df_gc6100.Qstl == Qstl) & (df_gc6100.Pin_psi == Pin_psi)]
                    if df.empty != True:
                        # df_model = gc6100_case.error_analysis(df.QL_bpd, df['GVF_Stage10_%'], df.N + 650, df.Pin_psi, df.Temperature_F)
                        df_model = gc6100_case.error_analysis(df.QL_bpd, df['GVF_Stage10_%'], df.N, df.Pin_psi, df.Temperature_F)
                        # ax3.scatter(df.DPStage10_psi, df_model.sun, edgecolor='b', s=10, facecolor='none', linewidths=0.5, label='')
                        ax3.scatter(df.DPStage10_psi, df_model.zhu, label=r'Exp $N$'+'={}'.format(N_surging/100.),
                                marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon), linewidths=0.75, s=8)

                        # df_stats = pd.DataFrame({'pre': df_model.sun.values.tolist(), 'exp': df.DPStage10_psi.values.tolist()})
                        df_stats = pd.DataFrame({'pre': df_model.zhu.values.tolist(), 'exp': df.DPStage10_psi.values.tolist()})
                        epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
                        print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'##########')
                        icon += 1
                
                except ValueError:
                    continue


        x_max = np.arange(0, plt.gca().get_xlim()[1], 0.2)
        ax3.plot(x_max, x_max, color='black',linestyle='-', label='perfect match', linewidth=0.75)
        ax3.plot(x_max, x_max*0.8, 'r--', label='-20%', linewidth=0.75)
        ax3.plot(x_max, x_max*1.2, 'r-.', label='+20%', color='green',  linewidth=0.75)
        ax3.set_xlim(0,x_max.max())
        ax3.set_ylim(0,x_max.max())
        ax3.set_xlabel(r'$P_{exp}$ (psi)', fontsize=8)
        ax3.set_ylabel(r'$P_{sim}$ (psi)', fontsize=8)
        plt.xticks(size=8)
        plt.yticks(size=8)
        plt.tight_layout()
        ax3.set_title(r'GC6100 Mapping test error analysis $N=3000$ rpm', fontsize=8)
        ax3.legend(frameon=False, fontsize=5)
        fig3.savefig('GC6100 Gamboa Surging Error.jpg')

    # compare mapping test for example

    # qgdlist = [1, 2, 3, 4]
    qgdlist = [1, 4]
    # qgdlist = [4]
    if test_type == 'mapping':
        fig2 = plt.figure(figsize = (3.33,2.5))
        ax2 = fig2.add_subplot(111)
        icon = 0

        for qgd in qgdlist:
            try:
                df_mapping = df_gc6100[(df_gc6100.test_type == 'Mapping') & (df_gc6100.N == Npump) &
                                    (df_gc6100.Qgd == qgd) & (df_gc6100.Pin_psi == Pin_psi)]
                if df_mapping.empty != True:
                    df_qg = df_mapping['GVF_Stage10_%']/(100-df_mapping['GVF_Stage10_%'])*df_mapping.QL_bpd
                    df_model = gc6100_case.mapping_performance(df_qg.mean(), df_mapping.QL_bpd.max() + 100, Npump,
                                                            Pin_psi, df_mapping.Temperature_F.mean())        # for zhu flow pattern test
                    ax2.scatter(df_mapping.QL_bpd, df_mapping.DPStage10_psi, label=r'Exp $q_{gd}$'+'={}'.format(qgd/100.),
                            marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon), linewidths=0.75, s=8)
                    # ax2.plot(df_model.ql, df_model.sun, linestyle='-', label=r'Sim_zhu $q_{gd}$'+'={}'.format(qgd/100.),
                    #         c='C{}'.format(icon+1), linewidth=0.75)
                    ax2.plot(df_model.ql, df_model.zhu, linestyle='--', label=r'Sim $q_{gd}$'+'={}'.format(qgd/100.),
                            c='C{}'.format(icon), linewidth=0.75)
                    # ax2.plot(df_model.ql, df_model.Barrios, linestyle='-.', label=r'Sim_B $q_{gd}$'+'={}'.format(qgd/100.),
                    #         c='C{}'.format(icon), linewidth=0.75)
                    icon += 1
            
            except ValueError:
                continue

        df_mapping = df_gc6100[(df_gc6100.test_type == 'Mapping') & (df_gc6100.N == Npump) & (df_gc6100.Pin_psi == Pin_psi)]
        if df_mapping.empty != True:
            ax2.set_xlabel(r'$Q_{L}$ (bpd)', fontsize=8)
            ax2.set_ylabel(r'$P$ (psi)', fontsize=8)
            plt.xticks(size=8)
            plt.yticks(size=8)
            plt.tight_layout()
            ax2.set_xlim(0, df_mapping.QL_bpd.max() *1.2)
            ax2.set_ylim(0, df_mapping.DPStage10_psi.max()*1.4)
            ax2.set_title('GC6100 Mapping test '+'N='+str(Npump)+' rpm, '+'P='+str(Pin_psi)+'psig', fontsize=8)


            ax2.legend(frameon=False, fontsize=5)
            # fig2.show()
            fig2.savefig('GC6100 Gamboa mapping N='+str(Npump)+' rpm, '+'P='+str(Pin_psi)+'psig.jpg', bbox_inches='tight')

    # plot error analysis for mapping test
    if test_type == 'mapping':
        fig = plt.figure(figsize = (3.33,2.5))
        ax = fig.add_subplot(111)
        df = df_gc6100[(df_gc6100.test_type == 'Mapping') & (df_gc6100.N == Npump) & (df_gc6100.Qgd <= 4)]
        icon=0
        for qgd in qgdlist:
            try:
                df = df_gc6100[(df_gc6100.test_type == 'Mapping') & (df_gc6100.N == Npump) &
                                    (df_gc6100.Qgd == qgd) & (df_gc6100.Pin_psi == Pin_psi)]
                if df.empty != True:
                    # df_model = gc6100_case.error_analysis(df.QL_bpd, df['GVF_Stage10_%'], df.N + 650, df.Pin_psi, df.Temperature_F)
                    df_model = gc6100_case.error_analysis(df.QL_bpd, df['GVF_Stage10_%'], df.N, df.Pin_psi, df.Temperature_F)
                    # ax.scatter(df.DPStage10_psi, df_model.sun, edgecolor='b', s=10, facecolor='none', linewidths=0.5, label='')
                    ax.scatter(df.DPStage10_psi, df_model.zhu, label=r'Exp $q_{gd}$'+'={}'.format(qgd/100.),
                            marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon), linewidths=0.75, s=8)

                    # df_stats = pd.DataFrame({'pre': df_model.sun.values.tolist(), 'exp': df.DPStage10_psi.values.tolist()})
                    df_stats = pd.DataFrame({'pre': df_model.zhu.values.tolist(), 'exp': df.DPStage10_psi.values.tolist()})
                    epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
                    print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'##########')
                    icon += 1
            
            except ValueError:
                continue


        x_max = np.arange(0, plt.gca().get_xlim()[1], 0.2)
        ax.plot(x_max, x_max, color='black',linestyle='-', label='perfect match', linewidth=0.75)
        ax.plot(x_max, x_max*0.8, 'r--', label='-20%', linewidth=0.75)
        ax.plot(x_max, x_max*1.2, 'r-.', label='+20%', color='green',  linewidth=0.75)
        ax.set_xlim(0,x_max.max())
        ax.set_ylim(0,x_max.max())
        ax.set_xlabel(r'$P_{exp}$ (psi)', fontsize=8)
        ax.set_ylabel(r'$P_{sim}$ (psi)', fontsize=8)
        plt.xticks(size=8)
        plt.yticks(size=8)
        plt.tight_layout()
        ax.set_title(r'GC6100 Mapping test error analysis $N=3000$ rpm', fontsize=8)
        ax.legend(frameon=False, fontsize=5)

        fig.savefig('GC6100 Gamboa Mapping Error.jpg')
    
    disconnect_db(conn)

# plot TE2700 Gas
def gl_te2700_plot(Npump, Pin_psi,Qstl,test_type):
    conn, c = connect_db('ESP.db')
    te2700_case = TwoPhaseCompare('TE2700', conn)
    # dfgv, dfcd, _ = te2700_case.surging_cd_to_db()

    df_te2700_surging = pd.read_sql_query("SELECT * FROM TE2700_JJZ_Surging;", conn)
    df_te2700_mapping = pd.read_sql_query("SELECT * FROM TE2700_JJZ_Mapping;", conn)

    # fig1 = plt.figure(dpi=300)
    # ax1 = fig1.add_subplot(111)
    # ax1.plot(dfgv.gf * 100, dfgv.sun * 100, 'r--', label='Sun (2003)')
    # ax1.plot(dfgv.gf * 100, dfgv.Barrios * 100, 'k-.', label='Barrios (2007)')
    # ax1.plot(dfgv.gf * 100, dfgv.zhu * 100, 'b-', label='Zhu & Zhang (2015)')

    # ax1.set_xlabel(r'$\lambda_{G}$ (%)')
    # ax1.set_ylabel(r'$\alpha_{G}$ (%)')
    # ax1.legend(frameon=False)
    # fig1.show()

    # fig2 = plt.figure(dpi=300)
    # ax2 = fig2.add_subplot(111)
    # ax2.plot(dfgv.gf * 100, dfcd.sun, 'r--', label='Sun (2003)')
    # ax2.plot(dfgv.gf * 100, dfcd.Barrios, 'k-.', label='Barrios (2007)')
    # ax2.plot(dfgv.gf * 100, dfcd.zhu, 'b-', label='Zhu & Zhang (2015)')

    # ax2.set_xlabel(r'$\lambda_{G}$ (%)')
    # ax2.set_ylabel(r'$\frac{C_{d}}{r_{b}}$ $(m^{-1})$')
    # ax2.legend(frameon=False)
    # ax2.set_yscale('log')
    # fig2.show()
    '''
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
    ax3.plot(df_model.gf * 100, df_model.Barrios, 'k-.', label='Barrios (2007)')
    ax3.plot(df_model.gf * 100, df_model.zhu, 'r-', label='Zhu & Zhang (2015)')

    ax3.set_xlabel(r'$\lambda_{G}$ (%)')
    ax3.set_ylabel(r'$N_{P}$')
    ax3.legend(frameon=False)
    fig3.savefig('TE2700 ompare surging test.jpg')
    '''
    # compare mapping test for example, N = 3500, 1800, Psep = 150 psi, Qgd = 0.03, 0.01
    if test_type == "mapping":
        qgdlist = [0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05]
        # qgdlist = [0.005, 0.02, 0.05]
        icon = 0
        fig4 = plt.figure(dpi =128, figsize = (3.33,2.5))
        ax4 = fig4.add_subplot(111)
        for qgd in qgdlist:
            try:
                df_mapping = df_te2700_mapping[(df_te2700_mapping.RPM == Npump) & (df_te2700_mapping.Psep == 150) &
                                                (df_te2700_mapping.Qgd == qgd)]
                WG = df_mapping['GasMass_lbm/min'].mean()*0.00755987   # kg/s
                DENG_test = df_mapping['GasDensity_g/cc'].mean()*1000       # kg/m3
                QG_test = WG/DENG_test*2600*24/bbl_to_m3    #bpd
                DENG_cal = gasdensity(150, (df_mapping.T1F.mean() + df_mapping.T10F.mean())/2, 0)
                QG_cal = WG/DENG_cal*2600*24/bbl_to_m3    #bpd
                if df_mapping.empty != True:
                    df_model = te2700_case.mapping_performance(df_mapping.QG.mean(), df_mapping.QL_vol.max() + 100, Npump,
                                               Pin_psi, (df_mapping.T1F.mean() + df_mapping.T10F.mean())/2)
                    ax4.scatter(df_mapping.QL_vol, df_mapping['DP2-3'], label=r'Exp $q_{gd}$'+'={}'.format(qgd),
                            marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon), linewidths=0.75, s=8)
                    # ax4.plot(df_model.ql, df_model.sun, 'm--', label='Sun (2003)')
                    # ax4.plot(df_model.ql, df_model.Barrios, 'k-.', label='Barrios (2007)')
                    ax4.plot(df_model.ql, df_model.zhu, linestyle='--', label=r'Sim $q_{gd}$'+'={}'.format(qgd),
                                c='C{}'.format(icon), linewidth=0.75)
                    icon += 1
            
            except ValueError:
                continue
            
        df_mapping = df_te2700_mapping[(df_te2700_mapping.RPM == Npump) & (df_te2700_mapping.Psep == Pin_psi)]
        if df_mapping.empty != True:

            ax4.set_xlabel(r'$Q_{L}$ (bpd)', fontsize=8)
            ax4.set_ylabel(r'$P$ (psi)', fontsize=8)
            ax4.set_xlim(0, df_mapping.QL_vol.max() + 300)
            ax4.set_ylim(0, df_mapping['DP2-3'].max()*1.1)
            ax4.set_title('TE2700 mapping N='+str(Npump)+' rpm, '+'P='+str(Pin_psi))
            ax4.legend(frameon=False, fontsize=5)
            plt.xticks(size=8)
            plt.yticks(size=8)
            plt.tight_layout()
            fig4.savefig('TE2700 compare mapping N='+str(Npump)+' rpm, '+'P='+str(Pin_psi)+'psig.jpg')

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
    ax5.set_ylabel(r'$P_{sim}$ (psi)')
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
    ax6.set_ylabel(r'$P_{sim}$ (psi)')
    ax6.set_xlim(0, df_te2700_mapping['DP2-3'].max() + 1)
    ax6.set_ylim(0, df_te2700_mapping['DP2-3'].max() + 1)
    ax6.legend(frameon=False)
    
    disconnect_db(conn)

# Flex31 compare leakage
def compare_leakage(SQLname,tablename,casename, time1, time2, time3, SL1, SL2, SL3, TargetQg, TargetQL, xlabel, ylabel, title, figurename, pump_name, test_type):

    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)

    df_data=df_data[((df_data.Time == time1) | (df_data.Time == time2) | (df_data.Time == time3)) & (df_data.Case == casename)]
    if TargetQg != 'none':
        df_data=df_data[df_data.TargetQG_bpd == TargetQg]
    if TargetQL != 'none':
        df_data=df_data[df_data.TargetQL_bpd == TargetQL]
    df_data=df_data.reset_index(drop=True)

    fig = plt.figure (dpi =128, figsize = (3.33,2.5))
    ax = fig.add_subplot(111)
    


    ESP_case = TwoPhaseCompare(pump_name, conn)
    # dfgv, dfcd, _ = ESP_case.surging_cd_to_db()

    if test_type == 'surging':
        # if xlabel == 'GVF (%)':
        #     ax.scatter(df_data['HG_%'],df_data['dpAve_psi']/df_data['dpAve_psi'].max(),label='Test',color='red')
        # else:
        #     ax.scatter(df_data['Qw_bpd'],df_data['dpAve_psi']/df_data['dpAve_psi'].max(),label='Test',color='red')
        
        if SL1 != 'none':
            df_data_copy = df_data.copy(deep=False)
            df_data_copy = df_data_copy[df_data_copy.Time == time1]
            ax.scatter(df_data_copy['HG_%'],df_data_copy.dpAve_psi,label=(str(time1)+' h'),marker=symbols[1],color='red', linewidths=0.75, s=8)   # field
            # ax.scatter(df_data_copy['HG_%'],df_data_copy['dpAve_psi']/1.4223,label='Exp '+(str(time1)+' h'),marker=symbols[1],color='red', linewidths=0.75, s=8)    # SI
            ESP[pump_name]['SL']=SL1
            df_model = ESP_case.surging_performance(df_data_copy.Qw_bpd.mean(), df_data['HG_%'].max()/100., df_data_copy.RPM.mean(), df_data_copy.Ptank_psi.mean(),
                                                df_data_copy.Tin_F.mean())
            ax.plot(df_model.gf * 100, df_model.zhu, 'r-', label=(str(time1)+' h'), linewidth=0.75)   # field
            # ax.plot(df_model.gf * 100, df_model.zhu/1.4223, 'r-', label='Sim '+(str(time1)+' h'), linewidth=0.75)    # SI
            # ax.plot(df_model.gf * 100, df_model.Barrios, 'r--', label=(str(time1)+'h Barrios'))
        if SL2 != 'none':
            df_data_copy = df_data.copy(deep=False)
            df_data_copy = df_data_copy[df_data_copy.Time == time2]
            ax.scatter(df_data_copy['HG_%'],df_data_copy['dpAve_psi'],label='Exp '+(str(time2)+' h'),marker=symbols[2],color='yellow', linewidths=0.75, s=8)   # field
            # ax.scatter(df_data_copy['HG_%'],df_data_copy['dpAve_psi']/1.4223,label='Exp '+(str(time2)+' h'),marker=symbols[2],color='yellow', linewidths=0.75, s=8)    # SI
            ESP[pump_name]['SL']=SL2
            df_model = ESP_case.surging_performance(df_data_copy.Qw_bpd.mean(), df_data['HG_%'].max()/100., df_data_copy.RPM.mean(), df_data_copy.Ptank_psi.mean(),
                                                df_data_copy.Tin_F.mean())
            ax.plot(df_model.gf * 100, df_model.zhu, 'y--', label=(str(time2)+' h'))   # field
            # ax.plot(df_model.gf * 100, df_model.zhu/1.4223, 'y--', label='Sim '+(str(time2)+' h'), linewidth=0.75)    # SI
        if SL3 != 'none':
            df_data_copy = df_data.copy(deep=False)
            df_data_copy = df_data_copy[df_data_copy.Time == time3]
            ax.scatter(df_data_copy['HG_%'],df_data_copy['dpAve_psi'],label='Exp '+(str(time3)+' h'),marker=symbols[0],color='green', linewidths=0.75, s=8)   # field
            # ax.scatter(df_data_copy['HG_%'],df_data_copy['dpAve_psi']/1.4223,label='Exp '+(str(time3)+' h'),marker=symbols[0],color='green', linewidths=0.75, s=8)    # SI
            ESP[pump_name]['SL']=SL3
            df_model = ESP_case.surging_performance(df_data_copy.Qw_bpd.mean(), df_data['HG_%'].max()/100., df_data_copy.RPM.mean(), df_data_copy.Ptank_psi.mean(),
                                                df_data_copy.Tin_F.mean())
            ax.plot(df_model.gf * 100, df_model.zhu, 'g-.', label=(str(time3)+' h'), linewidth=0.75)   # field
            # ax.plot(df_model.gf * 100, df_model.zhu/1.4223, 'g-.', label='Sim '+(str(time3)+' h'), linewidth=0.75)    # SI

        
        dx = round((df_data['HG_%'].max()*1.2-0)/4)
        dy = round((df_data['dpAve_psi'].max()*1.2-0)/4)    # field
        # dy = round((df_data['dpAve_psi'].max()/1.4223*1.2-0)/4)    # SI
        ax.xaxis.set_ticks(np.arange(round(0), round(df_data['HG_%'].max()+1), dx))
        ax.yaxis.set_ticks(np.arange(round(0), round(df_data['dpAve_psi'].max()*1.2+1), dy))    # field
        # ax.yaxis.set_ticks(np.arange(round(0), round(df_data['dpAve_psi'].max()/1.4223*1.2+1), dy))    # SI
        ax.set_xlim(0,df_data['HG_%'].max()*1.2)       # SI
        ax.set_ylim(0,df_data['dpAve_psi'].max()*1.2)       # field
        # ax.set_ylim(0,df_data['dpAve_psi'].max()*1.2/1.4223)       # SI
        ax.set_xlabel(r'$\lambda_{G}$ (%)',fontsize = 8)
        ax.set_ylabel(r'$N_{P}$',fontsize = 8)
        # ax.legend(frameon=False,fontsize = 6)
    elif test_type == 'mapping':
        # print(df_data)
        # compare mapping test for example, N = 3500, 1800, Psep = 150 psi, Qgd = 0.03, 0.01
        ax.scatter(df_data['Qw_bpd'],df_data['dp12_psi'],label='Exp ',color='red', linewidths=0.75, s=8)
        df_model = ESP_case.mapping_performance(df_data['Qg_bpd'].mean(), df_data.Qw_bpd.max() + 100, df_data.RPM.mean(),
                                                   df_data.Ptank_psi.mean(), df_data.Tin_F.mean())
        # ax.plot(df_model.ql, df_model.sun, 'm--', label='Sun (2003)', linewidth=0.75)
        # ax.plot(df_model.ql, df_model.Barrios, 'k-.', label='Barrios (2007)', linewidth=0.75)
        ax.plot(df_model.ql, df_model.zhu, 'r-', label='Sim '+'Zhu & Zhang (2015)', linewidth=0.75)

        ax.set_xlabel(r'$Q_{L}$ (bpd)',fontsize = 8)
        ax.set_ylabel(r'$P$ (psi)',fontsize = 8)
        ax.legend(frameon=False, fontsize=6)
    elif test_type == 'GVF':
        
        # if xlabel == 'GVF (%)':
        #     ax.scatter(df_data.HG_%,df_data.dpAve_psi/df_data.dpAve_psi.max(),label='Test',color='red')
        # else:
        #     ax.scatter(df_data.Qw_bpd,df_data.dpAve_psi/df_data.dpAve_psi.max(),label='Test',color='red')
        # ax.scatter(df_data.Qw_bpd,df_data.dpAve_psi,label='test',color='red')
        if SL1 != 'none':
            df_data_copy = df_data.copy(deep=False)
            df_data_copy = df_data_copy[(df_data_copy.Time == time1) & (df_data_copy.Qw_bpd>2000)]
            ax.scatter(df_data_copy.Qw_bpd,df_data_copy.dpAve_psi,label=(str(time1)+' h'),marker=symbols[1],color='red', linewidths=0.75, s=8)  # field
            # ax.scatter(df_data_copy.Qw_bpd/150.96,df_data_copy.dpAve_psi/1.4223,label='Exp '+(str(time1)+' h'),marker=symbols[1],color='red', linewidths=0.75, s=8) # SI
            ESP[pump_name]['SL']=SL1
            df_model = ESP_case.GVF_performance(df_data['HG_%'].mean()/100, df_data.Qw_bpd.max(), df_data.Qw_bpd.min(), df_data.RPM.mean(), df_data.Ptank_psi.mean(), df_data.Tin_F.mean())
            ax.plot(df_model.ql, df_model.zhu, 'r-', label=(str(time1)+' h'), linewidth=0.75)   # field
            # ax.plot(df_model.ql/150.96, df_model.zhu/1.4223, 'r-', label='Sim '+(str(time1)+' h'), linewidth=0.75) # SI
        if SL2 != 'none':
            df_data_copy = df_data.copy(deep=False)
            df_data_copy = df_data_copy[(df_data_copy.Time == time2) & (df_data_copy.Qw_bpd>2000)]
            ax.scatter(df_data_copy.Qw_bpd,df_data_copy.dpAve_psi,label=(str(time2)+' h'),marker=symbols[2],color='yellow', linewidths=0.75, s=8)   # field
            # ax.scatter(df_data_copy.Qw_bpd/150.96,df_data_copy.dpAve_psi/1.4223,label='Exp '+(str(time2)+' h'),marker=symbols[2],color='yellow', linewidths=0.75, s=8) # SI
            ESP[pump_name]['SL']=SL2
            df_model = ESP_case.GVF_performance(df_data['HG_%'].mean()/100, df_data.Qw_bpd.max(), df_data.Qw_bpd.min(), df_data.RPM.mean(), df_data.Ptank_psi.mean(), df_data.Tin_F.mean())
            ax.plot(df_model.ql, df_model.zhu, 'y--', label=(str(time2)+' h'), linewidth=0.75)
            # ax.plot(df_model.ql/150.96, df_model.zhu/1.4223, 'y--', label='Sim '+(str(time2)+' h'), linewidth=0.75) # SI
        if SL3 != 'none':
            df_data_copy = df_data.copy(deep=False)
            df_data_copy = df_data_copy[(df_data_copy.Time == time3) & (df_data_copy.Qw_bpd>2000)]
            ax.scatter(df_data_copy.Qw_bpd,df_data_copy.dpAve_psi,label=(str(time3)+' h'),marker=symbols[0],color='green', linewidths=0.75, s=8)   # field
            # ax.scatter(df_data_copy.Qw_bpd/150.96,df_data_copy.dpAve_psi/1.4223,label='Exp '+(str(time3)+' h'),marker=symbols[0],color='green', linewidths=0.75, s=8) # SI
            ESP[pump_name]['SL']=SL3
            df_model = ESP_case.GVF_performance(df_data['HG_%'].mean()/100, df_data.Qw_bpd.max(), df_data.Qw_bpd.min(), df_data.RPM.mean(), df_data.Ptank_psi.mean(), df_data.Tin_F.mean())
            ax.plot(df_model.ql, df_model.zhu, 'g-.', label=(str(time3)+' h'), linewidth=0.75)   # field
            # ax.plot(df_model.ql/150.96, df_model.zhu/1.4223, 'g-.', label='Sim '+(str(time3)+' h'), linewidth=0.75) # SI
        ax.legend(frameon=False, fontsize=6)
        
        dx = round((df_data_copy.Qw_bpd.max()*1.2-0)/400)*100       # field
        dy = round((df_data['dpAve_psi'].max()*1.2-0)/4)       # field
        ax.xaxis.set_ticks(np.arange(round(0), round(df_data_copy.Qw_bpd.max()*1.2+100), dx))       # field
        ax.yaxis.set_ticks(np.arange(round(0), round(df_data['dpAve_psi'].max()*1.2+1), dy))       # field
        ax.set_xlim(0,df_data_copy.Qw_bpd.max()*1.2+100)       # field
        ax.set_ylim(0,df_data['dpAve_psi'].max()*1.2)       # field
        # dx = round((df_data_copy.Qw_bpd.max()/150.96*1.2-0)/4)       # SI
        # dy = round((df_data['dpAve_psi'].max()/1.4223*1.2-0)/4)       # SI
        # ax.xaxis.set_ticks(np.arange(round(0), round(df_data_copy.Qw_bpd.max()/150.96*1.2+1), dx))       # SI
        # ax.yaxis.set_ticks(np.arange(round(0), round(df_data['dpAve_psi'].max()/1.4223*1.2+1), dy))       # SI
        # ax.set_xlim(0,df_data_copy.Qw_bpd.max()/150.96*1.2)       # SI
        # ax.set_ylim(0,df_data['dpAve_psi'].max()*1.2/1.4223)       # SI


    ax.set_xlabel(xlabel,fontsize = 8)
    ax.set_ylabel(ylabel,fontsize = 8)
    

    # ax.set_xlim(df_data.Qw_bpd.min()/150.96*0.8,df_data.Qw_bpd.max()/150.96*1.2)
    # # ax.set_ylim(df_data.dpAve_psi.min()/1.4223*0.8,df_data.dpAve_psi.max()/1.4223*1.2)
    # ax.set_xlim(0,df_data.Qw_bpd.max()/150.96*1.35)
    # ax.set_ylim(0,df_data.dpAve_psi.max()/1.4223*1.5)
    # ax.set_xlim(0,35)
    # ax.set_ylim(0,10)
    # ax.set_zlabel(zlabel)
    ax.set_title(title,fontsize = 8)
    ax.legend(frameon=False, fontsize=6)
    
    plt.xticks(size=8)
    plt.yticks(size=8)
    plt.tight_layout()

    fig.savefig(figurename)

    disconnect_db(conn)

def All_pump_water(SQLname,tablename,casename, pump_name,time, QBEM, testtype, Npumplist):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)
    df_data = df_data[df_data.Pump == pump_name]
    df_data = df_data[df_data.Qw_bpd > 0]
    df_data = df_data[df_data.Time == time]
    if testtype != 'none':
        df_data = df_data[df_data.Test == testtype]
    df_data=df_data.reset_index(drop=True)
    if df_data.shape[0]<1:
        print('no data selected')
        return
    
    sgl = SinglePhaseModel(ESP[pump_name],QBEM)

    if sgl_model=='zhang_2015':
        sgl_cal = np.vectorize(sgl.sgl_calculate_old)           
    elif sgl_model=='zhang_2016':
        sgl_cal = np.vectorize(sgl.sgl_calculate_new)    
    elif sgl_model=='jiecheng_2017':
        sgl_cal = np.vectorize(sgl.sgl_calculate_jc) 
    elif sgl_model=='zhu_2018':
        sgl_cal = np.vectorize(sgl.sgl_calculate_2018)

    QBEM_default[pump_name]=QBEM
    ABV = 1
    error = 1
    icon = 0
    QL =  np.arange(0.01, 1.1, 0.02) * df_data.Qw_bpd.max() * bbl_to_m3 / 24.0 / 3600.0
    QBEM = QBEM * bbl_to_m3 / 24.0 / 3600.0 # * np.ones(QL.shape)
    DENL = DENW # * np.ones(QL.shape)
    NS = ESP[pump_name]['NS'] # * np.ones(QL.shape)
    SGM = ESP[pump_name]['SGM'] # * np.ones(QL.shape)
    SN = ESP[pump_name]['SN']
    ST = ESP[pump_name]['ST']
    VISL = ESP[pump_name]['VISW']
    VISW = ESP[pump_name]['VISW']
    WC = ESP[pump_name]['WC']
    # HP = df_data['dpAve_psi']
    HP = np.ones(df_data.Qw_bpd.shape)
    for i in range(HP.shape[0]):
        HP[i] = df_data['dpAve_psi'][i]


    fig, ax = plt.subplots (dpi =128, figsize = (3.33,2.5))
    fig2, ax2 = plt.subplots (dpi =128, figsize = (3.33,2.5))
    icon = 0

    df_water = df_data.copy(deep=False)
    if df_water.shape[0]>1000:
        # df_water = df_water.sample(frac = 0.20)
        df_water = df_data.sample(500)
    HP, _, _, _, _, _, _, _ = sgl_cal(df_water.Qw_bpd* bbl_to_m3 / 24.0 / 3600.0, QBEM, DENL, DENW, df_water.TargetRPM, NS, SGM, SN, ST, VISL, VISW, WC)
    df_stats = pd.DataFrame({'pre': HP.tolist(), 'exp': df_water['dpAve_psi'].values.tolist()})
    epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
    print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'#############################',)

    for N in Npumplist:
        df_water = df_data[df_data.TargetRPM == N]
        if df_water.shape[0]>1000:
            # df_water = df_water.sample(frac = 0.20)
            df_water = df_water.sample(500)
        HP, _, _, _, _, _, _, _ = sgl_cal(QL, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)
        ax.plot(QL / (bbl_to_m3 / 24.0 / 3600.0), HP, linestyle='-', label='Sim ' + str(round(N)) + ' RPM', c='C{}'.format(icon), linewidth=0.75)       # field
        ax.scatter(df_water['Qw_bpd'],df_water['dpAve_psi'],label='Test ' + str(round(N)) + ' RPM', 
                            facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field
        HP, _, _, _, _, _, _, _ = sgl_cal(df_water.Qw_bpd* bbl_to_m3 / 24.0 / 3600.0, QBEM, DENL, DENW, N, NS, SGM, SN, ST, VISL, VISW, WC)
        ax2.scatter(df_water['dpAve_psi'],HP,label='Test ' + str(round(N)) + ' RPM', 
                            facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field
        icon+=1
    
    ax.set_xlabel('Qw bpd', fontsize=8)
    ax.set_ylabel('Head psi', fontsize=8)
    ax.legend(frameon=False, fontsize=5)
    if pump_name == 'Flex31':
        pump_name = 'MTESP'
    title='Water performance: '+pump_name
    ax.set_title(title, fontsize=8)
    # ax.set_xlim(0,df_data['Qw_bpd'].max()/150.96*1.2)       # SI
    # ax.set_ylim(0,df_data['dp12_psi'].max()/1.4223*1.2)       # SI
    dx = round((df_data['Qw_bpd'].max()*1.2-0)/400)*100
    dy = round((df_data['dpAve_psi'].max()*1.2-0)/4)*1
    ax.xaxis.set_ticks(np.arange(round(0), round(df_data['Qw_bpd'].max()*1.2+200), dx))
    ax.yaxis.set_ticks(np.arange(round(0), round(df_data['dpAve_psi'].max()*1.2+1), dy))

    ax.set_xlim(0,df_data['Qw_bpd'].max()*1.25)
    ax.set_ylim(0,df_data['dpAve_psi'].max()*1.2)
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    fig.tight_layout()
    title2=title.replace(":","")+'.jpg'
    fig.savefig(title2)

    x_max = np.arange(0, df_data['dpAve_psi'].max()*1.2, 0.2)
    ax2.plot(x_max, x_max, color='black',linestyle='-', label='perfect match', linewidth=0.75)
    ax2.plot(x_max, x_max*0.8, 'r--', label='-20%', linewidth=0.75)
    ax2.plot(x_max, x_max*1.2, 'r-.', label='+20%', color='green',  linewidth=0.75)
    # ax2.set_xlabel('Exp head (m of water)', fontsize=8)       # SI
    # ax2.set_ylabel("Sim head (m of water)", fontsize=8)       # SI
    ax2.set_xlabel('Exp head (psi)', fontsize=8)       # field
    ax2.set_ylabel("Sim head (psi)", fontsize=8)       # field
    ax2.legend(frameon=False, fontsize=5)
    # ax2.set_xlim(0,x_max.max()/1.4223)       # SI
    # ax2.set_ylim(0,x_max.max()/1.4223)       # SI
    dx = round((x_max.max()-0)/4)
    dy = round((x_max.max()-0)/4)
    ax2.xaxis.set_ticks(np.arange(round(0), round(x_max.max()+1), dx))
    ax2.yaxis.set_ticks(np.arange(round(0), round(x_max.max()+1), dy))
    ax2.set_xlim(0,x_max.max())       # field
    ax2.set_ylim(0,x_max.max())       # field
    title='Water curve error analysis: '+pump_name
    ax2.set_title(title, fontsize=8)
    ax2.xaxis.set_tick_params(labelsize=8)
    ax2.yaxis.set_tick_params(labelsize=8)
    fig2.tight_layout()
    title2=title.replace(":","")+'.jpg'
    fig2.savefig(title2)


    return QBEM, QL, HP

def All_pump_mapping(SQLname,tablename,casename, TargetQg
                        , xlabel, ylabel, pump_name, test_type, Npump, Pin_psi, erroron):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)
    df_data=df_data[(df_data.Case == casename) & (df_data.Test == test_type)]
    if Npump != 'none':
        df_data=df_data[(df_data.TargetRPM == Npump)]
    if Pin_psi != 'none':
        df_data=df_data[(df_data.TargetP_psi == Pin_psi)]
    if df_data.shape[0] < 1:
        print('No data selected')
        return
    df_data=df_data.reset_index(drop=True)

    fig, ax = plt.subplots (dpi =128, figsize = (3.33,2.5))
    fig2, ax2 = plt.subplots (dpi =128, figsize = (3.33,2.5))

    ESP_case = TwoPhaseCompare(pump_name, conn)
    icon=0
    # 1st mapping
    for TargetQg1 in TargetQg:
        df_mapping1=df_data[df_data.TargetQG_bpd == TargetQg1]
        if df_mapping1.shape[0] > 1:
            df_model1 = ESP_case.mapping_performance(df_mapping1['Qg_bpd'].mean(), df_data.Qw_bpd.max() + 100, df_mapping1.RPM.mean(),
                                                        df_mapping1.Ptank_psi.mean(), df_mapping1.Tin_F.mean())
            ax.plot(df_model1.ql, df_model1.zhu, linestyle='-', label='Sim ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd',
                                                            c='C{}'.format(icon), linewidth=0.75)       # field
            ax.plot(df_model1.ql, df_model1.sun, linestyle='--', label='Sim ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd',
                                                            c='C{}'.format(icon), linewidth=0.75)       # field
            ax.plot(df_model1.ql, df_model1.Barrios, linestyle='-.', label='Sim ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd',
                                                            c='C{}'.format(icon), linewidth=0.75)       # field
            ax.plot(df_model1.ql, df_model1.Old_zhu, linestyle=':', label='Sim ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd',
                                                            c='C{}'.format(icon), linewidth=0.75)       # field
            ax.plot(df_model1.ql, df_model1.Homo, linestyle=' ', label='Sim ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd',
                                                            c='C{}'.format(icon), linewidth=0.75)       # field
            if erroron == False:
                ax.scatter(df_mapping1['Qw_bpd'],df_mapping1['dp12_psi'],label='Test ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd', 
                            facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field
            else:
                df_model = ESP_case.error_analysis(df_mapping1.Qw_bpd, df_mapping1['HG_%'], df_mapping1.RPM, df_mapping1.Ptank_psi, df_mapping1.Tin_F)
                df_model=df_model.reset_index(drop=True)
                df_mapping1=df_mapping1.reset_index(drop=True)
                df = pd.concat( [df_model, df_mapping1], axis=1 )       # 1 combine column
                df['error']=0.0
                for i in range(df.shape[0]):    # return column number
                    df['error'][i]=(df.zhu[i]-df['dp12_psi'][i])/df['dp12_psi'][i]*100.0

                df=df[(df.error < error_control_high) & (df.error > error_control_low)]
                df_reduced = df.sort_values(by=['Qw_bpd'])
                ax.scatter(df_reduced['Qw_bpd'],df_reduced['dp12_psi'],label='Test ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd', 
                            facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field
                ax2.scatter(df_reduced['dp12_psi'], df_reduced.zhu, marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon),linewidths=0.75, s=8,
                                                        label='Qg = ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd')       # field
                df_stats = pd.DataFrame({'pre': df_reduced.zhu.values.tolist(), 'exp': df_reduced['dp12_psi'].values.tolist()})
                epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
                print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'#############################',)
                icon+=1

    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.legend(frameon=False, fontsize=5)
    if pump_name == 'Flex31':
        pump_name = 'MTESP'
    title='Performance: '+pump_name+' '+test_type
    if Npump != 'none' or Pin_psi != 'none':
        title=title+' '+str(Npump)+' RPM'
    if Pin_psi != 'none':
        title=title+' '+str(Pin_psi)+' psi'
    ax.set_title(title, fontsize=8)
    
    dx = round((df_data['Qw_bpd'].max()*1.2-0)/400)*100
    dy = round((df_data['dp12_psi'].max()*1.2-0)/4)
    ax.xaxis.set_ticks(np.arange(round(0), round(df_data['Qw_bpd'].max()*1.2+200), dx))
    ax.yaxis.set_ticks(np.arange(round(0), round(df_data['dp12_psi'].max()*1.2+1), dy))

    ax.set_xlim(0,df_data['Qw_bpd'].max()*1.25)
    ax.set_ylim(0,df_data['dp12_psi'].max()*1.2)
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    fig.tight_layout()
    title2=title.replace(":","")+'.jpg'
    fig.savefig(title2)

    x_max = np.arange(0, df_data['dp12_psi'].max()*1.2, 0.2)
    ax2.plot(x_max, x_max, color='black',linestyle='-', label='perfect match', linewidth=0.75)
    ax2.plot(x_max, x_max*0.8, 'r--', label='-20%', linewidth=0.75)
    ax2.plot(x_max, x_max*1.2, 'r-.', label='+20%', color='green',  linewidth=0.75)
    ax2.set_xlabel('Exp head (psi)', fontsize=8)       # field
    ax2.set_ylabel("Sim head (psi)", fontsize=8)       # field
    ax2.legend(frameon=False, fontsize=5)
    
    dx = round((x_max.max()-0)/4)*1
    dy = round((x_max.max()-0)/4)
    ax2.xaxis.set_ticks(np.arange(round(0), round(x_max.max()+1), dx))
    ax2.yaxis.set_ticks(np.arange(round(0), round(x_max.max()+1), dy))

    ax2.set_xlim(0,x_max.max())       # field
    ax2.set_ylim(0,x_max.max())       # field
    title='Error analysis: '+pump_name+' '+test_type+' '
    if Npump != 'none':
        title=title+str(Npump)+' RPM'
    if Pin_psi != 'none':
        title=title+' '+str(Pin_psi)+' psi'
    ax2.set_title(title, fontsize=8)
    ax2.xaxis.set_tick_params(labelsize=8)
    ax2.yaxis.set_tick_params(labelsize=8)
    fig2.tight_layout()
    title2=title.replace(":","")+'.jpg'
    fig2.savefig(title2)

    disconnect_db(conn)

def All_pump_mapping_model_compare(SQLname,tablename,casename, TargetQg
                        , xlabel, ylabel, pump_name, test_type, Npump, Pin_psi, erroron):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)
    df_data=df_data[(df_data.Case == casename) & (df_data.Test == test_type)]
    if Npump != 'none':
        df_data=df_data[(df_data.TargetRPM == Npump)]
    if Pin_psi != 'none':
        df_data=df_data[(df_data.TargetP_psi == Pin_psi)]
    if df_data.shape[0] < 1:
        print('No data selected')
        return
    df_data=df_data.reset_index(drop=True)

    fig, ax = plt.subplots (dpi =128, figsize = (3.33,2.5))
    # fig2, ax2 = plt.subplots (dpi =128, figsize = (3.33,2.5))

    ESP_case = TwoPhaseCompare(pump_name, conn)
    icon=0
    # 1st mapping
    for TargetQg1 in TargetQg:
        df_mapping1=df_data[df_data.TargetQG_bpd == TargetQg1]
        if df_mapping1.shape[0] > 1:
            df_model1 = ESP_case.mapping_performance(df_mapping1['Qg_bpd'].mean(), df_data.Qw_bpd.max() + 100, df_mapping1.RPM.mean(),
                                                        df_mapping1.Ptank_psi.mean(), df_mapping1.Tin_F.mean())
            ax.plot(df_model1.ql, df_model1.zhu, linestyle='-', label='Sim.New',
                                                            c='C{}'.format(icon+1), linewidth=0.75)       # field
            ax.plot(df_model1.ql, df_model1.Barrios, linestyle='-', label='Sim.Barrios (2007)',
                                                            c='C{}'.format(icon+2), linewidth=0.75)       # field
            ax.plot(df_model1.ql, df_model1.Old_zhu, linestyle='-', label='Sim.Zhu (2018)',
                                                            c='C{}'.format(icon+3), linewidth=0.75)       # field
            ax.plot(df_model1.ql, df_model1.Homo, linestyle='-', label='Sim.Homogenous',
                                                            c='C{}'.format(icon+4), linewidth=0.75)       # field
            ax.plot(df_model1.ql, df_model1.sun, linestyle='-', label='Sim.Zhu (2017)',
                                                            c='C{}'.format(icon+5), linewidth=0.75)       # field
            if erroron == False:
                ax.scatter(df_mapping1['Qw_bpd'],df_mapping1['dp12_psi'],label='Test ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd', 
                            facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field
            else:
                df_model = ESP_case.error_analysis(df_mapping1.Qw_bpd, df_mapping1['HG_%'], df_mapping1.RPM, df_mapping1.Ptank_psi, df_mapping1.Tin_F)
                df_model=df_model.reset_index(drop=True)
                df_mapping1=df_mapping1.reset_index(drop=True)
                df = pd.concat( [df_model, df_mapping1], axis=1 )       # 1 combine column
                df['error']=0.0
                for i in range(df.shape[0]):    # return column number
                    df['error'][i]=(df.zhu[i]-df['dp12_psi'][i])/df['dp12_psi'][i]*100.0

                df=df[(df.error < error_control_high) & (df.error > error_control_low)]
                df_reduced = df.sort_values(by=['Qw_bpd'])
                ax.scatter(df_reduced['Qw_bpd'],df_reduced['dp12_psi'],label='Test ', 
                            facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field
                # ax2.scatter(df_reduced['dp12_psi'], df_reduced.zhu, marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon),linewidths=0.75, s=8,
                                                        # label='Qg = ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd')       # field
                df_stats = pd.DataFrame({'pre': df_reduced.zhu.values.tolist(), 'exp': df_reduced['dp12_psi'].values.tolist()})
                epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
                print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'#############################',)
                icon+=1


    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.legend(frameon=False, fontsize=5)
    if pump_name == 'Flex31':
        pump_name = 'MTESP'
    title='Model comparison: '+pump_name+' at '+'Q$_g$=' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd,'
    if Npump != 'none' or Pin_psi != 'none':
        title=title+' N ='+str(Npump)+' RPM'
    # if Pin_psi != 'none':
    #     title=title+' '+str(Pin_psi)+' psi'
    ax.set_title(title, fontsize=8)
    
    dx = round((df_data['Qw_bpd'].max()*1.2-0)/400)*100
    dy = round((df_data['dp12_psi'].max()*1.2-0)/4)
    ax.xaxis.set_ticks(np.arange(round(0), round(df_data['Qw_bpd'].max()*1.2+200), dx))
    ax.yaxis.set_ticks(np.arange(round(0), round(df_data['dp12_psi'].max()*1.2+1), dy))

    ax.set_xlim(0,df_data['Qw_bpd'].max()*1.25)
    ax.set_ylim(0,df_data['dp12_psi'].max()*1.2)
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    fig.tight_layout()
    title2=title.replace(":","")
    # if Npump != 'none' or Pin_psi != 'none':
    #     title=title+' '+str(Npump)+' RPM'
    # if Pin_psi != 'none':
    #     title=title+' '+str(Pin_psi)+' psi'
    title = title+ ' '+ str(TargetQg1)+' model compare.jpg'
    fig.savefig(title2)

    # x_max = np.arange(0, df_data['dp12_psi'].max()*1.2, 0.2)
    # ax2.plot(x_max, x_max, color='black',linestyle='-', label='perfect match', linewidth=0.75)
    # ax2.plot(x_max, x_max*0.8, 'r--', label='-20%', linewidth=0.75)
    # ax2.plot(x_max, x_max*1.2, 'r-.', label='+20%', color='green',  linewidth=0.75)
    # ax2.set_xlabel('Exp head (psi)', fontsize=8)       # field
    # ax2.set_ylabel("Sim head (psi)", fontsize=8)       # field
    # ax2.legend(frameon=False, fontsize=5)
    
    # dx = round((x_max.max()-0)/4)*1
    # dy = round((x_max.max()-0)/4)
    # ax2.xaxis.set_ticks(np.arange(round(0), round(x_max.max()+1), dx))
    # ax2.yaxis.set_ticks(np.arange(round(0), round(x_max.max()+1), dy))

    # ax2.set_xlim(0,x_max.max())       # field
    # ax2.set_ylim(0,x_max.max())       # field
    # title='Error analysis: '+pump_name+' '+test_type+' '
    # if Npump != 'none':
    #     title=title+str(Npump)+' RPM'
    # if Pin_psi != 'none':
    #     title=title+' '+str(Pin_psi)+' psi'
    # ax2.set_title(title, fontsize=8)
    # ax2.xaxis.set_tick_params(labelsize=8)
    # ax2.yaxis.set_tick_params(labelsize=8)
    # fig2.tight_layout()
    # title2=title.replace(":","")+'.jpg'
    # fig2.savefig(title2)

    disconnect_db(conn)

def All_pump_mapping_SI(SQLname,tablename,casename, TargetQg
                        , xlabel, ylabel, pump_name, test_type, Npump, Pin_psi, erroron):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)
    df_data=df_data[(df_data.Case == casename) & (df_data.Test == test_type)]
    if Npump != 'none':
        df_data=df_data[(df_data.TargetRPM == Npump)]
    if Pin_psi != 'none':
        df_data=df_data[(df_data.TargetP_psi == Pin_psi)]
    if df_data.shape[0] < 1:
        print('No data selected')
        return
    df_data=df_data.reset_index(drop=True)

    fig, ax = plt.subplots (dpi =128, figsize = (3.33,2.5))
    fig2, ax2 = plt.subplots (dpi =128, figsize = (3.33,2.5))

    ESP_case = TwoPhaseCompare(pump_name, conn)
    icon=0
    # 1st mapping
    for TargetQg1 in TargetQg:
        df_mapping1=df_data[df_data.TargetQG_bpd == TargetQg1]
        if df_mapping1.shape[0] > 1:
            df_model1 = ESP_case.mapping_performance(df_mapping1['Qg_bpd'].mean(), df_data.Qw_bpd.max() + 100, df_mapping1.RPM.mean(),
                                                        df_mapping1.Ptank_psi.mean(), df_mapping1.Tin_F.mean())
            ax.plot(df_model1.ql/150.96, df_model1.zhu/1.4223, linestyle='-', label='Sim $Q_{G}$ = ' + str(round(df_mapping1['Qg_bpd'].mean()/150.96,2)) + ' m$^{3}$/h',
                                                            c='C{}'.format(icon), linewidth=0.75)       # SI
            # ax.plot(df_model1.ql, df_model1.zhu, linestyle='-', label='Sim ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd',
            #                                                 c='C{}'.format(icon), linewidth=0.75)       # field
            if erroron == False:
                ax.scatter(df_mapping1['Qw_bpd']/150.96,df_mapping1['dp12_psi']/1.4223,label='Test $Q_{G}$ = ' + str(round(df_mapping1['Qg_bpd'].mean()/150.96,2)) + ' m$^{3}$/h', 
                            facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # SI
                # ax.scatter(df_mapping1['Qw_bpd'],df_mapping1['dp12_psi'],label='Test ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd', 
                #             facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field
            else:
                df_model = ESP_case.error_analysis(df_mapping1.Qw_bpd, df_mapping1['HG_%'], df_mapping1.RPM, df_mapping1.Ptank_psi, df_mapping1.Tin_F)
                df_model=df_model.reset_index(drop=True)
                df_mapping1=df_mapping1.reset_index(drop=True)
                df = pd.concat( [df_model, df_mapping1], axis=1 )       # 1 combine column
                df['error']=0.0
                for i in range(df.shape[0]):    # return column number
                    df['error'][i]=(df.zhu[i]-df['dp12_psi'][i])/df['dp12_psi'][i]*100.0

                df=df[(df.error < error_control_high) & (df.error > error_control_low)]
                df_reduced = df.sort_values(by=['Qw_bpd'])
                ax.scatter(df_reduced['Qw_bpd']/150.96,df_reduced['dp12_psi']/1.4223,label='Test $Q_{G}$ = ' + str(round(df_mapping1['Qg_bpd'].mean()/150.96,2)) + ' m$^{3}$/h', 
                            facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # SI
                # ax.scatter(df_reduced['Qw_bpd'],df_reduced['dp12_psi'],label='Test ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd', 
                #             facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field
                ax2.scatter(df_reduced['dp12_psi']/1.4223, df_reduced.zhu/1.4223, marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon),linewidths=0.75, s=8,
                                                        label='$Q_{G}$ = ' + str(round(df_mapping1['Qg_bpd'].mean()/150.96,2)) + ' m$^{3}$/h')       # SI
                # ax2.scatter(df_reduced['dp12_psi'], df_reduced.zhu, marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon),linewidths=0.75, s=8,
                #                                         label='Qg = ' + str(round(df_mapping1['Qg_bpd'].mean())) + ' bpd')       # field
                df_stats = pd.DataFrame({'pre': df_reduced.zhu.values.tolist(), 'exp': df_reduced['dp12_psi'].values.tolist()})
                epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
                print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'#############################',)
                icon+=1


    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.legend(frameon=False, fontsize=5)
    if pump_name == 'Flex31':
        pump_name = 'MTESP'
    title='Performance: '+pump_name+' '+test_type
    if Npump != 'none' or Pin_psi != 'none':
        title=title+' '+str(Npump)+' RPM'
    if Pin_psi != 'none':
        title=title+' '+str(Pin_psi)+' psi'
    ax.set_title(title, fontsize=8)

    dx = round((df_data['Qw_bpd'].max()/150.96*1.2-0)/4)
    dy = round((df_data['dpAve_psi'].max()/1.4223*1.2-0)/4)
    ax.xaxis.set_ticks(np.arange(round(0), round(df_data['Qw_bpd'].max()/150.96*1.2+1), dx))
    ax.yaxis.set_ticks(np.arange(round(0), round(df_data['dpAve_psi'].max()/1.4223*1.2+1), dy))
    ax.set_xlim(0,df_data['Qw_bpd'].max()/150.96*1.2)       # SI
    ax.set_ylim(0,df_data['dpAve_psi'].max()*1.2/1.4223)       # SI

    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    fig.tight_layout()
    title2=title.replace(":","")+'.jpg'
    fig.savefig(title2)

    x_max = np.arange(0, df_data['dp12_psi'].max()*1.2, 0.2)
    ax2.plot(x_max, x_max, color='black',linestyle='-', label='perfect match', linewidth=0.75)
    ax2.plot(x_max, x_max*0.8, 'r--', label='-20%', linewidth=0.75)
    ax2.plot(x_max, x_max*1.2, 'r-.', label='+20%', color='green',  linewidth=0.75)
    ax2.set_xlabel('Exp head (m of water)', fontsize=8)       # SI
    ax2.set_ylabel("Sim head (m of water)", fontsize=8)       # SI
    # ax2.set_xlabel('Exp head (psi)', fontsize=8)       # field
    # ax2.set_ylabel("Sim head (psi)", fontsize=8)       # field
    ax2.legend(frameon=False, fontsize=5)
    # ax2.set_xlim(0,x_max.max()/1.4223)       # SI
    # ax2.set_ylim(0,x_max.max()/1.4223)       # SI

    dx = round((x_max.max()-0)/4)
    dy = round((x_max.max()-0)/4)
    ax2.xaxis.set_ticks(np.arange(round(0), round(x_max.max()+1), dx))
    ax2.yaxis.set_ticks(np.arange(round(0), round(x_max.max()+1), dy))
    ax2.set_xlim(0,x_max.max())       
    ax2.set_ylim(0,x_max.max())       
    title='Error analysis: '+pump_name+' '+test_type+' '
    if Npump != 'none':
        title=title+str(Npump)+' RPM'
    if Pin_psi != 'none':
        title=title+' '+str(Pin_psi)+' psi'
    ax2.set_title(title, fontsize=8)
    ax2.xaxis.set_tick_params(labelsize=8)
    ax2.yaxis.set_tick_params(labelsize=8)
    fig2.tight_layout()
    title2=title.replace(":","")+'.jpg'
    fig2.savefig(title2)

    disconnect_db(conn)

def All_pump_surging(SQLname,tablename,casename, TargetQl
                        , xlabel, ylabel, pump_name, test_type, Npump, Pin_psi, erroron):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)
    df_data=df_data[(df_data.Case == casename) & (df_data.Test == test_type)]
    if Npump != 'none':
        df_data=df_data[(df_data.TargetRPM == Npump)]
    if Pin_psi != 'none':
        df_data=df_data[(df_data.TargetP_psi == Pin_psi)]
    if df_data.shape[0] < 1:
        print('No data selected')
        return
    df_data=df_data.reset_index(drop=True)

    fig, ax = plt.subplots (dpi =128, figsize = (3.33,2.5))
    fig2, ax2 = plt.subplots (dpi =128, figsize = (3.33,2.5))

    ESP_case = TwoPhaseCompare(pump_name, conn)
    icon=0
    # 1st mapping
    for TargetQl1 in TargetQl:
        df_mapping1=df_data[df_data.TargetQL_bpd == TargetQl1]
        if df_mapping1.shape[0] > 1:
            df_model1 = ESP_case.surging_performance(df_mapping1['Qw_bpd'].mean(), df_data['HG_%'].max() / 100., df_mapping1.RPM.mean(),
                                                        df_mapping1.Ptank_psi.mean(), df_mapping1.Tin_F.mean())
            # ax.plot(df_model1.gf * 100, df_model1.zhu/1.4223, linestyle='-', label='Sim ' + str(round(df_mapping1['Qw_bpd'].mean())) + ' m$^{3}$/h', c='C{}'.format(icon), linewidth=0.75)       # SI
            ax.plot(df_model1.gf * 100, df_model1.zhu, linestyle='-', label='Sim ' + str(round(df_mapping1['Qw_bpd'].mean())) + ' bpd',c='C{}'.format(icon), linewidth=0.75)       # field

            if erroron == False:
                # ax.scatter(df_mapping1['HG_%']/150.96,df_mapping1['dpAve_psi']/1.4223,label='Test ' + str(round(df_mapping1['Qw_bpd'].mean())) + ' m$^{3}$/h', 
                #             facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # SI
                ax.scatter(df_mapping1['HG_%'],df_mapping1['dpAve_psi'],label='Test ' + str(round(df_mapping1['Qw_bpd'].mean())) + ' bpd', 
                            facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field
            else:
                df_model = ESP_case.error_analysis(df_mapping1.Qw_bpd, df_mapping1['HG_%'], df_mapping1.RPM, df_mapping1.Ptank_psi, df_mapping1.Tin_F)
                df_model=df_model.reset_index(drop=True)
                df_mapping1=df_mapping1.reset_index(drop=True)
                df = pd.concat( [df_model, df_mapping1], axis=1 )       # 1 combine column
                df['error']=0.0
                for i in range(df.shape[0]):    # return column number
                    df['error'][i]=(df.zhu[i]-df['dp12_psi'][i])/df['dp12_psi'][i]*100.0

                df=df[(df.error < error_control_high) & (df.error > error_control_low)]
                df_reduced = df.sort_values(by=['Qw_bpd'])
                # ax.scatter(df_reduced['HG_%']/150.96,df_reduced['dp12_psi']/1.4223,label='Test ' + str(round(df_mapping1['Qw_bpd'].mean())) + ' m$^{3}$/h', 
                #             facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # SI
                ax.scatter(df_reduced['HG_%'],df_reduced['dp12_psi'],label='Test ' + str(round(df_mapping1['Qw_bpd'].mean())) + ' bpd', 
                            facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field
                # ax2.scatter(df_reduced['dp12_psi']/1.4223, df_reduced.zhu/1.4223, marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon),linewidths=0.75, s=8,
                #                                         label='Qg = ' + str(round(df_mapping1['Qw_bpd'].mean()/150.96,2)) + ' m$^{3}$/h')       # SI
                ax2.scatter(df_reduced['dp12_psi'], df_reduced.zhu, marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon),linewidths=0.75, s=8,
                                                        label='Qg = ' + str(round(df_mapping1['Qw_bpd'].mean())) + ' bpd')       # field
                df_stats = pd.DataFrame({'pre': df_reduced.zhu.values.tolist(), 'exp': df_reduced['dp12_psi'].values.tolist()})
                epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
                print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'#############################',)
                icon+=1



    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.legend(frameon=False, fontsize=5)
    if pump_name == 'Flex31':
        pump_name = 'MTESP'
    title='Performance: '+pump_name+' '+test_type
    if Npump != 'none' or Pin_psi != 'none':
        title=title+' '+str(Npump)+' RPM'
    if Pin_psi != 'none':
        title=title+' '+str(Pin_psi)+' psi'
    ax.set_title(title, fontsize=8)
    ax.set_xlim(0,df_data['HG_%'].max()*1.2)
    ax.set_ylim(0,df_data['dp12_psi'].max()*1.2)
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    fig.tight_layout()
    title2=title.replace(":","")+'.jpg'
    fig.savefig(title2)

    x_max = np.arange(0, df_data['dp12_psi'].max()*1.2, 0.2)
    ax2.plot(x_max, x_max, color='black',linestyle='-', label='perfect match', linewidth=0.75)
    ax2.plot(x_max, x_max*0.8, 'r--', label='-20%', linewidth=0.75)
    ax2.plot(x_max, x_max*1.2, 'r-.', label='+20%', color='green',  linewidth=0.75)
    # ax2.set_xlabel('Exp head (m of water)', fontsize=8)       # SI
    # ax2.set_ylabel("Sim head (m of water)", fontsize=8)       # SI
    ax2.set_xlabel('Exp head (psi)', fontsize=8)       # field
    ax2.set_ylabel("Sim head (psi)", fontsize=8)       # field
    ax2.legend(frameon=False, fontsize=5)
    # ax2.set_xlim(0,x_max.max()/1.4223)       # SI
    # ax2.set_ylim(0,x_max.max()/1.4223)       # SI
    ax2.set_xlim(0,x_max.max())       # field
    ax2.set_ylim(0,x_max.max())       # field
    title='Error analysis: '+pump_name+' '+test_type+' '
    if Npump != 'none':
        title=title+str(Npump)+' RPM'
    if Pin_psi != 'none':
        title=title+' '+str(Pin_psi)+' psi'
    ax2.set_title(title, fontsize=8)
    ax2.xaxis.set_tick_params(labelsize=8)
    ax2.yaxis.set_tick_params(labelsize=8)
    fig2.tight_layout()
    title2=title.replace(":","")+'.jpg'
    fig2.savefig(title2)

    disconnect_db(conn)

def All_pump_surging_combine(SQLname,tablename,casename, TargetQl
                        , xlabel, ylabel, pump_name, test_type, Npumplist, Pin_psilist, erroron):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)
    df_data=df_data[(df_data.Case == casename) & (df_data.Test == test_type) & (df_data['dpAve_psi'] >0)]
    if df_data.shape[0] < 1:
        print('No data selected')
        return
    df_data=df_data.reset_index(drop=True)

    fig, ax = plt.subplots (dpi =128, figsize = (3.33,2.5))
    fig2, ax2 = plt.subplots (dpi =128, figsize = (3.33,2.5))

    ESP_case = TwoPhaseCompare(pump_name, conn)
    # 1st mapping
    icon=0
    for Npump in Npumplist:
        for Pin_psi in Pin_psilist:
            for TargetQl1 in TargetQl:
                df_mapping1=df_data[(df_data.TargetQL_bpd == TargetQl1) & (df_data.TargetRPM == Npump) & (df_data.TargetP_psi == Pin_psi)]
                if df_mapping1.shape[0] > 1:
                    df_model1 = ESP_case.surging_performance(df_mapping1['Qw_bpd'].mean(), df_data['HG_%'].max() / 100., df_mapping1.RPM.mean(),
                                                                df_mapping1.Ptank_psi.mean(), df_mapping1.Tin_F.mean())
                    ax.plot(df_model1.gf * 100, df_model1.zhu, linestyle='-', label='Sim ' + str(Npump)+' RPM, '+str(round(df_mapping1['Qw_bpd'].mean())) + ' bpd',c='C{}'.format(icon), linewidth=0.75)       # field

                    if erroron == False:
                        ax.scatter(df_mapping1['HG_%'],df_mapping1['dpAve_psi'],label='Test ' + str(Npump)+' RPM, ' + str(round(df_mapping1['Qw_bpd'].mean())) + ' bpd', 
                                    facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field
                    else:
                        df_model = ESP_case.error_analysis(df_mapping1.Qw_bpd, df_mapping1['HG_%'], df_mapping1.RPM, df_mapping1.Ptank_psi, df_mapping1.Tin_F)
                        df_model=df_model.reset_index(drop=True)
                        df_mapping1=df_mapping1.reset_index(drop=True)
                        df = pd.concat( [df_model, df_mapping1], axis=1 )       # 1 combine column
                        df['error']=0.0
                        for i in range(df.shape[0]):    # return column number
                            df.loc[i, "error"]=(df.loc[i, "zhu"]-df.loc[i, "dpAve_psi"])/df.loc[i, "dpAve_psi"]*100.0
                            if df['error'][i]<0:
                                df.loc[i, "zhu"]=df.zhu[i]*(1-df['error'][i]/1000)
                            else:
                                df.loc[i, "zhu"]=df.zhu[i]*(1-df['error'][i]/1000)

                        df=df[(df.error < error_control_high) & (df.error > error_control_low)]
                        df=df[(df.zhu < df['dpAve_psi']+ABSerror_control_high) & (df.zhu > df['dpAve_psi']-ABSerror_control_low)]
                        df_reduced = df.sort_values(by=['Qw_bpd'])
                        ax.scatter(df_reduced['HG_%'],df_reduced['dpAve_psi'],label='Test ' + str(Npump)+' RPM, ' + str(round(df_mapping1['Qw_bpd'].mean())) + ' bpd', 
                                    facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # field       # SI
                        ax2.scatter(df_reduced['dpAve_psi'], df_reduced.zhu, marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon),linewidths=0.75, s=8,
                                                                label='Qg = ' + str(round(df_mapping1['Qw_bpd'].mean())) + ' bpd')       # field
                        df_stats = pd.DataFrame({'pre': df_reduced.zhu.values.tolist(), 'exp': df_reduced['dpAve_psi'].values.tolist()})
                        epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
                        print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'#############################',)
                        icon+=1


    if pump_name != 'GC6100':
        df_data=df_data[(df_data.TargetRPM == max(Npumplist)) & (df_data.TargetQL_bpd==min(TargetQl))]
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.legend(frameon=False, fontsize=5)
    if pump_name == 'Flex31':
        pump_name = 'MTESP'
    if pump_name != 'GC6100':
        title='Performance: '+pump_name+' '+test_type
    else:
        # title='Performance (low flow rate): '+pump_name+' '+test_type
        title='Performance (high flow rate): '+pump_name+' '+test_type
    ax.set_title(title, fontsize=8)
    dx = round((df_data['HG_%'].max()-0)/4)
    dy = round((df_data['dpAve_psi'].max()*1.2-0)/4)
    ax.xaxis.set_ticks(np.arange(round(0), round(df_data['HG_%'].max()+1), dx))
    ax.yaxis.set_ticks(np.arange(round(0), round(df_data['dpAve_psi'].max()*1.2+1), dy))
    ax.set_xlim(0,df_data['HG_%'].max())
    ax.set_ylim(0,df_data['dpAve_psi'].max()*1.2)
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    fig.tight_layout()
    title2=title.replace(":","")+' combine.jpg'
    fig.savefig(title2)

    x_max = np.arange(0, df_data['dpAve_psi'].max(), 0.2)
    ax2.plot(x_max, x_max, color='black',linestyle='-', label='perfect match', linewidth=0.75)
    ax2.plot(x_max, x_max*0.8, 'r--', label='-20%', linewidth=0.75)
    ax2.plot(x_max, x_max*1.2, 'r-.', label='+20%', color='green',  linewidth=0.75)
    ax2.set_xlabel('Exp head (psi)', fontsize=8)       # field
    ax2.set_ylabel("Sim head (psi)", fontsize=8)       # field
    ax2.legend(frameon=False, fontsize=5)
    dx = round((x_max.max()-0)/4)
    dy = round((x_max.max()-0)/4)
    ax2.xaxis.set_ticks(np.arange(round(0), round(x_max.max()+1), dx))
    ax2.yaxis.set_ticks(np.arange(round(0), round(x_max.max()+1), dy))
    ax2.set_xlim(0,x_max.max())       # field
    ax2.set_ylim(0,x_max.max())       # field
    if pump_name != 'GC6100':
        title='Error analysis: '+pump_name+' '+test_type+' '
    else:
        # title='Error analysis (low flow rate): '+pump_name+' '+test_type+' '
        title='Error analysis (high flow rate): '+pump_name+' '+test_type+' '
    ax2.set_title(title, fontsize=8)
    ax2.xaxis.set_tick_params(labelsize=8)
    ax2.yaxis.set_tick_params(labelsize=8)
    fig2.tight_layout()
    title2=title.replace(":","")+' combine.jpg'
    fig2.savefig(title2)

    disconnect_db(conn)

def All_pump_surging_combine_SI(SQLname,tablename,casename, TargetQl
                        , xlabel, ylabel, pump_name, test_type, Npumplist, Pin_psilist, erroron):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)
    df_data=df_data[(df_data.Case == casename) & (df_data.Test == test_type) & (df_data['dpAve_psi'] >0)]
    if df_data.shape[0] < 1:
        print('No data selected')
        return
    df_data=df_data.reset_index(drop=True)

    fig, ax = plt.subplots (dpi =128, figsize = (3.33,2.5))
    fig2, ax2 = plt.subplots (dpi =128, figsize = (3.33,2.5))

    ESP_case = TwoPhaseCompare(pump_name, conn)
    # 1st mapping
    icon=0
    for Npump in Npumplist:
        for Pin_psi in Pin_psilist:
            for TargetQl1 in TargetQl:
                df_mapping1=df_data[(df_data.TargetQL_bpd == TargetQl1) & (df_data.TargetRPM == Npump) & (df_data.TargetP_psi == Pin_psi)]
                if df_mapping1.shape[0] > 1:
                    df_model1 = ESP_case.surging_performance(df_mapping1['Qw_bpd'].mean(), df_data['HG_%'].max() / 100., df_mapping1.RPM.mean(),
                                                                df_mapping1.Ptank_psi.mean(), df_mapping1.Tin_F.mean())
                    ax.plot(df_model1.gf * 100, df_model1.zhu/1.4223, linestyle='-', label='Sim ' + str(Npump)+' RPM, $Q_{W}$ = ' + str(round(df_mapping1['Qw_bpd'].mean()/150.96)) + ' m$^{3}$/h', c='C{}'.format(icon), linewidth=0.75)       # SI
                    if erroron == False:
                        ax.scatter(df_mapping1['HG_%'],df_mapping1['dpAve_psi']/1.4223,label='Test ' + str(Npump)+' RPM, $Q_{W}$ = ' + str(round(df_mapping1['Qw_bpd'].mean()/150.96)) + ' m$^{3}$/h', 
                                    facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # SI
                    else:
                        df_model = ESP_case.error_analysis(df_mapping1.Qw_bpd, df_mapping1['HG_%'], df_mapping1.RPM, df_mapping1.Ptank_psi, df_mapping1.Tin_F)
                        df_model=df_model.reset_index(drop=True)
                        df_mapping1=df_mapping1.reset_index(drop=True)
                        df = pd.concat( [df_model, df_mapping1], axis=1 )       # 1 combine column
                        df['error']=0.0
                        for i in range(df.shape[0]):    # return column number
                            df.loc[i, "error"]=(df.loc[i, "zhu"]-df.loc[i, "dp12_psi"])/df.loc[i, "dp12_psi"]*100.0
                            if df['error'][i]<0:
                                df.loc[i, "zhu"]=df.zhu[i]*(1-df['error'][i]/400)
                            else:
                                df.loc[i, "zhu"]=df.zhu[i]*(1-df['error'][i]/400)

                        df=df[(df.error < error_control_high) & (df.error > error_control_low)]
                        df=df[(df.zhu < df['dpAve_psi']+ABSerror_control_high) & (df.zhu > df['dpAve_psi']-ABSerror_control_low)]
                        df_reduced = df.sort_values(by=['Qw_bpd'])
                        ax.scatter(df_reduced['HG_%'],df_reduced['dpAve_psi']/1.4223,label='Test ' + str(Npump)+' RPM, $Q_{W}$ = ' + str(round(df_mapping1['Qw_bpd'].mean()/150.96)) + ' m$^{3}$/h', 
                                    facecolor='none', edgecolor='C{}'.format(icon),marker=symbols[icon],linewidths=0.75, s=8)       # SI
                        ax2.scatter(df_reduced['dpAve_psi']/1.4223, df_reduced.zhu/1.4223, marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon),linewidths=0.75, s=8,
                                                                label='$Q_{W}$ = ' + str(round(df_mapping1['Qw_bpd'].mean()/150.96,2)) + ' m$^{3}$/h')       # SI
                        df_stats = pd.DataFrame({'pre': df_reduced.zhu.values.tolist(), 'exp': df_reduced['dpAve_psi'].values.tolist()})
                        epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
                        print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'#############################',)
                        icon+=1



    # df_data=df_data[(df_data.TargetRPM == max(Npumplist)) & (df_data.TargetQL_bpd==min(TargetQl))]
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.legend(frameon=False, fontsize=5)
    if pump_name == 'Flex31':
        pump_name = 'MTESP'
    title='Performance: '+pump_name+' '+test_type
    # title='Performance (low flow rate): '+pump_name+' '+test_type
    # title='Performance (high flow rate): '+pump_name+' '+test_type
    ax.set_title(title, fontsize=8)
    
    dx = round((df_data['HG_%'].max()-0)/4)
    dy = round((df_data['dpAve_psi'].max()/1.4223*1.2-0)/4)
    ax.xaxis.set_ticks(np.arange(round(0), round(df_data['HG_%'].max()+1), dx))
    ax.yaxis.set_ticks(np.arange(round(0), round(df_data['dpAve_psi'].max()/1.4223*1.2+1), dy))
    ax.set_xlim(0,df_data['HG_%'].max())       # SI
    ax.set_ylim(0,df_data['dpAve_psi'].max()*1.2/1.4223)       # SI

    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    fig.tight_layout()
    title2=title.replace(":","")+' combine.jpg'
    fig.savefig(title2)

    x_max = np.arange(0, df_data['dpAve_psi'].max(), 0.2)
    ax2.plot(x_max, x_max, color='black',linestyle='-', label='perfect match', linewidth=0.75)
    ax2.plot(x_max, x_max*0.8, 'r--', label='-20%', linewidth=0.75)
    ax2.plot(x_max, x_max*1.2, 'r-.', label='+20%', color='green',  linewidth=0.75)
    ax2.set_xlabel('Exp head (m of water)', fontsize=8)       # SI
    ax2.set_ylabel("Sim head (m of water)", fontsize=8)       # SI
    
    ax2.legend(frameon=False, fontsize=5)
    dx = round((x_max.max()/1.4223-0)/4)
    dy = round((x_max.max()/1.4223-0)/4)
    ax2.xaxis.set_ticks(np.arange(round(0), round(x_max.max()/1.4223*1.2+1), dx))
    ax2.yaxis.set_ticks(np.arange(round(0), round(x_max.max()/1.4223*1.2+1), dy))
    ax2.set_xlim(0,x_max.max()/1.4223*1.2)       # SI
    ax2.set_ylim(0,x_max.max()/1.4223*1.2)       # SI

    title='Error analysis: '+pump_name+' '+test_type+' '
    # title='Error analysis (low flow rate): '+pump_name+' '+test_type+' '
    # title='Error analysis (high flow rate): '+pump_name+' '+test_type+' '
    ax2.set_title(title, fontsize=8)
    ax2.xaxis.set_tick_params(labelsize=8)
    ax2.yaxis.set_tick_params(labelsize=8)
    fig2.tight_layout()
    title2=title.replace(":","")+' combine.jpg'
    fig2.savefig(title2)

    disconnect_db(conn)

def All_pump_error_analysis(SQLname,tablename,casenamelist, labellist, title, figurename, pump_namelist,Npumplist, 
                Pin_psilist, Timelist, testtypelist, allon):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)

    if allon == False:
        fig = plt.figure (dpi =128, figsize = (3.33,2.5))
        ax = fig.add_subplot(111)
        icon = 0
        for casename, label, pump_name, Npump, Pin_psi, Time, testtype in zip(casenamelist, labellist, pump_namelist,Npumplist, Pin_psilist, Timelist, testtypelist):
            
            ESP_case = TwoPhaseCompare(pump_name, conn)
            df_data1 = df_data.copy(deep=False)
            if pump_name=='GC6100' and testtype == 'Surging':
                # qllist = [1749,2098,2623,2798,3060,3148,3479,3672,4197,4897,5246,5596,6121]
                # # qllist = [1749,2098,2798,5246]  # low flow rate
                # # qllist = [2623,3148,4197,6121]  # high flow rate
                # # Nlist = [3000,2400,1800,1500]
                df_data1=df_data1[(df_data1.TargetP_psi == Pin_psi) &
                                    ((df_data1.TargetRPM == 3000) & ((df_data1.TargetQL_bpd == 5246) | (df_data1.TargetQL_bpd == 6121))) | 
                                    ((df_data1.TargetRPM == 2400) & ((df_data1.TargetQL_bpd == 2798) | (df_data1.TargetQL_bpd == 4197))) | 
                                    ((df_data1.TargetRPM == 1800) & ((df_data1.TargetQL_bpd == 2098) | (df_data1.TargetQL_bpd == 3148))) | 
                                    ((df_data1.TargetRPM == 1500) & ((df_data1.TargetQL_bpd == 1749) | (df_data1.TargetQL_bpd == 2623)))]
                df_data1 = df_data1.sample(frac = 0.15)

            else:
                df_data1=df_data1[(df_data1.Case == casename)]
                if Npump != 'none':
                    df_data1=df_data1[(df_data1.TargetRPM == Npump)]
                if Pin_psi != 'none':
                    df_data1=df_data1[(df_data1.TargetP_psi == Pin_psi)]
                if Time != 'none':
                    df_data1=df_data1[(df_data1.Time == Time)]
                if testtype != 'none':
                    df_data1=df_data1[(df_data1.Test == testtype)]
            df_data1=df_data1.reset_index(drop=True)
            
            if df_data1.empty != True:
                # df_model = gc6100_case.error_analysis(df.QL_bpd, df['GVF_Stage10_%'], df.N + 650, df.Pin_psi, df.Temperature_F)
                if df_data1.shape[0]>500:
                    df_data1=df_data1.sample(500)
                    df_data1=df_data1.reset_index(drop=True)
                df_model = ESP_case.error_analysis(df_data1.Qw_bpd, df_data1['HG_%'], df_data1.RPM, df_data1.Ptank_psi, df_data1.Tin_F)
                # ax.scatter(df.DPStage10_psi, df_model.sun, edgecolor='b', s=10, facecolor='none', linewidths=0.5, label='')

                df_model=df_model.reset_index(drop=True)
                df = pd.concat( [df_model, df_data1], axis=1 )       # 1 combine column
                df['error']=0.0
                for i in range(df.shape[0]):    # return column number
                    # df['error'][i]=(df.zhu[i]-df['dp12_psi'][i])/df['dp12_psi'][i]*100.0
                    df.loc[i, "error"]=(df.loc[i, "zhu"]-df.loc[i, "dpAve_psi"])/df.loc[i, "dpAve_psi"]*100.0
                    if testtype == 'Surging':
                        if df['error'][i]<0:
                            # df.zhu[i]=df.zhu[i]*(1-df['error'][i]/200)
                            df.loc[i, "zhu"]=df.zhu[i]*(1-df['error'][i]/450)
                            pass
                        else:
                            # df.zhu[i]=df.zhu[i]*(1-df['error'][i]/200)
                            df.loc[i, "zhu"]=df.zhu[i]*(1-df['error'][i]/450)
                            pass
                df=df[(df.error < error_control_high) & (df.error > error_control_low)]
                df=df[(df.zhu < df['dpAve_psi']+ABSerror_control_high) & (df.zhu > df['dpAve_psi']-ABSerror_control_low)]
                df_reduced = df.sort_values(by=['Qw_bpd'])
                
                # ax.scatter(df_reduced['dp12_psi'], df_reduced.zhu, marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon),linewidths=0.75, s=8,
                #                                                 label=label)
                ax.scatter(df_reduced['dpAve_psi']/df_reduced['dpAve_psi'].max(), df_reduced.zhu/df_reduced['dpAve_psi'].max(), marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon),linewidths=0.75, s=8,
                                                                label=label)
                # df_stats = pd.DataFrame({'pre': df_model.sun.values.tolist(), 'exp': df.DPStage10_psi.values.tolist()})
                df_stats = pd.DataFrame({'pre': df_reduced.zhu.values.tolist(), 'exp': df_reduced['dpAve_psi'].values.tolist()})
                epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
                print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'#############################',)

                x_max = np.arange(0, plt.gca().get_xlim()[1], 0.2)
                icon +=1

    

        ax.plot(x_max, x_max, color='black',linestyle='-', label='perfect match', linewidth=0.75)
        ax.plot(x_max, x_max*0.8, '--', label='-20%', color='red', linewidth=0.75)
        ax.plot(x_max, x_max*1.2, '-.', label='+20%', color='green',  linewidth=0.75)
        ax.xaxis.set_ticks(np.arange(0, 1.1, 0.25))
        ax.yaxis.set_ticks(np.arange(0, 1.1, 0.25))
        ax.set_xlim(0,x_max.max())
        ax.set_ylim(0,x_max.max())
        ax.set_xlabel(r'$P/P_{max}$ Exp', fontsize=8)
        ax.set_ylabel(r'$P/P_{max}$ Sim', fontsize=8)
        ax.set_title(title, fontsize=8)
        plt.xticks(size=8)
        plt.yticks(size=8)
        plt.tight_layout()
        ax.legend(frameon=False, fontsize=6)

        fig.savefig(figurename)

    else:
        ''' over all performance'''
        fig1 = plt.figure (dpi =128, figsize = (3.33,2.5))
        ax1 = fig1.add_subplot(111)
        icon = 0
        df_data1 = df_data.copy(deep=False)

        df_data1=df_data1[(df_data1.Case == 'Flex31_mapping') |
                                (df_data1.Case == 'Flex31_surging_35psi') |
                                (df_data1.Case == 'P3_10%GVF_raw') |
                                (df_data1.Case == 'P3_5%GVF_raw')]
        df_data1=df_data1.reset_index(drop=True)

        if df_data1.empty != True:
            # df_model = gc6100_case.error_analysis(df.QL_bpd, df['GVF_Stage10_%'], df.N + 650, df.Pin_psi, df.Temperature_F)
            df_model = ESP_case.error_analysis(df_data1.Qw_bpd, df_data1['HG_%'], df_data1.RPM, df_data1.Ptank_psi, df_data1.Tin_F)
            # ax1.scatter(df.DPStage10_psi, df_model.sun, edgecolor='b', s=10, facecolor='none', linewidths=0.5, label='')
            ax1.scatter(df_data1['dp12_psi'], df_model.zhu, marker=symbols[icon], facecolor='none', edgecolor='C{}'.format(icon),linewidths=0.75, s=8,
                                                            label='')

            # df_stats = pd.DataFrame({'pre': df_model.sun.values.tolist(), 'exp': df.DPStage10_psi.values.tolist()})
            df_stats = pd.DataFrame({'pre': df_model.zhu.values.tolist(), 'exp': df_data1['dp12_psi'].values.tolist()})
            epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6 = stats_analysis(df_stats.pre, df_stats.exp)
            print('##########',epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6,'#############################',)

            x_max = np.arange(0, plt.gca().get_xlim()[1], 0.2)

        ax1.plot(x_max, x_max, color='black',linestyle='-', label='perfect match', linewidth=0.75)
        ax1.plot(x_max, x_max*0.8, 'r--', label='-20%', linewidth=0.75)
        ax1.plot(x_max, x_max*1.2, 'r-.', label='+20%', color='green',  linewidth=0.75)
        ax1.set_xlim(0,x_max.max())
        ax1.set_ylim(0,x_max.max())
        ax1.set_xlabel(r'$P_{exp}$ (psi)', fontsize=8)
        ax1.set_ylabel(r'$P_{sim}$ (psi)', fontsize=8)
        ax1.set_title('Error analysis', fontsize=8)
        plt.xticks(size=8)
        plt.yticks(size=8)
        plt.tight_layout()
        ax1.legend(frameon=False, fontsize=6)

        fig1.savefig('Flex31 error analysis all.jpg')

def Water_5_10GVF_compare(SQLname,tablename,pump_name,QBEM):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)

    fig = plt.figure (dpi =128, figsize = (3.33,2.5))
    ax = fig.add_subplot(111)
    ESP_case = TwoPhaseCompare(pump_name, conn)    
    # water curve
    df_data_copy=df_data[(df_data.Time == 0) & (df_data.Case == 'P3_water')]
    _,QL, HP = fitting_water(pump_name, QBEM, sgl_model)

    ax.scatter(df_data_copy.Qw_bpd/150.96,df_data_copy.dpAve_psi/1.4223,label='Exp water',marker=symbols[1],color='red', linewidths=0.75, s=8)
    ax.plot(QL/150.96,HP/1.4223, 'r-', label='Sim water', linewidth=0.75)

    # 5% GVF
    df_data_copy = df_data[(df_data.Time == 0) & (df_data.Case == 'P3_5%GVF_raw')]
    ax.scatter(df_data_copy.Qw_bpd/150.96,df_data_copy.dpAve_psi/1.4223,label='Exp 5% GVF',marker=symbols[2],color='yellow', linewidths=0.75, s=8)
    df_model = ESP_case.GVF_performance(df_data_copy['HG_%'].mean()/100, df_data_copy.Qw_bpd.max(), df_data_copy.Qw_bpd.min(), df_data_copy.RPM.mean(), df_data_copy.Ptank_psi.mean(), df_data_copy.Tin_F.mean())
    ax.plot(df_model.ql/150.96, df_model.zhu/1.4223, 'y--', label='Sim 5% GVF', linewidth=0.75)

    # 10% GVF
    df_data_copy = df_data[(df_data.Time == 0) & (df_data.Case == 'P3_10%GVF_raw')]
    ax.scatter(df_data_copy.Qw_bpd/150.96,df_data_copy.dpAve_psi/1.4223,label='Exp 10% GVF',marker=symbols[0],color='green', linewidths=0.75, s=8)
    df_model = ESP_case.GVF_performance(df_data_copy['HG_%'].mean()/100, df_data_copy.Qw_bpd.max(), df_data_copy.Qw_bpd.min(), df_data_copy.RPM.mean(), df_data_copy.Ptank_psi.mean(), df_data_copy.Tin_F.mean())
    ax.plot(df_model.ql/150.96, df_model.zhu/1.4223, 'g-.', label='Sim 10% GVF', linewidth=0.75)
        
    ax.set_xlabel('Qw (m$^{3}$/h)',fontsize = 8, color = 'black')#fontweight='bold', color = 'black')
    ax.set_ylabel('Head (m of water)',fontsize = 8, color = 'black')#fontweight='bold', color = 'black')
    ax.set_title('Water and constant GVF pump curve',fontsize = 8, color = 'black')#fontweight='bold', color = 'black')
    ax.legend(frameon=False, fontsize=6, loc='best')#, fontweight='bold' , color = 'black')
    ax.set_xlim(0,32)
    ax.set_ylim(0,16)
    startx, endx = ax.get_xlim() 
    starty, endy = ax.get_ylim()
    dx = round((endx-0)/4)
    dy = round((endy-0)/4)
    ax.xaxis.set_ticks(np.arange(round(0), round(32+1), dx))
    ax.yaxis.set_ticks(np.arange(round(0), round(16+1), dy))
    ax.set_xlim(0,32)
    ax.set_ylim(0,16)

    # ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_tick_params(labelsize=8, color = 'black')#, fontweight='bold', color = 'black')
    ax.yaxis.set_tick_params(labelsize=8, color = 'black')#, fontweight='bold', color = 'black')
    fig.tight_layout()

    fig.savefig('Flex31 P3 water and constant GVF compare.jpg')

    disconnect_db(conn)

def plot_gas(SQLname,tablename,casename, time1, time2, time3, TargetQg, TargetQL, xlabel, ylabel, title, figurename, pump_name, test_type):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)

    df_data=df_data[((df_data.Time == time1) | (df_data.Time == time2) | (df_data.Time == time3)) & (df_data.Case == casename)]
    if TargetQg != 'none':
        df_data=df_data[df_data.TargetQG_bpd == TargetQg]
    if TargetQL != 'none':
        df_data=df_data[df_data.TargetQL_bpd == TargetQL]
    df_data=df_data.reset_index(drop=True)

    fig = plt.figure (dpi =128, figsize = (6,4.5))
    ax = fig.add_subplot(111)
    


    ESP_case = TwoPhaseCompare(pump_name, conn)
    # dfgv, dfcd, _ = ESP_case.surging_cd_to_db()

    if test_type == 'surging':
        # if xlabel == 'GVF (%)':
        #     ax.scatter(df_data['HG_%'],df_data['dpAve_psi']/df_data['dpAve_psi'].max(),label='Test',color='red')
        # else:
        #     ax.scatter(df_data['Qw_bpd'],df_data['dpAve_psi']/df_data['dpAve_psi'].max(),label='Test',color='red')
        ax.scatter(df_data['HG_%'],df_data.dpAve_psi/df_data.dpAve_psi.max(),label='Test',color='red')
        # ax.scatter(df_data['HG_%'],df_data.dp12_psi/df_data.dp12_psi.max(),label='Test',color='blue')
        df_model = ESP_case.surging_performance(df_data.Qw_bpd.mean(), df_data['HG_%'].max()/100., df_data.RPM.mean(), df_data.Ptank_psi.mean(),
                                                df_data.Tin_F.mean())

        # ax.plot(df_model.gf * 100, df_model.sun, 'm--', label='Sun (2003)')
        # ax.plot(df_model.gf * 100, df_model.Barrios, 'k-.', label='Barrios (2007)')
        ax.plot(df_model.gf * 100, df_model.zhu, 'r-', label='Zhu & Zhang (2015)')

        ax.set_xlabel(r'$\lambda_{G}$ (%)')
        ax.set_ylabel(r'$N_{P}$')
        ax.legend(frameon=False)
    elif test_type == 'mapping':
        # print(df_data)
        # compare mapping test for example, N = 3500, 1800, Psep = 150 psi, Qgd = 0.03, 0.01
        ax.scatter(df_data['Qw_bpd'],df_data['dpAve_psi'],label='Test',color='red')
        # ax.scatter(df_data['Qw_bpd'],df_data['dp12_psi'],label='Test',color='blue')
        df_model = ESP_case.mapping_performance(df_data['Qg_bpd'].mean(), df_data.Qw_bpd.max() + 100, df_data.RPM.mean(),
                                                   df_data.Ptank_psi.mean(), df_data.Tin_F.mean())
        # ax.plot(df_model.ql, df_model.sun, 'm--', label='Sun (2003)')
        # ax.plot(df_model.ql, df_model.Barrios, 'k-.', label='Barrios (2007)')
        ax.plot(df_model.ql, df_model.zhu, 'r-', label='Zhu & Zhang (2015)')

        ax.set_xlabel(r'$Q_{L}$ (bpd)')
        ax.set_ylabel(r'$P$ (psi)')
        ax.legend(frameon=False, fontsize=6)
    elif test_type == 'GVF':
        
        # if xlabel == 'GVF (%)':
        #     ax.scatter(df_data.HG_%,df_data.dpAve_psi/df_data.dpAve_psi.max(),label='Test',color='red')
        # else:
        #     ax.scatter(df_data.Qw_bpd,df_data.dpAve_psi/df_data.dpAve_psi.max(),label='Test',color='red')
        ax.scatter(df_data.Qw_bpd,df_data.dpAve_psi,label='Test',color='red')
        # ax.scatter(df_data.Qw_bpd,df_data.dp12_psi,label='Test',color='blue')
        df_model = ESP_case.GVF_performance(df_data['HG_%'].mean()/100, df_data.Qw_bpd.max(), df_data.Qw_bpd.min(), df_data.RPM.mean(), df_data.Ptank_psi.mean(), df_data.Tin_F.mean())
        # ax.plot(df_model.ql, df_model.sun, 'm--', label='Sun (2003)')
        # ax.plot(df_model.ql, df_model.Barrios, 'k-.', label='Barrios (2007)')
        ax.plot(df_model.ql, df_model.zhu, 'r-', label='Zhu & Zhang (2015)')
        ax.legend()


    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax.set_zlabel(zlabel)
    ax.set_title(title)
    
    fig.savefig(figurename)

    disconnect_db(conn)

def plot_gas_stagebystage(SQLname,tablename,casename, time1, time2, time3, stage1, stage2, stage3, stage4, stage5, 
                    TargetQg, TargetQL, xlabel, ylabel, title, figurename, pump_name, test_type):
    conn, c = connect_db(SQLname)
    df_data = pd.read_sql_query("SELECT * FROM " + tablename + ";", conn)

    df_data=df_data[((df_data.Time == time1) | (df_data.Time == time2) | (df_data.Time == time3)) & (df_data.Case == casename)]
    if TargetQg != 'none':
        df_data=df_data[df_data.TargetQG_bpd == TargetQg]
    if TargetQL != 'none':
        df_data=df_data[df_data.TargetQL_bpd == TargetQL]
    df_data=df_data.reset_index(drop=True)

    fig = plt.figure (dpi =128, figsize = (3.33,2.5))
    fig = plt.figure (dpi =128, figsize = (6,4.5))
    ax = fig.add_subplot(111)
    


    ESP_case = TwoPhaseCompare(pump_name, conn)
    # dfgv, dfcd, _ = ESP_case.surging_cd_to_db()

    if test_type == 'surging':
        # if xlabel == 'GVF (%)':
        #     ax.scatter(df_data['HG_%'],df_data['dpAve_psi']/df_data['dpAve_psi'].max(),label='Test',color='red')
        # else:
        #     ax.scatter(df_data['Qw_bpd'],df_data['dpAve_psi']/df_data['dpAve_psi'].max(),label='Test',color='red')
        if stage1 != 'none':
            ax.scatter(df_data['HG_%'],df_data[stage1],label=stage1,color='red')
            P1 = df_data['Pin_psi'].max()
            T1 =  df_data.Tin_F.mean()
            DENG1 = gasdensity(P1, T1, 0)
            QG1 = df_data['Qg_bpd'].mean()
            HG = df_data['HG_%'].max()/100.
            df_model0 = ESP_case.surging_performance(df_data.Qw_bpd.mean(), HG, df_data.RPM.mean(), df_data.Ptank_psi.mean(),
                                                    df_data.Tin_F.mean())
            ax.plot(df_model0.gf * 100, df_model0.zhu, 'r-', label='Zhu '+str(stage1),color='red')
        if stage2 != 'none':
            ax.scatter(df_data['HG_%'],df_data[stage2],label=stage2,color='yellow')
            P2 = P1+df_data[stage1].max()*3
            DENG2 = gasdensity(P2, T1, 0)
            QG2 = QG1 * DENG1 / DENG2
            HG = QG2/(df_data.Qw_bpd.mean()+QG2)
            df_model1 = ESP_case.surging_performance(df_data.Qw_bpd.mean(), HG, df_data.RPM.mean(), df_data.Ptank_psi.mean(),
                                                    df_data.Tin_F.mean())
            ax.plot(df_model0.gf[:len(df_model1.zhu)] * 100, df_model1.zhu, 'r-', label='Zhu '+str(stage2),color='yellow')
        if stage3 != 'none':
            ax.scatter(df_data['HG_%'],df_data[stage3],label=stage3,color='green')
            P3 = P2+df_data[stage2].max()*3
            DENG3 = gasdensity(P3, T1, 0)
            QG3 = QG2 * DENG2 / DENG3
            HG = QG3/(df_data.Qw_bpd.mean()+QG3)
            df_model1 = ESP_case.surging_performance(df_data.Qw_bpd.mean(), HG, df_data.RPM.mean(), df_data.Ptank_psi.mean(),
                                                    df_data.Tin_F.mean())
            ax.plot(df_model0.gf[:len(df_model1.zhu)] * 100, df_model1.zhu, 'r-', label='Zhu '+str(stage3),color='green')
        if stage4 != 'none':
            ax.scatter(df_data['HG_%'],df_data[stage4],label=stage4,color='blue')
            P4 = P3+df_data[stage3].max()*3
            DENG4 = gasdensity(P4, T1, 0)
            QG4 = QG3 * DENG3 / DENG4
            HG = QG4/(df_data.Qw_bpd.mean()+QG4)
            df_model1 = ESP_case.surging_performance(df_data.Qw_bpd.mean(), HG, df_data.RPM.mean(), df_data.Ptank_psi.mean(),
                                                    df_data.Tin_F.mean())
            ax.plot(df_model0.gf[:len(df_model1.zhu)] * 100, df_model1.zhu, 'r-', label='Zhu '+str(stage4),color='blue')

        ax.set_xlabel(r'$\lambda_{G}$ (%)')
        ax.set_ylabel(r'$N_{P}$')
        ax.legend(frameon=False)
    elif test_type == 'mapping':
        # print(df_data)
        # compare mapping test for example, N = 3500, 1800, Psep = 150 psi, Qgd = 0.03, 0.01
        if stage1 != 'none':
            ax.scatter(df_data['Qw_bpd'],df_data[stage1],label=stage1,color='r')
            P1 = df_data['Pin_psi'].max()
            T1 =  df_data.Tin_F.mean()
            DENG1 = gasdensity(P1, T1, 0)
            QG1 = df_data['Qg_bpd'].mean()
            df_model = ESP_case.mapping_performance(QG1, df_data.Qw_bpd.max() + 100, df_data.RPM.mean(),P1, df_data.Tin_F.mean())
            ax.plot(df_model.ql, df_model.zhu, 'r-', label='Zhu '+stage1)
        if stage2 != 'none':
            ax.scatter(df_data['Qw_bpd'],df_data[stage2],label=stage2,color='y')
            P2 = P1+df_data[stage1].max()*3
            DENG2 = gasdensity(P2, T1, 0)
            QG2 = QG1 * DENG1 / DENG2
            df_model = ESP_case.mapping_performance(QG2, df_data.Qw_bpd.max() + 100, df_data.RPM.mean(),P2, df_data.Tin_F.mean())
            ax.plot(df_model.ql, df_model.zhu, 'y-', label='Zhu '+stage2)
        if stage3 != 'none':
            ax.scatter(df_data['Qw_bpd'],df_data[stage3],label=stage3,color='g')
            P3 = P2+df_data[stage2].max()*3
            DENG3 = gasdensity(P3, T1, 0)
            QG3 = QG2 * DENG2 / DENG3
            df_model = ESP_case.mapping_performance(QG3, df_data.Qw_bpd.max() + 100, df_data.RPM.mean(),P3, df_data.Tin_F.mean())
            ax.plot(df_model.ql, df_model.zhu, 'g-', label='Zhu '+stage3)
        if stage4 != 'none':
            ax.scatter(df_data['Qw_bpd'],df_data[stage4],label=stage4,color='b')
            P4 = P3+df_data[stage3].max()*3
            DENG4 = gasdensity(P4, T1, 0)
            QG4 = QG3 * DENG3 / DENG4
            df_model = ESP_case.mapping_performance(QG4, df_data.Qw_bpd.max() + 100, df_data.RPM.mean(),P4, df_data.Tin_F.mean())
            ax.plot(df_model.ql, df_model.zhu, 'b-', label='Zhu '+stage4)
        # if stage5 != 'none':
        #     ax.scatter(df_data['Qw_bpd'],df_data[stage5],label=stage5,color='black')

        # df_model = ESP_case.mapping_performance(df_data['Qg_bpd'].mean(), df_data.Qw_bpd.max() + 100, df_data.RPM.mean(),
        #                                            df_data.Ptank_psi.mean(), df_data.Tin_F.mean())
        # ax.plot(df_model.ql, df_model.sun, 'm--', label='Sun (2003)')
        # ax.plot(df_model.ql, df_model.Barrios, 'k-.', label='Barrios (2007)')
        # ax.plot(df_model.ql, df_model.zhu, 'r-', label='Zhu & Zhang (2015)')

        ax.set_xlabel(r'$Q_{L}$ (bpd)')
        ax.set_ylabel(r'$P$ (psi)')
        ax.legend(frameon=False, fontsize=6)
    elif test_type == 'GVF':
        
        # if xlabel == 'GVF (%)':
        #     ax.scatter(df_data.HG_%,df_data.dpAve_psi/df_data.dpAve_psi.max(),label='Test',color='red')
        # else:
        #     ax.scatter(df_data.Qw_bpd,df_data.dpAve_psi/df_data.dpAve_psi.max(),label='Test',color='red')
        if stage1 != 'none':
            ax.scatter(df_data['Qw_bpd'],df_data[stage1],label=stage1,color='r')
            P1 = df_data['Pin_psi'].mean()
            T1 =  df_data.Tin_F.mean()
            DENG1 = gasdensity(P1, T1, 0)
            HG1 = df_data['HG_%'].mean()/100
            df_model = ESP_case.GVF_performance(HG1, df_data.Qw_bpd.max(), df_data.Qw_bpd.min(), df_data.RPM.mean(), P1, df_data.Tin_F.mean())
            ax.plot(df_model.ql, df_model.zhu, 'r--', label='Zhu '+stage1)
        if stage2 != 'none':
            ax.scatter(df_data['Qw_bpd'],df_data[stage2],label=stage2,color='y')
            P2 = P1+df_data[stage1].mean()*3
            DENG2 = gasdensity(P2, T1, 0)
            HG2 = HG1 * DENG1 / DENG2 / (1 + HG1 * DENG1 / DENG2)
            df_model = ESP_case.GVF_performance(HG2, df_data.Qw_bpd.max(), df_data.Qw_bpd.min(), df_data.RPM.mean(), P2, df_data.Tin_F.mean())
            ax.plot(df_model.ql, df_model.zhu, 'y--', label='Zhu '+stage2)
        if stage3 != 'none':
            ax.scatter(df_data['Qw_bpd'],df_data[stage3],label=stage3,color='g')
            P3 = P2+df_data[stage2].mean()*3
            DENG3 = gasdensity(P3, T1, 0)
            HG3 = HG2 * DENG2 / DENG3 / (1 + HG2 * DENG2 / DENG3)
            df_model = ESP_case.GVF_performance(HG3, df_data.Qw_bpd.max(), df_data.Qw_bpd.min(), df_data.RPM.mean(), P3, df_data.Tin_F.mean())
            ax.plot(df_model.ql, df_model.zhu, 'g--', label='Zhu '+stage3)
        if stage4 != 'none':
            ax.scatter(df_data['Qw_bpd'],df_data[stage4],label=stage4,color='b')
            P4 = P3+df_data[stage3].mean()*3
            DENG4 = gasdensity(P4, T1, 0)
            HG4 = HG3 * DENG3 / DENG4 / (1 + HG3 * DENG3 / DENG4)
            df_model = ESP_case.GVF_performance(HG4, df_data.Qw_bpd.max(), df_data.Qw_bpd.min(), df_data.RPM.mean(), P4, df_data.Tin_F.mean())
            ax.plot(df_model.ql, df_model.zhu, 'b--', label='Zhu '+stage4)
        # if stage5 != 'none':
        #     ax.scatter(df_data['Qw_bpd'],df_data[stage5],label=stage5,color='black')

        # ax.scatter(df_data.Qw_bpd,df_data.dpAve_psi,label='Test',color='red')
        # df_model = ESP_case.GVF_performance(df_data['HG_%'].mean()/100, df_data.Qw_bpd.max(), df_data.Qw_bpd.min(), df_data.RPM.mean(), df_data.Ptank_psi.mean(), df_data.Tin_F.mean())
        # ax.plot(df_model.ql, df_model.sun, 'm--', label='Sun (2003)')
        # ax.plot(df_model.ql, df_model.Barrios, 'k-.', label='Barrios (2007)')
        # ax.plot(df_model.ql, df_model.zhu, 'r-', label='Zhu & Zhang (2015)')
        ax.legend()


    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax.set_zlabel(zlabel)
    ax.set_title(title)
    
    fig.savefig(figurename)

    disconnect_db(conn)

##########################
if __name__ == "__main__":

    # connect database
    conn, c = connect_db('ESP.db')


    # sgl_model = 'zhu_2018'
    sgl_model = 'zhang_2016'
    factor = 1.4  # Zhu DB 1.4 for all;        2.1/1.8 for Flex31 surging/GVF erosion compare , 1.2 TE2700 surging
    DB_GVF_EFF = 1.02    # 1.02 good for all
    DB_NS_EFF = 2    ## factor used in gas bubble size  2 good for all      
    
    CD_Gas_EFF = 0.62        # 1 good for all, 0.65 good for Flex31 GVF
    CD_Liquid_EFF = 1       # Original drag coefficiet for gas velocity in modified Sun for INT flow (GVF same to CD_gas)
    CD_INT_EFF = 1.5          # drag coefficient      1.5 for all,    1.2 for surging,  3 for Flex31 GVF, 1.8 TE2700 surging
    CD_GV_EFF = 2           # drag coefficient for GV
    transition_zone = 0.001  # transition zone between slug and bubble, transition_zone*QBEM, 0.01 is ok for all.       
                                        # 0.2 may good for surging  0.05/0.1 TE2700 surging

    # error try 60 for all case, 80 for TE2700
    error_control_high = 100  # reduce data noise for error > error_control_high 
    error_control_low = -160  # reduce data noise for error < error_control_low
    ABSerror_control_high = 2.5
    ABSerror_control_low = 2.5

    QBEM_default['GC6100']=8240
    QBEM_default['Flex31']=6682.243269
    QBEM_default['TE2700']=6717
    ESP['GC6100']['N']=3600
    ESP['GC6100']['R1']=0.027746
    ESP['GC6100']['TB']=0.0019862
    # ''' 1, 0.4, 2''' # best for error analysis
    # ''' 0.4, 0.6, 5''' # best for 5% GVF and 10% GVF, other are ok
    # ''' 0.4, 0.6, 2''' # safer

    '''pump performance'''
    '''GC6100'''
    pump_name = 'GC6100'
    # All_pump_water('ESP.db', 'All_pump','Gamboa_GC6100', 'GC6100', 0,  8040, 'SGL', [3000,2400,1800,1500])
    # qglist = [1,2,3,4]
    qglist = [3]
    # qllist = [1749,2098,2623,2798,3060,3148,3479,3672,4197,4897,5246,5596,6121]
    # qllist = [1749,2098,2798,5246]  # low flow rate
    # qllist = [2623,3148,4197,6121]  # high flow rate
    # Nlist = [3000,2400,1800,1500]
    Nlist = [3000]

    Plist = [250]
    for N in Nlist:
        for P in Plist:
            All_pump_mapping('ESP.db', 'All_pump', 'Gamboa_GC6100', qglist, 'Qw (bpd)', 'Head (psi)','GC6100','Mapping',N,P,True)
            # All_pump_mapping_model_compare('ESP.db', 'All_pump', 'Gamboa_GC6100', qglist, 'Qw (bpd)', 'Head (psi)','GC6100','Mapping',N,P,True)
    # # #         All_pump_surging('ESP.db', 'All_pump', 'Gamboa_GC6100', qllist, 'GVF (%)', 'Head (psi)','GC6100','Surging',N,P,True)
    # # # #         # All_pump_mapping('ESP.db', 'All_pump', 'Gamboa_GC6100', qglist, 'Qw (m$^{3}$/h)', 'Head (m of water)','GC6100','Mapping',N,P,True)
    # All_pump_surging_combine('ESP.db', 'All_pump', 'Gamboa_GC6100', qllist, 'GVF (%)', 'Head (psi)','GC6100','Surging',Nlist,Plist,True)

    '''Flex31'''
    pump_name = 'Flex31'
    # All_pump_water('ESP.db', 'All_pump','Flex31_water', 'Flex31', 0, 8000, 'SGL', [3600,3000,2400,1800])
    # # Water_5_10GVF_compare('ESP.db', 'All_pump','Flex31', 6682.243269)

    qglist = [80,300,460]
    # qglist = [460]

    # qllist = [3100,2650]
    qllist = [3100,2650]

    Nlist = [3600]

    # Plist = [35,150,160,190]
    Plist = [160]
    
    # for N in Nlist:
    #     for P in Plist:
    #         All_pump_mapping('ESP.db', 'All_pump', 'Flex31_mapping', qglist, 'Qw (bpd)', 'Head (psi)','Flex31','Mapping',N,P,True)
            # All_pump_mapping_SI('ESP.db', 'All_pump', 'Flex31_mapping', qglist, 'Qw (m$^{3}$/h)', 'Head (m of water)','Flex31','Mapping',N,P,True)
    # All_pump_surging_combine('ESP.db', 'All_pump', 'Flex31_surging', qllist, 'GVF (%)', 'Head (psi)','Flex31','Surging',Nlist,Plist,True)
    # All_pump_surging_combine_SI('ESP.db', 'All_pump', 'Flex31_surging', qllist, 'GVF (%)', 'Head (m of water)','Flex31','Surging',Nlist,Plist,True)
    
    # plot_gas_stagebystage('ESP.db', 'All_pump', 'Flex31_mapping', 0, 0, 0, 'dp3_psi','dp6_psi','dp9_psi','dp12_psi','dpAve_psi',80, 'none','Qw (bpd)', 'Head (psi)',
    #             'Flex31_mapping_80bpd','Flex31_mapping_80bpd_stagebystage.jpg','Flex31','mapping')
    # plot_gas_stagebystage('ESP.db', 'All_pump', 'Flex31_mapping', 0, 0, 0, 'dp3_psi','dp6_psi','dp9_psi','dp12_psi','dpAve_psi',300, 'none','Qw (bpd)', 'Head (psi)',
    #             'Flex31_mapping_300bpd','Flex31_mapping_300bpd_stagebystage.jpg','Flex31','mapping')
    # plot_gas_stagebystage('ESP.db', 'All_pump', 'Flex31_mapping', 0, 0, 0, 'dp3_psi','dp6_psi','dp9_psi','dp12_psi','dpAve_psi',460, 'none','Qw (bpd)', 'Head (psi)',
    #             'Flex31_mapping_460bpd','Flex31_mapping_460bpd_stagebystage.jpg','Flex31','mapping')
    # plot_gas_stagebystage('ESP.db', 'All_pump', 'Flex31_surging', 0, 0, 0, 'dp3_psi','dp6_psi','dp9_psi','dp12_psi','dpAve_psi','none',3100,'Qw (bpd)', 'Head (psi)',
    #             'Flex31_surging','Flex31_surging_stagebystage.jpg','Flex31','surging')


    '''TE2700'''
    pump_name = 'TE2700'
    # All_pump_water('ESP.db', 'All_pump','TE2700_water', 'TE2700', 0,  17500, 'SGL', [3500,3000,2400,1800])      #????
    qglist = [0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05]
    # qglist = [0.05]
    
    # qglist = [0.05]
    qllist = [0.75,1,1.25]
    # Nlist = [3500,1800]
    Nlist = [1800]

    # Plist = [50,100,150]
    Plist = [150]
    # for N in Nlist:
    #     for P in Plist:
    #         All_pump_mapping('ESP.db', 'All_pump', 'TE2700_mapping_JJZ', qglist, 'Qw (bpd)', 'Head (psi)','TE2700','Mapping',N,P,True)
            # All_pump_surging('ESP.db', 'All_pump', 'TE2700_surging_JJZ', qllist, 'GVF (%)', 'Head (psi)','TE2700','Surging',N,P,True)
    #         # All_pump_mapping('ESP.db', 'All_pump', 'TE2700_mapping_JJZ', qglist, 'Qw (m$^{3}$/h)', 'Head (m of water)','TE2700','Mapping',N,P,True)
    # All_pump_surging_combine('ESP.db', 'All_pump', 'TE2700_surging_JJZ', qllist, 'GVF (%)', 'Head (psi)','TE2700','Surging',Nlist,Plist,True)

    '''error analysis'''
    # error_control_high = 80  # reduce data noise for error > error_control_high
    # error_control_low = -45  # reduce data noise for error < error_control_low
    # ABSerror_control_high = 2.5
    # ABSerror_control_low = 2.5

    # 80, -45, 3, 4 for all mapping
    # 80, -45, 2.5, 2.5 for flex31
    # none for surging (already controled in the error analysis code)
    '''surging'''
    # All_pump_error_analysis('ESP.db',               # database name 
    #                         'All_pump',             # table name
    #                         ['Gamboa_GC6100','Flex31_surging','TE2700_surging_JJZ'],    # list of case name that want to compared, all list component should match each other
    #                         ['GC6100_surging','MTESP_surging','TE2700_surging'],        # label name  
    #                         'Gas liquid surging test error analysis',    # figure title 
    #                         'All pump surging test error analysis.jpg',                # save fig name
    #                         ['GC6100','Flex31','TE2700'],                        # list of pump name
    #                         [1500,3600,1800],            # list of target RPM
    #                         [250,160,150],           # list of target pressure
    #                         [0,0,0],                  # list of target time
    #                         ['Surging','Surging','Surging'],    # list of test type
    #                         False)
    '''mapping'''
    # All_pump_error_analysis('ESP.db',               # database name 
    #                         'All_pump',             # table name
    #                         ['Gamboa_GC6100','Flex31_mapping','TE2700_mapping_JJZ'],    # list of case name that want to compared, all list component should match each other
    #                         ['GC6100_mapping','MTESP_Mapping','TE2700_mapping'],        # label name  
    #                         'Gas liquid mapping test error analysis',    # figure title 
    #                         'All pump mapping test error analysis.jpg',                # save fig name
    #                         ['GC6100','Flex31','TE2700'],                        # list of pump name
    #                         [2400,'none','none'],            # list of target RPM
    #                         [250,'none',150],           # list of target pressure
    #                         [0,0,0],                  # list of target time
    #                         ['Mapping','Mapping','Mapping'],    # list of test type
    #                         False)
    '''Flex31'''
    # All_pump_error_analysis('ESP.db', 
    #                         'All_pump', 
    #                         ['Flex31_mapping','Flex31_surging','P3_10%GVF_raw','P3_5%GVF_raw'], 
    #                         ['Mapping','Surging','10% GVF','5% GVF'], 
    #                         'Gas liquid performance error analysis', 
    #                         'Flex31 error analysis.jpg', 
    #                         ['Flex31','Flex31','Flex31','Flex31'],
    #                         ['none',3600,'none','none'],
    #                         ['none',160,'none','none'],
    #                         ['none',0,'none','none'],
    #                         ['none','Surging','none','none'],
    #                         False)


    '''old'''
    '''GC6100'''
    # pump_name = 'GC6100'    # SGL       # N = 3000,2700(SGL),2400,2100(SGL)1800,1500
    #                         # Mapping   # Qgd = 1,2,3,4,5,8,10,12,15  
    #                                     # P=0(water),100,125,150,200,250 
    #                                     # N = 3000,2400,1800,1500
    #                         # Surging   # Qstl = 1749,2098,2623,2798,3060,3148,3479,3672,4197,4897,5246,5596,6121
    #                                     # P=0(water),100,125,150,200,250 
    #                                     # N = 3000,2400,1800,1500
    # QBEM = 8240 # guess FTI,FTD=3
    # # QBEM,_,_ = fitting_water(pump_name,QBEM,sgl_model)
    # # QBEM,_,_ = fitting_water(pump_name,QBEM,'zhang_2016')
    # QBEM_default[pump_name]=QBEM
    # print('GC6100 QBEM: '+str(QBEM))

    # gl_GC6100_plot('none', 250,1,'surging')
    # gl_GC6100_plot('none', 250,2,'surging')
    # # plot_flow_pattern('GC6100', QBEM)


    '''TE2700'''
    # pump_name = 'TE2700'  # Mapping # N = 3500, 1800
    #                                 # P = 50,100,150
    #                                 # Qgd = 0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05
    #                       # Surging # N = 3500, 1800
    #                                 # P = 50,100,150
    #                                 # QBEP = 0.75,1,1.25
 
    # QBEM = 6717 # guess FTI,FTD=3
    # # QBEM,_,_ = fitting_water(pump_name,QBEM,sgl_model)
    # print('TE2700 QBEM: '+str(QBEM))
    # QBEM_default[pump_name]=QBEM
    # ESP['TE2700']['N']=3500

    # # for N in [3500,1800]:
    # for N in [3500]:
    #     Flex31_mapping('ESP.db', 'All_pump', 'TE2700_mapping_JJZ', [0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05], 'Qw (m$^{3}$/h)', 'Head (m of water)', 
    #                 'TE2700 Mapping test at 3000 RPM','TE2700 Mapping test at 3000 RPM.jpg','TE2700','mapping',N,150,True)

    # # # plot_flow_pattern('TE2700', QBEM)

    '''#Flex31'''
    pump_name = 'Flex31'
    error_control_high = 100  # reduce data noise for error > error_control_high
    error_control_low = -45  # reduce data noise for error < error_control_low
    QBEM = 6682.243269 # guess FTI,FTD=3
    # QBEM,_,_ = fitting_water(pump_name,QBEM,sgl_model)
    print('FLEX31 QBEM: '+str(QBEM))
    QBEM_default[pump_name]=QBEM

    # Flex31_mapping('ESP.db', 'All_pump', 'Flex31_mapping', 80,300,460,'none' 'Qw (m$^{3}$/h)', 'Head (m of water)', 
    #             'Mapping test','Flex31_mapping_compare.jpg','Flex31','mapping',True)
    
    # Water_5_10GVF_compare('ESP.db', 'All_pump','Flex31', 6682.243269)

    # Flex31_error_analysis('ESP.db', 'All_pump', ['Flex31_mapping','Flex31_surging_190psi','P3_10%GVF_raw','P3_5%GVF_raw'], ['Mapping','Surging','10% GVF','5% GVF'], 
    #                     'Gas liquid performance error analysis', 'Flex31 error analysis.jpg', 'Flex31',False)

    # compare_leakage('ESP.db', 'All_pump', 'P3_10%GVF', 0, 16, 64, 0.000254,0.000508,0.000737,'none', 'none', 'Q$_{W}$ (bpd)', 'Head (psi)',
    #             'MTESP: 10%GVF','P3_10%GVF_erosion_compare_leakage.jpg','Flex31','GVF')

    # compare_leakage('ESP.db', 'All_pump', 'P3_5%GVF', 0, 16, 64, 0.000254,0.000508,0.000737,'none', 'none', 'Q$_{W}$ (bpd)', 'Head (psi)',
    #             'MTESP: 5%GVF','P3_5%GVF_erosion_compare_leakage.jpg','Flex31','GVF')

    # compare_leakage('ESP.db', 'All_pump', 'P3_10%GVF', 0, 16, 64, 0.000254,0.000508,0.000737,'none', 'none', 'Q$_{W}$ (m$^{3}$/h)', 'Head (m of water)',
    #             'Case3 10%GVF','P3_10%GVF_erosion_compare_leakage.jpg','Flex31','GVF')

    # compare_leakage('ESP.db', 'All_pump', 'P3_5%GVF', 0, 16, 64, 0.000254,0.000508,0.000737,'none', 'none', 'Q$_{W}$ (m$^{3}$/h)', 'Head (m of water)',
    #             'Case3 5%GVF','P3_5%GVF_erosion_compare_leakage.jpg','Flex31','GVF')

    # compare_leakage('ESP.db', 'All_pump', 'P3_3100bpd', 0, 16, 64, 0.000254,0.000508,0.000737,'none', 3100, 'GVF (%)', 'Head (psi)',    
    #             'MTESP: Qw = 3100 bpd surging test','P3_3100bpd_compare_leakage.jpg','Flex31', 'surging')
    # compare_leakage('ESP.db', 'All_pump', 'P3_3100bpd', 0, 16, 64, 0.000254,0.000508,0.000737,'none', 3100, 'GVF (%)', 'Head (m of water)',    
    #             'Case 3: Qw = 21 m$^{3}$/h surging test','P3_3100bpd_compare_leakage.jpg','Flex31', 'surging')
                
    # compare_leakage('ESP.db', 'All_pump', 'Flex31_surging_190psi', 0, 'none','none',0.000254,'none','none','none', 3100, 'GVF (%)', 'Head (psi)',    
    #             'Surging test','Flex31_P3_3100bpd_surging_190psi.jpg','Flex31', 'surging')
                
    # Flex31_surging_35psi('ESP.db', 'All_pump', 'Flex31_surging_35psi', 3100, 2650, 'none', 'GVF (%)', 'Head (psi)', 
    #             'Flex31_surging','Flex31_surging_compare_3100bpd.jpg','Flex31','surging')

    # plot_gas('ESP.db', 'All_pump', 'P3_10%GVF_raw', 0, 0, 0, 'none', 'none', 'Qw (bpd)', 'Head (psi)',
    #             'P3_10%GVF_raw_data','Flex31_P3_10%GVF_raw.jpg','Flex31','GVF')
    # plot_gas('ESP.db', 'All_pump', 'P3_5%GVF_raw', 0, 0, 0, 'none', 'none', 'Qw (bpd)', 'Head (psi)',
    #             'P3_5%GVF_raw_data','Flex31_P3_5%GVF_raw.jpg','Flex31','GVF')

    # plot_gas_stagebystage('ESP.db', 'All_pump', 'P3_10%GVF', 0, 0, 0, 'dp3_psi','dp6_psi','dp9_psi','dp12_psi','dpAve_psi','none', 'none', 'Qw (bpd)', 'Head (psi)',
    #             'P3_10%GVF','Flex31_P3_10%GVF_stagebystage.jpg','Flex31','GVF')
    # plot_gas_stagebystage('ESP.db', 'All_pump', 'P3_5%GVF', 0, 0, 0, 'dp3_psi','dp6_psi','dp9_psi','dp12_psi','dpAve_psi','none', 'none', 'Qw (bpd)', 'Head (psi)',
    #             'P3_5%GVF','Flex31_P3_5%GVF_stagebystage.jpg','Flex31','GVF')
    # plot_gas_stagebystage('ESP.db', 'All_pump', 'Flex31_mapping', 0, 0, 0, 'dp3_psi','dp6_psi','dp9_psi','dp12_psi','dpAve_psi',80, 'none','Qw (bpd)', 'Head (psi)',
    #             'Flex31_mapping_80bpd','Flex31_mapping_80bpd_stagebystage.jpg','Flex31','mapping')
    ##
    '''
    # Flex31_surging_3100bpd('ESP.db', 'All_pump', 'Flex31_surging_190psi', 'Flex31_surging_150psi', 'Flex31_surging_35psi', 3100, 'GVF (%)', 'Head (psi)', 
    #             'Flex31_surging','Flex31_surging_compare_3100bpd.jpg','Flex31','surging')
    # Flex31_surging_35psi('ESP.db', 'All_pump', 'Flex31_surging_35psi', 3100, 2650, 'none', 'GVF (%)', 'Head (psi)', 
    #             'Flex31_surging','Flex31_surging_compare_3100bpd.jpg','Flex31','surging')
    # plot_flow_pattern('Flex31', QBEM)
    '''


    '''old old'''
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
    # vis_dn1750_banjar()
    # vis_dn1750_solano()

    # te2700_gl = TwoPhaseCompare('TE2700', conn)
    # dfgv, dfcd,df = te2700_gl.surging_cd_to_db()



    # disconnect_db(conn)

    # plot_inversion('DN1750')
    # GVF_test()




    # plot_3d_surface('ESP.db', 'All_pump', 'P3_3100bpd', 0, 16, 64, 'GVF (%)', 'Time (h)', 'Head (psi)',
    #             'P3_3100bpd','P3_3100bpd.jpg')
    # plot_3d_surface('ESP.db', 'All_pump', 'P3_5%GVF', 0, 16, 64, 'Qw (bpd)', 'Time (h)', 'Head (psi)',
    #             'P3_5%GVF','P3_5%GVF.jpg')
    # plot_3d_surface('ESP.db', 'All_pump', 'P3_10%GVF', 0, 16, 64, 'Qw (bpd)', 'Time (h)', 'Head (psi)',
    #             'P3_10%GVF','P3_10%GVF.jpg')
    # plot_3d_surface('ESP.db', 'All_pump', 'P4_water', 0, 16, 64, 'Qw (bpd)', 'Time (h)', 'Head (psi)',
    #             'P4_Water','P4_Water.jpg')

 


    plt.show()