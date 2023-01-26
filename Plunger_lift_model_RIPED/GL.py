# -*- coding: utf-8 -*-
"""
Reference:  Unified Model for Gas-Liquid Pipe Flow Via Slug Dynamics
            #TUALP 3P 2012v1 * May 24, 2012
Developer: Jianjun Zhu
Date: Jan 12, 2019
"""

import numpy as np

# global variables
FEC = 0.90      # maximum entrainment fraction in gas core
HLSC = 0.24     # maximum liquid holdup in slug body
G = 9.81        # gravitational acceleration (m/s2)
PI = np.pi      # ratio of the circumference of a circle to its diameter
STAW = 0.0731   # water surface tension against air (N/m)
DENA = 1.2      # density of air at atmospheric pressure (kg/m3)
E1 = 0.0005     # tolerance for iterations
FIC = 0.0142    # critical interfacial friction factor
AXP = 0.        # cross sectional area of pipe (m2)

"""
This global is set anywhere when the input data has triggered a possible error. A non-zero value triggers a write of the
input data to the dump file for later analysis of the error
"""
Ierr = 0
IFFM = 0        # interfacial friction indicator

# get_fff - get fanning friction factor
def get_fff(re, ed):
    """
    :param re:  Reynolds number
    :param ed:  relative wall roughness
    :return:    Colebrook frictin factor correlation
    """
    lo = 800.
    hi = 1200.
    Fl = 0.
    Fh = 0.

    if re < hi: Fl = 16.0 / re
    if re > lo:
        AA = (2.457 * np.log(1.0 / ((7.0 / re) ** 0.9 + 0.27 * ed))) ** 16.0
        BB = (37530.0 / re) ** 16.0
        Fh = 2.0 * ((8.0 / re) ** 12.0 + 1.0 / (AA + BB) ** 1.5) ** (1.0 / 12.0)

    if re < lo:
        return Fl
    elif re > hi:
        return Fh
    else:
        return (Fh * (re - lo) + Fh * (hi - re)) / (hi - lo)


# wet_frc_biberg - Calculate wetted wall fraction assuming flat interface based on Biberg (1999) approximate method
def wet_fra_biberg(hx):
    third = 1.0 / 3.0
    eps = 1.E-06
    a = 0.
    h = 0.
    if (abs(hx) >= 0) and (abs(hx) <= eps):                 a = 0
    if (abs(hx) > eps) and (abs(hx) < (1 - eps)):           h = abs(hx)
    if (abs(hx) > (1. - eps)) and (abs(hx) <= (1. + eps)):  a = 1.
    if abs(hx) > (1. + eps):                                h = 1. / abs(hx)

    if (h > 0) and (h < 1):
        a = h + 0.533659 * (1.0 - 2.0 * h + h ** third - (1.0 - h) ** third) - h * (1.0 - h) * (1.0 - 2.0 * h) \
            * (1.0 + 4.0 * (h * h + (1.0 - h) ** 2.0)) / 200.0 / PI

    if a > 1.:                  a = 1.
    if (a < 0.01) and (a < hx): a = hx

    return a


# Calculate Wetted wall fraction according to Zhang and Sarica, SPEJ (2009)
def wet_fra_zhang_sarica(ang, deng, denl, vc, vf, d, hlf, th0):
    th = 1.0

    if abs(ang) >= 85. * PI / 180.:
        return th
    else:
        FRGF = deng * (vc - vf) ** 2 / (denl - deng) / G / d / np.cos(ang)
        FROT = 0.7 * np.sqrt((1.0 - hlf) / hlf)
        RFR = FRGF / FROT

        if RFR > 20: return th

        YO = -d * np.sin(PI * th0) ** 3 / (3.0 * PI * hlf)
        YI = d * 0.25 * (hlf ** 1.2 - 1.0)
        RAY = YO / YI

        if RAY < 0.000001: RAY = 0.000001

        COO = np.log(abs(RAY))
        RYD = abs(YI - YO) / d

        if RYD < 0.000001: return th

        xx = COO * RFR ** 1.4

        if abs(xx) > 600: return th

        YY = YO / np.exp(xx)

        if (RFR > 20) or (RYD < 0.000001): return th

        th = (1.0 + th0) / 2.0 + (1.0 - th0) * np.tan(PI * (2.0 * YY - YI - YO) / (3.0 * (YI - YO))) / 3.464

        if (th > 1.) or (YY > YI): th = 1.0
    return th


# Single Phase Flow Calculation
# SGL  =  calculates pressure gradient for single phase flow of liquid or gas
def SGL(D, ED, ANG, P, DEN, V, VIS):
    """
    :param D: pipe diameter (m)
    :param ED: relative pipe wall roughness
    :param ANG: angle of pipe from hor (rad)
    :param P: pressure (Pa)
    :param DEN: density (kg/m3)
    :param V: velocity (m/s)
    :param VIS: viscosity (Pa-s)
    :param PGT: total pressure gradient
    :param PGF: friction pressure gradient
    :param PGG: gravity pressure gradient
    :return: FF, PGT, PGF, PGG, PGA
    """
    # Calculate elevation pressure gradient
    PGG = -DEN * np.sin(ANG) * G

    # Calculate frictional pressure gradient
    RE = abs(D * DEN * V / VIS)
    FF = get_fff(RE, ED)
    PGF = -2.0 * FF * DEN * V * V / D

    # Calculate acceleration pressure gradient
    if DEN <= 400:
        EKK = DEN * V * V / P
        ICRIT = 0
        if EKK > 0.95:
            ICRIT = 1
        if ICRIT == 1:
            EKK = 0.95
        PGT = (PGG + PGF) / (1.0 - EKK)
        PGA = PGT * EKK
    else:
        PGA = 0.0

    PGT = (PGG + PGF + PGA)

    return FF, PGT, PGF, PGG, PGA


# DISLUG - calculates the superficial liquid velocity on the boundary between dispersed bubble flow and slug flow
#          with a given superficial gas velocity
#       -- slug/bubbly transition, constant vsg
def DISLUG(D, ED, ANG, VSG, DENL, DENG, VISL, STGL):
    """
    :param D: pipe diameter (m)
    :param ED: relative pipe wall roughness
    :param ANG: angle of pipe from hor (rad)
    :param VSG: gas superficial velocity (m/s)
    :param DENL: liquid density (kg/m3)
    :param DENG: gas density (kg/m3)
    :param VISL: liquid viscosity (Pa-s)
    :param STGL: gas viscosity (Pa-s)
    :return: VDB, Ierr
    """
    global Ierr, AXP, IFFM
    Ierr = 0
    CC = 1.25 - 0.5 * abs(np.sin(ANG))
    VMC = VSG / (1.0 - HLSC)

    # guess a VDB and VM
    VDB1 = 2.0
    VM = 2 * VSG

    for icon in np.arange(1, 501):
        VM = VDB1 + VSG
        HLB = VDB1 / VM
        DENM = (1.0 - HLB) * DENG + HLB * DENL
        REM = abs(DENM * D * VM / VISL)
        FM = get_fff(REM, ED)

        if REM > 5000.0:
            VMN = np.sqrt(abs(((1.0 / HLB - 1.0) * 6.32 * CC * np.sqrt(abs(DENL - DENG) * G * STGL)) / (FM * DENM)))
        else:
            VMN = np.sqrt(abs(((1.0 / HLB - 1.0) * 6.32 * CC * np.sqrt(abs(DENL - DENG) * G * STGL))
                              * np.sqrt(5000.0 / REM) / (FM * DENM)))

        if VMN < VMC: VMN = VMC

        ABM = abs((VMN - VM) / VM)

        if ABM < E1: break

        VM = (VMN + 4.0 * VM) / 5.0
        VDB1 = VM - VSG

        if icon > 500: Ierr = 4

    VDB = VM - VSG

    return VDB, Ierr


# BUSLUG - calculates the superficial gas velocity on the boundary between slug flow and bubbly flow
#          with a given superficial liquid velocity (for near vertical upward flow, >60 deg, and large D)
#       -- slug/bubbly transition, constant vsl (for near vertical upward flow, >60 deg, and large D)
def BUSLUG(D, ANG, VSL):
    """
    :param D: pipe diameter (m)
    :param ANG: angle of pipe from hor (rad)
    :param VSL: liquid superficial velocity (m/s)
    :return: VBU
    """
    VO = (0.54 * np.cos(ANG) + 0.35 * np.sin(ANG)) * np.sqrt(G * D)
    HGC = 0.25
    VBU = VSL * HGC / (1.0 - HGC) + VO * HGC
    return VBU


# STSLUG  =  calculates the superficial liquid velocity on the boundary between slug flow and stratified
#            (or annular)  flow with a given superficial gas velocity (for horizontal and downward flow)
#       -- slug/stratified-annular transition, constant vsg (for horizontal and downward flow)
def STSLUG(D, ED, ANG, VSG, DENL, DENG, VISL, VISG, STGL, IFFM):
    """
    :param D: pipe diameter (m)
    :param ED: relative pipe wall roughness
    :param ANG: angle of pipe from hor (rad)
    :param VSG: gas superficial velocity (m/s)
    :param DENL: liquid density (kg/m3)
    :param DENG: gas density (kg/m3)
    :param VISL: liquid viscosity (Pa-s)
    :param VISG: gas viscosity (Pa-s)
    :param IFFM: interfacial friction indicator
    :return: VST, Ierr
    """
    global Ierr, AXP
    CS = (32.0 * np.cos(ANG) ** 2 + 16.0 * np.sin(ANG) ** 2) * D
    CC = 1.25 - 0.5 * abs(np.sin(ANG))
    VDB, Ierr = DISLUG(D, ED, ANG, VSG, DENL, DENG, VISL, STGL)

    # Guess a VST
    VST = 0.5
    VM = VST + VSG
    HLS = 1.0 / (1.0 + abs(VM / 8.66) ** 1.39)

    if HLS < HLSC: HLS = HLSC

    FE = 0.0
    HLF = VST / VM
    VF = VM
    VC = VM
    REMX = 5000.0
    RESG = abs(DENG * VSG * D / VISG)
    WEB = abs(DENG * VSG * VSG * D / STGL)
    FRO = abs(np.sqrt(G * D) / VSG)
    VSGT = 5.0 * np.sqrt(DENA / DENG)

    # Initial interfacial absolute roughness
    EAI = D / 7.0
    FI = FIC
    FI1, FI2 = 0, 0

    for icon in np.arange(501):
        if VST > VDB:
            break

        # Entrainment fraction according to Oliemans et al.'s (1986) correlation
        RESL = abs(DENL * VST * D / VISL)
        CCC = 0.003 * WEB ** 1.8 * FRO ** 0.92 * RESL ** 0.7 * (DENL / DENG) ** 0.38 * \
              (VISL / VISG) ** 0.97 / RESG ** 1.24
        FEN = CCC / (1.0 + CCC)

        if FEN > FEC: FEN = FEC
        FE = (FEN + 9.0 * FE) / 10.0

        # Translational velocity according to Nicklin (1962), Bendiksen (1984) and Zhang et al. (2000)
        if REMX < 2000.0:
            VAV = 2.0 * VM
        elif REMX > 4000.0:
            VAV = 1.3 * VM
        else:
            VAV = (2.0 - 0.7 * (REMX - 2000.0) / 2000.0) * VM

        VT = VAV + (0.54 * np.cos(ANG) + 0.35 * np.sin(ANG)) * np.sqrt(G * D * abs(DENL - DENG) / DENL)
        HLFN = ((HLS * (VT - VM) + VST) * (VSG + VST * FE) - VT * VST * FE) / (VT * VSG)

        if HLFN <= 0.0: HLFN = abs(HLFN)
        if HLFN >= 1.0: HLFN = 1.0 / HLFN

        HLF = (HLFN + 9.0 * HLF) / 10.0
        HLC = (1.0 - HLF) * VST * FE / (VM - VST * (1.0 - FE))

        if HLC < 0.0: HLC = 0.0
        AF = HLF * AXP
        AC = (1.0 - HLF) * AXP

        # Calculate wetted wall fraction
        TH0 = wet_fra_biberg(HLF)

        # Wetted wall fraction according to Zhang and Sarica, SPEJ (2011)
        TH = wet_fra_zhang_sarica(ANG, DENG, DENL, VC, VF, D, HLF, TH0)

        # Wetted perimeters
        SF = PI * D * TH
        SC = PI * D * (1.0 - TH)
        AB = D * D * (PI * TH - np.sin(2.0 * TH * PI) / 2.0) / 4.0
        SI = (SF * (AB - AF) + D * np.sin(PI * TH) * AF) / AB

        # The hydraulic diameters
        DF = 4.0 * AF / (SF + SI)
        THF = 2.0 * AF / (SF + SI)
        DC = 4.0 * AC / (SC + SI)
        VC = (VM - VST * (1.0 - FE)) / (1.0 - HLF)

        # Reynolds numbers
        DENC = (DENL * HLC + DENG * (1.0 - HLF - HLC)) / (1.0 - HLF)
        REF = abs(DENL * VF * DF / VISL)
        REC = abs(DENG * VC * DC / VISG)

        # Friction factors
        FF = get_fff(REF, ED)
        FC = get_fff(REC, ED)

        # Interfacial friction factor:
        # Stratified flow interfacial friction factor
        if D <= 0.127:
            if IFFM == 1:
                FI1 = FC * (1.0 + 15.0 * abs(2.0 * THF / D) ** 0.5 * (VSG / VSGT - 1.0))
            else:
                FI1 = FC * (1.0 + 21.0 * abs(THF / D) ** 0.72 * abs(VSG / VSGT - 1.0) ** 0.8)
        else:
            # Interfacial friction factor according to Baker et al. (1988)
            WEE = DENG * VF * VF * EAI / STGL
            VIS = VISL * VISL / (DENL * STGL * EAI)
            WVM = WEE * VIS

            if WVM <= 0.005:
                EAI = 34.0 * STGL / (DENG * VF * VF)
            else:
                EAI = 170.0 * STGL * WVM ** 0.3 / (DENG * VF * VF)

            if EAI > THF / 4.0: EAI = THF / 4.0

            EDI = EAI / D
            FI2 = get_fff(REC, EDI)

        FIN = (0.127 * FI1 / D + FI2) / (1.0 + 0.127 / D)
        if FIN < FC: FIN=FC
        if FIN > 0.1: FIN = 0.1
        FI = (FIN + 9.0 * FI) / 10.0

        ABCD = (SC * FC * DENG * VC * abs(VC) / (2.0 * AC) + SI * FI * DENG * (VC - VF) * abs(VC - VF)
                * (1.0 / AF + 1.0 / AC) / 2.0 - (DENL - DENC) * G * np.sin(ANG)) * AF * 2.0 / (SF * FF * DENL)
        if ABCD < 0:
            VFN = VF * 0.9
        else:
            VFN = np.sqrt(ABCD)

        ABU = abs((VFN - VF) / VF)

        if ABU <= E1:
            VST = VFN * HLF / (1.0 - FE)
            break

        VF = (VFN + 9.0 * VF) / 10.0
        VST = VF * HLF / (1.0 - FE)
        VM = VST + VSG
        DPEX = (DENL * (VM - VF) * (VT - VF) * HLF + DENC * (VM - VC) * (VT - VC) * (1.0 - HLF)) * D / CS / 4.0
        REMX = abs(D * VM * DENL / VISL)
        FM = get_fff(REMX, ED)
        DPSL = FM * DENL * VM * VM / 2.0
        DPAL = DPSL + DPEX

        if REMX < 5000.: DPAL = DPAL * REMX / 5000.0

        AD = DPAL / (3.16 * CC * np.sqrt(STGL * abs(DENL - DENG) * G))
        HLSN = 1.0 / (1.0 + AD)

        if HLSN < HLSC: HLSN = HLSC

        HLS = (HLSN + 4.0 * HLS) / 5.0

        if icon >= 500: Ierr = 2

    if VST > VDB: VST = VDB

    return VST, Ierr


# ANSLUG - calculates the superficial gas velocity on the boundary between slug flow and annular (or stratified)
#          flow with a given superficial liquid velocity (for upward flow)
#       -- slug/stratified-annular transition, constant vsl (for upward flow)
def ANSLUG(D, ED, ANG, VSL, DENL, DENG, VISL, VISG, STGL, IFFM):
    """
    :param D: pipe diameter (m)
    :param ED: relative pipe wall roughness
    :param ANG: angle of pipe from hor (rad)
    :param VSL: liquid superficial velocity (m/s)
    :param DENL: liquid density (kg/m3)
    :param DENG: gas density (kg/m3)
    :param VISL: liquid viscosity (Pa-s)
    :param VISG: gas viscosity (Pa-s)
    :param IFFM: interfacial friction indicator
    :return: VAN, Ierr
    """
    global Ierr, AXP

    Ierr = 0
    CS = (32.0 * np.cos(ANG) ** 2 + 16.0 * np.sin(ANG) ** 2) * D
    CC = 1.25 - 0.5 * abs(np.sin(ANG))
    V24 = 19.0 * VSL / 6.0

    # guess a VAN
    VAN = 10.
    VM = VSL + VAN
    HLS = 1.0 / (1.0 + abs(VM / 8.66) ** 1.39)

    if HLS < HLSC: HLS = HLSC

    FE = 0.0
    HLF = VSL / VM
    VC = VAN / (1.0 - HLF)
    VF = VSL / HLF
    FI = FIC
    REMX = 5000.0
    RESL = DENL * VSL * D / VISL
    VSGT = 5.0 * np.sqrt(DENA / DENG)
    ABCD = V24 ** 2

    # Initial interfacial absolute roughness
    EAI = D / 7.0
    FI1, FI2 = 0, 0

    for icon in np.arange(1001):
        # Entrainment fraction according to Oliemans et al's (1986) correlation
        WEB = DENG * VAN * VAN * D / STGL
        FRO = np.sqrt(G * D) / VAN
        RESG = DENG * VAN * D / VISG
        CCC = 0.003 * WEB ** 1.8 * FRO ** 0.92 * RESL ** 0.7 * (DENL / DENG) ** 0.38 * \
              (VISL / VISG) ** 0.97 / RESG ** 1.24
        FEN = CCC / (1.0 + CCC)

        if FEN > 0.75: FEN = 0.75

        FE = FEN

        # Translational velocity according to Nicklin (1962), Bendiksen (1984) and Zhang et al. (2000)
        if REMX < 2000.0:
            VAV = 2.0 * VM
        elif REMX > 4000.0:
            VAV = 1.3 * VM
        else:
            VAV = (2.0 - 0.7 * (REMX - 2000.0) / 2000.0) * VM

        VT = VAV + (0.54 * np.cos(ANG) + 0.35 * np.sin(ANG)) * np.sqrt(G * D * abs(DENL - DENG) / DENL)
        HLFN = ((HLS * (VT - VM) + VSL) * (VAN + VSL * FE) - VT * VSL * FE) / (VT * VAN)

        if HLFN <= 0.0: HLFN = abs(HLFN)
        if HLFN >= 1.0: HLFN = 1.0 / HLFN

        HLF = HLFN
        HLC = (1.0 - HLF) * VSL * FE / (VM - VSL * (1.0 - FE))

        if HLC < 0.0: HLC = 0.0

        AF = HLF * AXP
        AC = (1.0 - HLF) * AXP

        # Calculate wet wall fraction
        TH0 = wet_fra_biberg(HLF)

        # Wetted wall fraction according to Zhang and Sarica, SPEJ (2011)
        TH = wet_fra_zhang_sarica(ANG, DENG, DENL, VC, VF, D, HLF, TH0)

        # Wet perimeters
        SF = PI * D * TH
        SC = PI * D * (1.0 - TH)
        AB = D * D * (PI * TH - np.sin(2.0 * TH * PI) / 2.0) / 4.0
        SI = (SF * (AB - AF) + D * np.sin(PI * TH) * AF) / AB

        # The hydraulic diameters
        DF = 4.0 * AF / (SF + SI)
        THF = 2.0 * AF / (SF + SI)
        DC = 4.0 * AC / (SC + SI)
        VFN = VSL * (1.0 - FE) / HLF
        VF = (VFN + 9.0 * VF) / 10.0

        # Reynolds numbers
        DENC = (DENL * HLC + DENG * (1.0 - HLF - HLC)) / (1.0 - HLF)
        REF = abs(DENL * VF * DF / VISL)
        REC = abs(DENG * VC * DC / VISG)

        # Frictional factors
        FF = get_fff(REF, ED)
        FC = get_fff(REC, ED)

        if D <= 0.127:
            # Interfacial friction factor (stratified) according to Andritsos et al. (1987)
            # Modified by Zhang (2001)
            if IFFM == 1:
                FI1 = FC * (1.0 + 15.0 * abs(2.0 * THF / D) ** 0.5 * (VAN / VSGT - 1.0))
            else:
                # Use Fan's correlation (2005)
                FI1 = FC * (1.0 + 21.0 * abs(THF / D) ** 0.72 * abs(VAN / VSGT - 1.0) ** 0.8)

        else:
            # Interfacial friction factor according to Baker et al. (1988)
            WEE = DENG * VF * VF * EAI / STGL
            VIS = VISL * VISL / (DENL * STGL * EAI)
            WVM = WEE * VIS

            if WVM <= 0.005:
                EAI = 34.0 * STGL / (DENG * VF * VF)
            else:
                EAI = 170.0 * STGL * WVM ** 0.3 / (DENG * VF * VF)

            if EAI > THF / 4.0: EAI = THF / 4.0

            EDI = EAI / D
            FI2 = get_fff(REC, EDI)

        FIN = (0.127 * FI1 / D + FI2) / (1.0 + 0.127 / D)

        if FIN < FC: FIN = FC
        if FIN > 0.1: FIN = 0.1
        FI = (FIN + 9.0 * FI) / 10.0

        ABCDN = (SF * FF * DENL * VF * VF / (2.0 * AF) - SC * FC * DENG * VC * VC / (2.0 * AC)
                 + (DENL - DENC) * G * np.sin(ANG)) * 2.0 / (SI * FI * DENG * (1.0 / AF + 1.0 / AC))

        ABCD = (ABCDN + 9.0 * ABCD) / 10.0

        if ABCD < 0.0:
            VCN = VC * 0.9
        else:
            VCN = np.sqrt(ABCD) + VF

        if VCN < V24: VCN = V24

        ABU = abs((VCN - VC) / VC)
        VC = VCN

        if ABU < E1: break

        VAN = VC * (1.0 - HLF) - VSL * FE

        if VAN < 0: VAN = -VAN

        VM = VSL + VAN
        DPEX = (DENL * (VM - VF) * (VT - VF) * HLF + DENC * (VM - VC) * (VT - VC) * (1.0 - HLF)) * D / CS / 4.0
        REMX = abs(D * VM * DENL / VISL)
        FM = get_fff(REMX, ED)
        DPSL = FM * DENL * VM * VM / 2.0
        DPAL = DPSL + DPEX

        if REMX < 5000: DPAL = DPAL * np.sqrt(REMX / 5000.0)

        AD = DPAL / (3.16 * CC * np.sqrt(STGL * abs(DENL - DENG) * G))
        HLSN = 1.0 / (1.0 + AD)

        if HLSN < HLSC: HLSN = HLSC

        HLS = HLSN

        if icon >= 1000: Ierr = 3

    VAN = VC * (1.0 - HLF) - VSL * FE

    return VAN, Ierr


# DBFLOW - calculates pressure gradient and liquid holdup for dispersed bubble flow (without bubble rise velocity)
def DBFLOW(D, ED, ANG, VSG, VSL, DENL, DENG, VISL):
    """
    :param D: pipe diameter (m)
    :param ED: relative pipe wall roughness
    :param ANG: angle of pipe from hor (rad)
    :param VSG: gas superficial velocity (m/s)
    :param DENL: liquid density (kg/m3)
    :param DENG: gas density (kg/m3)
    :param VISL: liquid viscosity (Pa-s)
    :param VISG: gas viscosity (Pa-s)
    :return: FM, HL, PGT, PGA, PGF, PGG
    """
    VM = VSG + VSL

    # Calculate liquid holdup
    HL = VSL / (VSG + VSL)
    DENM = (1.0 - HL) * DENG + HL * DENL
    DENS = DENM + HL * (DENL - DENM) / 3.0
    REM = abs(DENS * D * VM / VISL)
    FM = get_fff(REM, ED)
    PGF = -2.0 * FM * DENS * VM ** 2 / D
    PGG = -G * DENM * np.sin(ANG)
    PGA = 0.0
    PGT = (PGF + PGG + PGA)
    return FM, HL, PGT, PGA, PGF, PGG


# BUFLOW - calculates pressure gradient and liquid holdup for bubbly flow (with bubble rise velocity Vo)
def BUFLOW(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, STGL):
    """
    :param D: pipe diameter (m)
    :param ED: relative pipe wall roughness
    :param ANG: angle of pipe from hor (rad)
    :param VSG: gas superficial velocity (m/s)
    :param DENL: liquid density (kg/m3)
    :param DENG: gas density (kg/m3)
    :param VISL: liquid viscosity (Pa-s)
    :param STGL: gas-liquid surface tension (N/m)
    :return: FM, HL, PGT, PGA, PGF, PGG
    """
    VM = VSG + VSL
    VO = 1.53 * abs(G * (DENL - DENG) * STGL / DENL / DENL) ** 0.25 * np.sin(ANG)

    # Calculate liquid holdup
    if abs(ANG) < 10. * PI / 180.:
        HL = VSL / (VSG + VSL)
    else:
        HL = (np.sqrt(abs((VM - VO) ** 2 + 4.0 * VSL * VO)) - VSG - VSL + VO) / (2.0 * VO)

    # CALCULATE PRESSURE GRADIENTS
    DENM = (1.0 - HL) * DENG + HL * DENL
    DENS = DENM + HL * (DENL - DENM) / 3.0
    REM = abs(DENS * D * VM / VISL)

    FM = get_fff(REM, ED)
    PGF = -2.0 * FM * DENS * VM ** 2 / D
    PGG = -G * DENM * np.sin(ANG)
    PGA = 0.
    PGT = (PGF + PGG + PGA)
    return FM, HL, PGT, PGA, PGF, PGG


# ITFLOW - calculates pressure gradient, liquid holdup and slug characteristics for intermittent flow
def ITFLOW(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL):
    """
    :param D: pipe diameter (m)
    :param ED: relative pipe wall roughness
    :param ANG: angle of pipe from hor (rad)
    :param VSG: gas superficial velocity (m/s)
    :param DENL: liquid density (kg/m3)
    :param DENG: gas density (kg/m3)
    :param VISL: liquid viscosity (Pa-s)
    :param VISG: gas viscosity (Pa-s)
    :param STGL: gas-liquid surface tension (N/m)
    :return: HL, FF, PGT, PGA, PGF, PGG, FGL, HLF, CU, CS, CF, VF, VC, FQN, RSU, HLS, HLS, ICON, IFGL
    """

    global Ierr
    Ierr = 0
    VM = VSL + VSG
    HLS = 1.0 / (1.0 + abs(VM / 8.66) ** 1.39)

    if HLS < HLSC: HLS = HLSC

    # Translational velocity according to Nicklin (1962), Bendiksen (1984) and Zhang et al. (2000)
    REMM = abs(DENL * D * VM / VISL)

    if REMM < 2000:
        VAV = 2.0 * VM
    elif REMM > 4000:
        VAV = 1.3 * VM
    else:
        VAV = (2.0 - 0.7 * (REMM - 2000.0) / 2000.0) * VM

    VT = VAV + (0.54 * np.cos(ANG) + 0.35 * np.sin(ANG)) * np.sqrt(G * D * abs(DENL - DENG) / DENL)

    # slug length
    CS = (32.0 * np.cos(ANG) ** 2 + 16.0 * np.sin(ANG) ** 2) * D
    CE = 1.25 - 0.5 * abs(np.sin(ANG))
    VSGT = 5.0 * np.sqrt(DENA / DENG)
    VMT = (VM + VT) / 2.0

    # Guess CU and CF
    CU = CS * VM / VSL
    CF = CU - CS
    HLF = (VT - VM) / VT
    VF = VSL / 2.0
    VFN1 = VSL / 2.0
    VC = VM
    FF = 0.1
    FI = FIC

    # Initial interfacial absolute roughness
    EAI = D / 7.0
    FGL = 'INT'
    IFGL = 4
    FI1, FI2 = 0, 0
    icon = 1

    for icon in np.arange(1, 2001):
        # Overall liquid holdup
        HL = (HLS * (VT - VM) + VSL) / VT
        HLF = HLS * (VT - VM) / (VT - VF)

        if HLF > HL:  HLF = 0.99 * HL

        VCN = (VM - HLF * VF) / (1.0 - HLF)

        if VCN < 0.0:   VCN = -VCN
        if VCN > VT:    VCN = VT

        VC = VCN
        CFN = (HLS - HL) * CS / (HL - HLF)
        ABU = abs((CFN - CF) / CF)
        CF = (CFN + 19.0 * CF) / 20.0
        AF = HLF * AXP
        AC = (1.0 - HLF) * AXP

        # Slug liquid holdup
        DPEX = (DENL * (VM - VF) * (VT - VF) * HLF + DENG * (VM - VC) * (VT - VC) * (1.0 - HLF)) * D / CS / 4.0
        REM = abs(DENL * VM * D / VISL)
        FM = get_fff(REM, ED)
        DPSL = FM * DENL * VM * VM / 2.0
        DPAL = DPSL + DPEX

        if REM < 5000:  DPAL = DPAL * (REM / 5000.0)

        AD = DPAL / (3.16 * CE * np.sqrt(STGL * abs(DENL - DENG) * G))
        HLSN = 1.0 / (1.0 + AD)

        if HLSN < HLSC: HLSN = HLSC

        HLS = HLSN

        if (REM < 1500) and (HLS < 0.6):   HLS = 0.6

        if (VSL / VM > HLS) or (abs(CF) < D):
            FGL = "D-B"
            IFGL = 2
            break

        # Wetted wall fraction assuming flat film surface
        TH0 = wet_fra_biberg(HLF)

        # Wetted wall fraction according to Zhang and Sarica, SPEJ (2011)
        TH = wet_fra_zhang_sarica(ANG, DENG, DENL, VC, VF, D, HLF, TH0)

        # Wetted perimeters
        SF = PI * D * TH
        SC = PI * D * (1.0 - TH)
        AB = D * D * (PI * TH - np.sin(2.0 * TH * PI) / 2.0) / 4.0
        SI = (SF * (AB - AF) + D * np.sin(PI * TH) * AF) / AB

        # The hydraulic diameters
        DF = 4.0 * AF / (SF + SI)
        THF = 2.0 * AF / (SF + SI)
        DC = 4.0 * AC / (SC + SI)

        # Frictional factors
        REF = abs(DENL * VF * DF / VISL)
        REC = abs(DENG * VC * DC / VISG)
        FFN = get_fff(REF, ED)
        FC = get_fff(REC, ED)
        FF = (FFN + 9.0 * FF) / 10.0
        VSGF = VC * (1.0 - HLF)

        if D <= 0.127:
            if IFFM == 1:
                # Interfacial friction factor (stratified) according to Andritsos et al. (1987) Modified by Zhang (2001)
                FI1 = FC * (1.0 + 15.0 * abs(2.0 * THF / D) ** 0.5 * (VSGF / VSGT - 1.0))
            else:
                # Use Fan's correlation (2005)
                FI1 = FC * (1.0 + 21.0 * abs(THF / D) ** 0.72 * abs(VSGF / VSGT - 1.0) ** 0.8)
        else:

            # Interfacial friction factor according to Baker et al. (1988)
            WEE = DENG * VF * VF * EAI / STGL
            VIS = VISL * VISL / (DENL * STGL * EAI)
            WVM = WEE * VIS

            if WVM <= 0.005:
                EAI = 34.0 * STGL / (DENG * VF * VF)
            else:
                EAI = 170.0 * STGL * WVM ** 0.3 / (DENG * VF * VF)

            if EAI > THF / 4.0: EAI = THF / 4.0

            EDI = EAI / D
            FI2 = get_fff(REC, EDI)

        FRS = (0.127 * FI1 / D + FI2) / (1.0 + 0.127 / D)

        # Interfacial friction factor (annular) according to Ambrosini et al. (1991)
        REG = abs(VC * DENG * D / VISG)
        WED = DENG * VC * VC * D / STGL
        FIF = 0.046 / REG ** 0.2
        SHI = abs(FI * DENG * (VC - VF) ** 2 / 2.0)
        THFO = THF * np.sqrt(abs(SHI * DENG)) / VISG
        FRA = FIF * (1.0 + 13.8 * (THFO - 200.0 * np.sqrt(DENG / DENL)) * WED ** 0.2 / REG ** 0.6)
        FRA1 = get_fff(REC, THF / D)

        if FRA > FRA1:  FRA = FRA1

        RTH = SF / THF
        FIN = (50.0 * FRS / RTH + FRA) / (1.0 + 50.0 / RTH)

        if FIN < FC:    FIN = FC
        if FIN > 0.1:   FIN = 0.1

        FI = (FIN + 9.0 * FI) / 10.0

        # Calculate film length CF using the combined momentum equation
        FSL = (DENL * (VM - VF) * (VT - VF) - DENG * (VM - VC) * (VT - VC)) / CF
        ABCD = (FSL + SC * FC * DENG * VC * abs(VC) / 2.0 / AC +
                SI * FI * DENG * (VC - VF) * abs(VC - VF) / 2.0 * (1.0 / AF + 1.0 / AC) -
                (DENL - DENG) * G * np.sin(ANG)) * 2.0 * AF / (SF * FF * DENL)

        if ABCD > 0:
            VFN = np.sqrt(ABCD)
            if VFN > VM:    VFN = VM
        else:
            VFN = -np.sqrt(-ABCD)
            if VFN < -VM:   VFN = -VM

        ABV = abs((VFN - VF) / VF)
        VF = (VFN + 19.0 * VF) / 20.0

        if (ABU < E1) or (ABV < E1): break
        if icon >= 2000:
            Ierr = 1
            break

    # Slug unit length
    CU = CF + CS
    DENM = DENL * HLS + DENG * (1.0 - HLS)
    DENS = DENM + HLS * (DENL - DENM) / 3.0
    RES = abs(DENS * VM * D / VISL)
    FS = get_fff(RES, ED)
    FQN = VT / CU  # slug frequency
    RSU = CS / CU  # slug to slug unit length ratio

    # Pressure gradient in slug
    FOS = RSU * (DENL * (VM - VF) * (VT - VF) * HLF + DENG * (VM - VC) * (VT - VC) * (1.0 - HLF)) / CS
    DPS = -FS * DENS * VM * VM * 2.0 / D - DENM * G * np.sin(ANG) - FOS

    # Pressure gradient in film
    FOF = FOS * CS / CF
    DPF = FOF - SF * FF * DENL * VF * abs(VF) / (2.0 * AXP) - DENG * FC * SC * \
          VC * abs(VC) / (2.0 * AXP) - (DENL * HLF + DENG * (1.0 - HLF)) * G * np.sin(ANG)

    # Gravitational pressure gradient
    PGG = -(DENL * HL + DENG * (1.0 - HL)) * G * np.sin(ANG)

    # Frictional pressure gradient
    PGF = -((FS * DENS * (VF + HLS * (VM - VF)) ** 2.0 * 2.0 / D) * CS +
                  (SF * FF * DENL * VF * abs(VF) / (2.0 * AXP) +
                  SC * FC * DENG * VC * abs(VC) / (2.0 * AXP)) * CF) / CU

    # Acceleration pressure gradient
    PGA = 0.0

    # Total pressure gradient
    PGT = PGG + PGF + PGA

    return HL, FF, PGT, PGA, PGF, PGG, FGL, HLF, CU, CS, CF, VF, VC, FQN, RSU, HLS, icon, IFGL


# SAFLOW  =  calculates pressure gradient and liquid holdup for stratified or annular flow
def SAFLOW(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL, P):
    """
    :param D: pipe diameter (m)
    :param ED: relative pipe wall roughness
    :param ANG: angle of pipe from hor (rad)
    :param VSG: gas superficial velocity (m/s)
    :param DENL: liquid density (kg/m3)
    :param DENG: gas density (kg/m3)
    :param VISL: liquid viscosity (Pa-s)
    :param VISG: gas viscosity (Pa-s)
    :param STGL: gas-liquid surface tension (N/m)
    :param P: pressure (P)
    :return: HL, FE, FF, PGT, PGA, PGF, PGG, FGL, HLF, VF, SF, THF, ICON, IFGL
    """
    global Ierr
    Ierr = 0
    FI = FIC
    # Initial interfacial absolute roughness
    EAI = D / 7.0
    FGL = 'STR'
    IFGL = 5
    FI1, FI2 = 0, 0
    ICON = 1
    FF = 0

    # Entrainment fraction according to Oliemans et al's (1986) correlation
    RESG = abs(DENG * VSG * D / VISG)
    WEB = abs(DENG * VSG * VSG * D / STGL)
    FRO = abs(VSG / np.sqrt(G * D))
    RESL = abs(DENL * VSL * D / VISL)
    CCC = 0.003 * WEB ** 1.8 * FRO ** (-0.92) * RESL ** 0.7 * (DENL / DENG) ** 0.38 \
          * (VISL / VISG) ** 0.97 / RESG ** 1.24
    FE = CCC / (1.0 + CCC)

    if FE > FEC: FE = FEC

    VSGT = 5.0 * np.sqrt(DENA / DENG)

    # Guess a film velocity
    VF = VSL
    VC = VSL + VSG
    HLF = VSL / (VSL + VSG)
    ABCD = 5.0
    TH = 0.5

    for ICON in np.arange(1, 2001):
        if VF > 0:
            HLFN = VSL * (1.0 - FE) / VF
        else:
            HLFN = 1.2 * HLF

        if HLFN >= 1.:
            HLFN = 1.0 - 1.0 / HLFN

        HLF = (HLFN + 19.0 * HLF) / 20.0

        if VF > 0:
            VCN = (VSG + FE * VSL) / (1.0 - HLF)
        else:
            VCN = (VSG + VSL - VF * HLF) / (1.0 - HLF)

        VC = (VCN + 9.0 * VC) / 10.0
        AF = HLF * AXP
        AC = (1.0 - HLF) * AXP

        if VF > 0:
            HLC = VSL * FE / VC
        else:
            HLC = (VSL - VF * HLF) / VC

        if HLC < 0: HLC = 0

        DENC = (DENL * HLC + DENG * (1.0 - HLF - HLC)) / (1.0 - HLF)

        # Wetted wall fraction assuming flat film surface
        TH0 = wet_fra_biberg(HLF)

        # Wetted wall fraction according to Zhang and Sarica, SPEJ (2011)
        TH = wet_fra_zhang_sarica(ANG, DENG, DENL, VC, VF, D, HLF, TH0)

        # Wetted perimeters
        SF = PI * D * TH
        SC = PI * D * (1.0 - TH)
        AB = D * D * (PI * TH - np.sin(2.0 * TH * PI) / 2.0) / 4.0
        SI = (SF * (AB - AF) + D * np.sin(PI * TH) * AF) / AB

        # The hydraulic diameters
        DF = 4.0 * AF / (SF + SI)
        THF = 2.0 * AF / (SF + SI)
        DC = 4.0 * AC / (SC + SI)

        # Frictional factors
        REF = abs(DENL * VF * DF / VISL)
        REC = abs(DENG * VC * DC / VISG)
        FFN = get_fff(REF, ED)
        FC = get_fff(REC, ED)
        FF = (FFN + 9.0 * FF) / 10.0
        VSGF = VC * (1.0 - HLF)

        if D <= 0.127:
            if IFFM == 1:
                # Interfacial friction factor (stratified) according to Andritsos et al. (1987) Modified by Zhang (2001)
                FI1 = FC * (1.0 + 15.0 * abs(2.0 * THF / D) ** 0.5 * (VSGF / VSGT - 1.0))
            else:
                # Use Fan's correlation (2005)
                FI1 = FC * (1.0 + 21.0 * abs(THF / D) ** 0.72 * abs(VSGF / VSGT - 1.0) ** 0.8)
        else:

            # Interfacial friction factor according to Baker et al. (1988)
            WEE = DENG * VF * VF * EAI / STGL
            VIS = VISL * VISL / (DENL * STGL * EAI)
            WVM = WEE * VIS

            if WVM <= 0.005:
                EAI = 34.0 * STGL / (DENG * VF * VF)
            else:
                EAI = 170.0 * STGL * WVM ** 0.3 / (DENG * VF * VF)

            if EAI > THF / 4.0: EAI = THF / 4.0

            EDI = EAI / D
            FI2 = get_fff(REC, EDI)

        FRS = (0.127 * FI1 / D + FI2) / (1.0 + 0.127 / D)

        # Interfacial friction factor (annular) according to Ambrosini et al. (1991)
        REG = abs(VC * DENG * D / VISG)
        WED = DENG * VC * VC * D / STGL
        FIF = 0.046 / REG ** 0.2
        SHI = abs(FI * DENG * (VC - VF) ** 2 / 2.0)
        THFO = THF * np.sqrt(abs(SHI * DENG)) / VISG
        FRA = FIF * (1.0 + 13.8 * (THFO - 200.0 * np.sqrt(DENG / DENL)) * WED ** 0.2 / REG ** 0.6)
        FRA1 = get_fff(REC, THF / D)

        if FRA > FRA1:  FRA = FRA1

        RTH = SF / THF
        FIN = (50.0 * FRS / RTH + FRA) / (1.0 + 50.0 / RTH)

        if FIN < FC:    FIN = FC
        if FIN > 0.1:   FIN = 0.1

        FI = (FIN + 9.0 * FI) / 10.0

        ABCD = (SC * FC * DENG * VC * abs(VC) / 2.0 / AC +
                SI * FI * DENG * (VC - VF) * abs(VC - VF) / 2.0 * (1.0 / AF + 1.0 / AC) -
                (DENL - DENG) * G * np.sin(ANG)) * 2.0 * AF / (SF * FF * DENL)

        if ABCD > 0:
            VFN = np.sqrt(ABCD)
        else:
            VFN = 0.95 * VF

        ABU = abs((VFN - VF) / VF)
        VF = (VFN + 9.0 * VF) / 10.0
        # print(VF, FI, ABU)

        if ABU < E1: break

        if ICON >= 2000:
            Ierr = 6
            break

    # Total pressure gradient due to friction
    PGF = -SF * FF * DENL * VF * abs(VF) / (2.0 * AXP)-SC * FC * DENG * VC * abs(VC) / (2.0 * AXP)
    PGG = -(DENL * HLF + DENC * (1.0 - HLF)) * G * np.sin(ANG)

    # total pressure gradient
    PGT = (PGF + PGG) / (1.0 - DENG * VC * VSG / (P * (1.0 - HLF)))

    # Total pressure gradient due to acceleration
    PGA = PGT - PGF - PGG

    # liquid holdup
    HL = HLF + HLC

    return HL, FE, FF, PGT, PGA, PGF, PGG, FGL, HLF, VF, SF, THF, ICON, IFGL


def ITFLOW_ETC(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL, P, PGT, PGF, PGG, PGA, FGL, IFGL, VF,
               VC, VT, HLF, SL, SI, THF, TH, TH0, FE, FF, FI, HLS, CU, CF, FQN, RSU, ICON, extra_in):

    # ----Initialize
    VM = VSL + VSG  # mixture velocity
    FGL = 'INT'
    IFGL = 4
    SL = D * PI

    HL, FF, PGT, PGA, PGF, PGG, FGL, HLF, CU, CS, CF, VF, VC, FQN, RSU, HLS, ICON, IFGL = \
        ITFLOW(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL)

    if FGL == 'D-B':
        SL = D * PI
        FF, HL, PGT, PGA, PGF, PGG = DBFLOW(D, ED, ANG, VSG, VSL, DENL, DENG, VISL)
        THF = D / 2.0
        TH = 1.0
        RSU = 0.0
        HLF = HL
        VF = VM
        VC = 0.0
        CF = 0.0
        CU = 0.0
        FE = 0.0
        FQN = 0.0
        SI = 0.0
        FI = 0.0
    elif FGL == 'STR':
        HL, FE, FF, PGT, PGA, PGF, PGG, FGL, HLF, VF, SF, THF, ICON, IFGL = \
            SAFLOW(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL, P)
        HLS = HL
        SL = SF

    return CU, CF, FE, FI, FQN, HLF, HLS, HL, PGT, PGA, PGF, PGG, RSU, SL, SI, TH, TH0, VF, VC, VT, FGL, IFGL, ICON


def GAL(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL, P, extra_in):
    """
    :param D: pipe diameter (m)
    :param ED: relative pipe wall roughness
    :param ANG: angle of pipe from hor (deg)
    :param VSL: liquid superficial velocity (m/s)
    :param VSG: gas superficial velocity (m/s)
    :param DENL: liquid density (kg/m3)
    :param DENG: gas density (kg/m3)
    :param VISL: liquid viscosity (Pa-s)
    :param VISG: gas viscosity (Pa-s)
    :param STGL: liquid surface tension (N/m)
    :param P: pressure (Pa)
    :param extra_in: additional input vector/list
    :return: FE, FI, FQN, HLF, HLS, HL, PGT, PGA, PGF, PGG, RSU, SL, SI, TH, TH0, VF, VC, VT, FGL, IFGL, ICON
    """
    global Ierr, AXP
    ANG = ANG * np.pi / 180.    # convert deg to rad

    FE, FI, FQN, HLF, HLS, PGT, PGA, PGF, PGG, RSU, SL, SI, TH, TH0, VF, VC, VT, FGL, \
    IFGL, ICON, THF, FF, CU, CF, HL = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '', 0, 0, 0, 0, 0, 0, 0

    # ----Initialize
    IFFM = extra_in[0]      # set the interfacial friction factor method
    Ierr = 0
    VM = VSL + VSG          # mixture velocity
    AXP = PI * D * D / 4.0  # cross sectional area of the pipe
    ICON = 0

    # - --------------------------
    # Check for single phase flow
    # - --------------------------

    ENS = VSL / VM
    HLS = ENS

    if ENS >= 0.99999:          # liquid
        FGL = 'LIQ'
        IFGL = 1
        HL = 1.0
        SL = D * PI
        FF, PGT, PGF, PGG, PGA = SGL(D, ED, ANG, P, DENL, VSL, VISL)
        THF = D / 2.0
        TH = 1.0

    elif ENS <= 0.0000001:      # gas
        FGL = 'GAS'
        IFGL = 7
        HL = 0.0
        SL = 0.0
        FF, PGT, PGF, PGG, PGA = SGL(D, ED, ANG, P, DENL, VSL, VISL)
        TH = 0.0
        THF = 0.0

    else:
        # - --------------------------
        # Check INT - D-B transition boundary
        # - --------------------------
        FGL = 'N-A'
        IFGL = 0
        if ENS > 0.36:
            VDB, Ierr = DISLUG(D, ED, ANG, VSG, DENL, DENG, VISL, STGL)
            if VSL > VDB:
                FGL = 'D-B'
                IFGL = 2
                SL = D * PI
                FF, HL, PGT, PGA, PGF, PGG = DBFLOW(D, ED, ANG, VSG, VSL, DENL, DENG, VISL)
                TH = 1.0
                THF = D / 2.0
                RSU = 0.0
                HLF = HL
                VF = VM
                VC = 0.0
                CF = 0.0
                CU = 0.0
                FE = 0.0
                FQN = 0.0
                SI = 0.0
                FI = 0.0

        if FGL != 'D-B':
            if ANG <= 0:        # downhill or horizontal
                # - --------------------------
                # Check I-SA transition boundary for downward flow (mostly I-S)
                # - --------------------------
                VST, Ierr = STSLUG(D, ED, ANG, VSG, DENL, DENG, VISL, VISG, STGL, IFFM)
                if VSL < VST:
                    HL, FE, FF, PGT, PGA, PGF, PGG, FGL, HLF, VF, SF, THF, ICON, IFGL = \
                        SAFLOW(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL, P)
                else:
                    CU, CF, FE, FI, FQN, HLF, HLS, HL, PGT, PGA, PGF, PGG, RSU, SL, SI, \
                    TH, TH0, VF, VC, VT, FGL, IFGL, ICON = \
                    ITFLOW_ETC(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL, P, PGT, PGF, PGG, PGA, FGL, IFGL, VF,
                               VC, VT, HLF, SL, SI, THF, TH, TH0, FE, FF, FI, HLS, CU, CF, FQN, RSU, ICON, extra_in)

            else:               # uphill
                VAN, Ierr = ANSLUG(D, ED, ANG, VSL, DENL, DENG, VISL, VISG, STGL, IFFM)
                if VSG > VAN:
                    FGL = 'ANN'
                    IFGL = 6
                    SL = D * PI
                    # print(FGL)
                    HL, FE, FF, PGT, PGA, PGF, PGG, FGL, HLF, VF, SF, THF, ICON, IFGL = \
                        SAFLOW(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL, P)
                    SL = SF
                    HLS = HL
                else:
                    # - --------------------------
                    # Check I-BU transition boundary
                    # - --------------------------
                    CKD = (DENL * DENL * G * D * D / (abs(DENL - DENG) * STGL)) ** 0.25
                    if CKD <= 4.36:
                        CU, CF, FE, FI, FQN, HLF, HLS, HL, PGT, PGA, PGF, PGG, RSU, SL, SI, \
                        TH, TH0, VF, VC, VT, FGL, IFGL, ICON = \
                            ITFLOW_ETC(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL, P, PGT, PGF, PGG, PGA, FGL,
                                       IFGL, VF, VC, VT, HLF, SL, SI, THF, TH, TH0, FE, FF,
                                       FI, HLS, CU, CF, FQN, RSU, ICON, extra_in)
                    else:
                        VBU = BUSLUG(D, ANG, VSL)
                        if (VSG < VBU) and (ANG > 60 * PI / 180):
                            FGL = 'BUB'
                            IFGL = 3
                            SL = D * PI
                            FM, HL, PGT, PGA, PGF, PGG = BUFLOW(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, STGL)
                            TH = 1.0
                            THF = D / 2.0
                            RSU = 0.0
                            HLF = HL
                            VF = VM
                            VC = 0.0
                            CF = 0.0
                            CU = 0.0
                            FE = 0.0
                            FQN = 0.0
                            SI = 0.0
                            FI = 0.0
                        else:
                            CU, CF, FE, FI, FQN, HLF, HLS, HL, PGT, PGA, PGF, PGG, RSU, SL, SI, \
                            TH, TH0, VF, VC, VT, FGL, IFGL, ICON = \
                                ITFLOW_ETC(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL, P,
                                           PGT, PGF, PGG, PGA, FGL, IFGL, VF, VC, VT, HLF, SL, SI,
                                           THF, TH, TH0, FE, FF, FI, HLS, CU, CF, FQN, RSU, ICON, extra_in)

    return CU, CF, FE, FI, FQN, HLF, HLS, HL, PGT, PGA, PGF, PGG, RSU, SL, SI, TH, TH0, VF, VC, VT, FGL, IFGL, ICON


if __name__ == "__main__":
    D = 0.1
    ED = 0
    ANG = 0
    VSL = 0.00078
    VSG = 0.21
    DENL = 850
    DENG = 18
    VISL = 0.001
    VISG = 0.000018
    STGL = 0.025
    P = 2075502
    extra_in = [1]

    CU, CF, FE, FI, FQN, HLF, HLS, HL, PGT, PGA, PGF, PGG, RSU, SL, SI, TH, TH0, VF, VC, VT, FGL, IFGL, ICON = \
        GAL(D, ED, ANG, VSL, VSG, DENL, DENG, VISL, VISG, STGL, P, extra_in)

    print(PGT, PGA, PGF, PGG, HL)

