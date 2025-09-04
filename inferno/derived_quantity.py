import pynbody
import numpy as np
from . import imf
from . import asplund

asplund = asplund.Asplund()

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def logAge(snap):
    return np.log10(snap.st['age'].in_units('yr'))

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def metal_H(snap, Y=0.246):
    Z = 2.09*snap.st['metal_O'] + 1.06*snap.st['metal_Fe']
    return (1.0-Y)*(1.0-Z)

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def mform(snap):
    code_units_to_Msol = 4921203503291205.0
    value = pynbody.array.SimArray(snap.st['birth_mass'] * code_units_to_Msol,
                                    units='Msol')
    return value

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def mform_resampling(snap,seed=920926):
    value = pynbody.array.SimArray(np.zeros(len(snap.st)),units='Msol')
    np.random.seed(seed) # Ensure same sampling everytime for reproducability
    
    star = (snap.st['tag'] > 1)
    value[star] = imf.resample_stellar_mass(snap.st['mform'].in_units('Msol')[star])
     
    spop = (snap.st['tag'] == 1)
    value[spop] = snap.st['mform'][spop]
    
    return value

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def sampling_mass(snap):
    code_units_to_Msol = 4921203503291205.0
    value = pynbody.array.SimArray(snap.st['extra_pvar'] * code_units_to_Msol,
                                    units='Msol')
    return value

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def metal_Z(snap):
    Z = 2.09*snap.st['metal_O'] + 1.06*snap.st['metal_Fe']
    return Z

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def metal_MH(snap):
    Z = 2.09*snap.st['metal_O'] + 1.06*snap.st['metal_Fe']
    MH = np.log10(Z/0.0207)
    return MH

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def FeH(snap):
    Fe = snap.st['metal_Fe']
    H = snap.st['metal_H']
    FeH = asplund.relative(Fe, H, 'Fe', 'H')
    return FeH

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def CFe(snap):
    Fe = snap.st['metal_Fe']
    C = snap.st['metal_C']
    CFe = asplund.relative(C, Fe, 'C', 'Fe')
    return CFe

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def NFe(snap):
    Fe = snap.st['metal_Fe']
    N = snap.st['metal_N']
    NFe = asplund.relative(N, Fe, 'N', 'Fe')
    return NFe

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def OFe(snap):
    Fe = snap.st['metal_Fe']
    O = snap.st['metal_O']
    OFe = asplund.relative(O, Fe, 'O', 'Fe')
    return OFe

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def MgFe(snap):
    Fe = snap.st['metal_Fe']
    Mg = snap.st['metal_Mg']
    MgFe = asplund.relative(Mg, Fe, 'Mg', 'Fe')
    return MgFe

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def AlFe(snap):
    Fe = snap.st['metal_Fe']
    Al = snap.st['metal_Al']
    AlFe = asplund.relative(Al, Fe, 'Al', 'Fe')
    return AlFe

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def SiFe(snap):
    Fe = snap.st['metal_Fe']
    Si = snap.st['metal_Si']
    SiFe = asplund.relative(Si, Fe, 'Si', 'Fe')
    return SiFe

@pynbody.snapshot.ramses.RamsesSnap.derived_quantity
def EuFe(snap):
    Fe = snap.st['metal_Fe']
    Eu = snap.st['metal_Eu']
    EuFe = asplund.relative(Eu, Fe, 'Eu', 'Fe')
    return EuFe
