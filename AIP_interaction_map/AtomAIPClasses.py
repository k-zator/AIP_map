

class Atom():
    "Describes Atom and its AIPs for the pairing algorithm"

    def __init__(self, index, atom_type, xyz):
        self.index = index
        self.type = atom_type
        self.xyz = xyz
        if self.type[0] in ["H", "O", "N", "S"] \
                and "no_lp" not in self.type \
                and "soft" not in self.type \
                and "pl3" not in self.type \
                and ".O2" not in self.type:
            self.polar = True
        else:
            self.polar = False

    def set_sasa35(self, value):
        self.sasa35 = value

    def set_sasa70(self, value):
        self.sasa70 = value

    def set_sasa35f(self, value):
        self.sasa35f = value

    def set_sasa70f(self, value):
        self.sasa70f = value

    def set_sa(self, value):
        self.sa = value

class AIP():
    "Describes AIP and its atom owner for the pairing algoritm"

    def __init__(self, index, value, atom_owner, atom_type, xyz, isosurface, fraction, dual=False):
        self.index = index
        self.value = value
        self.atom = atom_owner
        self.type = atom_type
        self.xyz = xyz
        self.isosurface = isosurface
        self.fraction = fraction
        if self.isosurface == 0.0300 or self.type in ["H.N", "H.O", "H.N.group2"]:
            self.polar = True
        else:
            self.polar = False

        if dual > 0:
            self.dual = True
            self.valueMax = dual
        else:
            self.dual = False
            
    def set_sasa_b35(self, value):
        self.sasa_b35 = value

    def set_sasa_b70(self, value):
        self.sasa_b70 = value

    def set_sasa_f(self, value):
        self.sasa_f = value
