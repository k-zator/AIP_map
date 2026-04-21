import functools

from simtk import unit
from simtk.openmm import app


_ORIGINAL_ADD_EXTRA_PARTICLES = None


def _call_original_add_extra_particles(modeller, forcefield, ignore_external_bonds):
    try:
        return _ORIGINAL_ADD_EXTRA_PARTICLES(
            modeller,
            forcefield,
            ignoreExternalBonds=ignore_external_bonds,
        )
    except TypeError as exc:
        if "ignoreExternalBonds" not in str(exc):
            raise

    try:
        return _ORIGINAL_ADD_EXTRA_PARTICLES(
            modeller,
            forcefield,
            ignore_external_bonds,
        )
    except TypeError as exc:
        if not ignore_external_bonds:
            raise
        if "positional" not in str(exc) and "takes" not in str(exc):
            raise

    return _ORIGINAL_ADD_EXTRA_PARTICLES(modeller, forcefield)


def _bonded_to_atom(topology):
    bonded_to_atom = {}
    for atom in topology.atoms():
        bonded_to_atom[atom.index] = set()
    for atom1, atom2 in topology.bonds():
        bonded_to_atom[atom1.index].add(atom2.index)
        bonded_to_atom[atom2.index].add(atom1.index)
    return bonded_to_atom


def _get_template_match(forcefield, residue, bonded_to_atom, ignore_external_bonds):
    template_signatures = getattr(forcefield, "_templateSignatures", None)
    try:
        return forcefield._getResidueTemplateMatches(
            residue,
            bonded_to_atom,
            template_signatures,
            ignoreExternalBonds=ignore_external_bonds,
            ignoreExtraParticles=True,
        )
    except TypeError as exc:
        if "ignoreExtraParticles" not in str(exc) and "ignoreExternalBonds" not in str(exc):
            raise
    except Exception:
        return None

    try:
        return forcefield._getResidueTemplateMatches(
            residue,
            bonded_to_atom,
            template_signatures,
            ignoreExternalBonds=ignore_external_bonds,
        )
    except TypeError as exc:
        if "ignoreExternalBonds" not in str(exc):
            raise
    except Exception:
        return None

    try:
        return forcefield._getResidueTemplateMatches(
            residue,
            bonded_to_atom,
            template_signatures,
            ignore_external_bonds,
            True,
        )
    except TypeError:
        pass
    except Exception:
        return None

    try:
        return forcefield._getResidueTemplateMatches(
            residue,
            bonded_to_atom,
            template_signatures,
            ignore_external_bonds,
        )
    except TypeError:
        pass
    except Exception:
        return None

    try:
        return forcefield._getResidueTemplateMatches(
            residue,
            bonded_to_atom,
            template_signatures,
        )
    except Exception:
        return None


def _fallback_add_extra_particles(modeller, forcefield, ignore_external_bonds):
    topology = modeller.topology
    positions = modeller.positions
    bonded_to_atom = _bonded_to_atom(topology)

    new_topology = app.Topology()
    box_vectors = topology.getPeriodicBoxVectors()
    if box_vectors is not None:
        new_topology.setPeriodicBoxVectors(box_vectors)

    atom_map = {}
    new_positions = []
    position_unit = positions.unit
    added_particle = False

    for chain in topology.chains():
        new_chain = new_topology.addChain(chain.id)
        for residue in chain.residues():
            match = _get_template_match(forcefield, residue, bonded_to_atom, ignore_external_bonds)
            template = None
            matched_template_indices = set()
            if match is not None:
                template, matched_template_indices = match
                matched_template_indices = set(matched_template_indices)

            new_residue = new_topology.addResidue(
                residue.name,
                new_chain,
                residue.id,
                residue.insertionCode,
            )

            anchor_position = None
            for atom in residue.atoms():
                new_atom = new_topology.addAtom(atom.name, atom.element, new_residue, atom.id)
                atom_map[atom] = new_atom
                atom_position = positions[atom.index]
                new_positions.append(atom_position)
                if anchor_position is None:
                    anchor_position = atom_position

            if template is None:
                continue

            missing_template_indices = [
                index
                for index in range(len(template.atoms))
                if index not in matched_template_indices
            ]
            if any(template.atoms[index].element is not None for index in missing_template_indices):
                return False

            for index in missing_template_indices:
                template_atom = template.atoms[index]
                new_topology.addAtom(template_atom.name, None, new_residue)
                new_positions.append(anchor_position)
                added_particle = True

    if not added_particle:
        return False

    for atom1, atom2 in topology.bonds():
        new_topology.addBond(atom_map[atom1], atom_map[atom2])

    modeller.topology = new_topology
    modeller.positions = unit.Quantity(new_positions, position_unit)
    return True


def apply_openmm_compat_patches():
    global _ORIGINAL_ADD_EXTRA_PARTICLES

    if getattr(app.Modeller.addExtraParticles, "__module__", "") == __name__:
        return

    _ORIGINAL_ADD_EXTRA_PARTICLES = app.Modeller.addExtraParticles

    @functools.wraps(_ORIGINAL_ADD_EXTRA_PARTICLES)
    def _compat_add_extra_particles(self, forcefield, ignoreExternalBonds=False):
        try:
            return _call_original_add_extra_particles(
                self,
                forcefield,
                ignoreExternalBonds,
            )
        except ValueError as exc:
            if not _fallback_add_extra_particles(self, forcefield, ignoreExternalBonds):
                raise exc
            return None

    app.Modeller.addExtraParticles = _compat_add_extra_particles