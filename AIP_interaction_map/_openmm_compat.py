import importlib
import functools

from simtk import unit
from simtk.openmm import app


def _import_openmm_module(*module_names):
    last_error = None
    for module_name in module_names:
        try:
            return importlib.import_module(module_name)
        except ImportError as exc:
            last_error = exc
    raise last_error


_FORCEFIELD_MODULE = _import_openmm_module(
    "simtk.openmm.app.forcefield",
    "openmm.app.forcefield",
)
_TOPOLOGY_MODULE = _import_openmm_module(
    "simtk.openmm.app.topology",
    "openmm.app.topology",
)
compiled = _import_openmm_module(
    "simtk.openmm.app.internal.compiled",
    "openmm.app.internal.compiled",
)
_createResidueSignature = _FORCEFIELD_MODULE._createResidueSignature
Residue = _TOPOLOGY_MODULE.Residue


_ORIGINAL_ADD_EXTRA_PARTICLES = None


def _normalize_template_match(match):
    if match is None:
        return None

    try:
        template, matched_template_indices = match
    except (TypeError, ValueError):
        return None

    if template is None or matched_template_indices is None:
        return None

    return template, set(matched_template_indices)


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


def _bonded_to_atom_without_extra_particles(topology):
    bonded_to_atom = {}
    for atom in topology.atoms():
        bonded_to_atom[atom.index] = set()
    for atom1, atom2 in topology.bonds():
        if atom1.element is None or atom2.element is None:
            continue
        bonded_to_atom[atom1.index].add(atom2.index)
        bonded_to_atom[atom2.index].add(atom1.index)
    return bonded_to_atom


def _copy_template_without_extra_particles(template):
    template_without_extra_particles = app.ForceField._TemplateData(template.name)
    full_index_by_stripped_index = {}
    stripped_index_by_full_index = {}

    for full_index, atom in enumerate(template.atoms):
        if atom.element is None:
            continue

        stripped_index = len(template_without_extra_particles.atoms)
        full_index_by_stripped_index[stripped_index] = full_index
        stripped_index_by_full_index[full_index] = stripped_index

        stripped_atom = app.ForceField._TemplateAtomData(atom.name, atom.type, atom.element)
        stripped_atom.externalBonds = atom.externalBonds
        template_without_extra_particles.atoms.append(stripped_atom)

    for atom1, atom2 in template.bonds:
        if atom1 not in stripped_index_by_full_index or atom2 not in stripped_index_by_full_index:
            continue

        stripped_atom1 = stripped_index_by_full_index[atom1]
        stripped_atom2 = stripped_index_by_full_index[atom2]
        template_without_extra_particles.bonds.append((stripped_atom1, stripped_atom2))
        template_without_extra_particles.atoms[stripped_atom1].bondedTo.append(stripped_atom2)
        template_without_extra_particles.atoms[stripped_atom2].bondedTo.append(stripped_atom1)

    for atom_index in template.externalBonds:
        if atom_index in stripped_index_by_full_index:
            template_without_extra_particles.externalBonds.append(
                stripped_index_by_full_index[atom_index]
            )

    return template_without_extra_particles, full_index_by_stripped_index


def _build_templates_without_extra_particles(forcefield):
    templates_without_extra_particles = {}
    for template in getattr(forcefield, "_templates", {}).values():
        if not any(atom.element is None for atom in template.atoms):
            continue
        templates_without_extra_particles[template] = _copy_template_without_extra_particles(template)
    return templates_without_extra_particles


def _create_residue_without_extra_particles(residue):
    atoms_without_extra_particles = [atom for atom in residue.atoms() if atom.element is not None]
    if len(atoms_without_extra_particles) == len(list(residue.atoms())):
        return residue

    residue_without_extra_particles = Residue(
        residue.name,
        residue.index,
        residue.chain,
        residue.id,
        residue.insertionCode,
    )
    residue_without_extra_particles._atoms = atoms_without_extra_particles
    return residue_without_extra_particles


def _match_residue_to_template(residue, template, bonded_to_atom, ignore_external_bonds):
    try:
        return compiled.matchResidueToTemplate(
            residue,
            template,
            bonded_to_atom,
            ignore_external_bonds,
        )
    except TypeError:
        return compiled.matchResidueToTemplate(
            residue,
            template,
            bonded_to_atom,
            ignore_external_bonds,
            False,
        )


def _get_template_match_without_extra_particles(
    forcefield,
    residue,
    bonded_to_atom_without_extra_particles,
    ignore_external_bonds,
    templates_without_extra_particles,
):
    signature = _createResidueSignature([atom.element for atom in residue.atoms()])
    residue_without_extra_particles = _create_residue_without_extra_particles(residue)

    for template in getattr(forcefield, "_templateSignatures", {}).get(signature, ()):
        if template not in templates_without_extra_particles:
            continue

        template_without_extra_particles, full_index_by_stripped_index = templates_without_extra_particles[template]
        matches = _match_residue_to_template(
            residue_without_extra_particles,
            template_without_extra_particles,
            bonded_to_atom_without_extra_particles,
            ignore_external_bonds,
        )
        if matches is None:
            continue

        return template, {
            full_index_by_stripped_index[match_index]
            for match_index in matches
        }

    return None


def _get_template_match(forcefield, residue, bonded_to_atom, ignore_external_bonds):
    template_signatures = getattr(forcefield, "_templateSignatures", None)
    try:
        return _normalize_template_match(forcefield._getResidueTemplateMatches(
            residue,
            bonded_to_atom,
            template_signatures,
            ignoreExternalBonds=ignore_external_bonds,
            ignoreExtraParticles=True,
        ))
    except TypeError as exc:
        if "ignoreExtraParticles" not in str(exc) and "ignoreExternalBonds" not in str(exc):
            raise
    except Exception:
        return None

    try:
        return _normalize_template_match(forcefield._getResidueTemplateMatches(
            residue,
            bonded_to_atom,
            template_signatures,
            ignoreExternalBonds=ignore_external_bonds,
        ))
    except TypeError as exc:
        if "ignoreExternalBonds" not in str(exc):
            raise
    except Exception:
        return None

    try:
        return _normalize_template_match(forcefield._getResidueTemplateMatches(
            residue,
            bonded_to_atom,
            template_signatures,
            ignore_external_bonds,
            True,
        ))
    except TypeError:
        pass
    except Exception:
        return None

    try:
        return _normalize_template_match(forcefield._getResidueTemplateMatches(
            residue,
            bonded_to_atom,
            template_signatures,
            ignore_external_bonds,
        ))
    except TypeError:
        pass
    except Exception:
        return None

    try:
        return _normalize_template_match(forcefield._getResidueTemplateMatches(
            residue,
            bonded_to_atom,
            template_signatures,
        ))
    except Exception:
        return None


def _fallback_add_extra_particles(modeller, forcefield, ignore_external_bonds):
    topology = modeller.topology
    positions = modeller.positions
    bonded_to_atom = _bonded_to_atom(topology)
    bonded_to_atom_without_extra_particles = _bonded_to_atom_without_extra_particles(topology)
    templates_without_extra_particles = _build_templates_without_extra_particles(forcefield)

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
            if match is None:
                match = _get_template_match_without_extra_particles(
                    forcefield,
                    residue,
                    bonded_to_atom_without_extra_particles,
                    ignore_external_bonds,
                    templates_without_extra_particles,
                )
            template = None
            matched_template_indices = set()
            if match is not None:
                template, matched_template_indices = match

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