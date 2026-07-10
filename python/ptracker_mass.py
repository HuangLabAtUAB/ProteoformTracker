"""Thin batch wrapper around pyteomics.mass, sourced once via reticulate.

Exists so building the reference-proteome mass index (tens of thousands of
sequences) pays one R<->Python round-trip for the whole batch, with the
per-sequence loop running in Python, instead of one round-trip per sequence.
"""

from pyteomics import mass


def batch_calculate_mass(sequences, average=False):
    return [mass.calculate_mass(sequence=s, average=average) for s in sequences]


def batch_ion_mass(sequences, ion_type, average=False):
    return [mass.calculate_mass(sequence=s, ion_type=ion_type, average=average) for s in sequences]
