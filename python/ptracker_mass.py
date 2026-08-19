"""Thin batch wrapper around pyteomics.mass, sourced once via reticulate.

Exists so building the reference-proteome mass index (tens of thousands of
sequences) pays one R<->Python round-trip for the whole batch, with the
per-sequence loop running in Python, instead of one round-trip per sequence.
"""

from pyteomics import mass

try:
    import IsoSpecPy
except ImportError:
    IsoSpecPy = None


def batch_calculate_mass(sequences, average=False):
    return [mass.calculate_mass(sequence=s, average=average) for s in sequences]


def batch_ion_mass(sequences, ion_type, average=False):
    return [mass.calculate_mass(sequence=s, ion_type=ion_type, average=average) for s in sequences]


def isotope_pattern(sequence, prob_to_cover=0.9999, bin_width=1.0):
    """Exact-composition isotope pattern for a bare amino-acid sequence,
    binned to nominal-mass peaks.

    Deliberately NOT using the Averagine model (an idealized "average
    residue" composition scaled to a target mass) -- Averagine exists for
    the case where the sequence is unknown (e.g. deconvoluting a raw
    intact-mass measurement). Every caller of this function already has
    the real sequence, so the real elemental composition
    (pyteomics.mass.Composition, summing each residue's own formula) is
    used directly instead -- strictly more accurate, same cost.

    IsoSpecPy (Lacki et al. 2017, Anal Chem 89:3346) enumerates the exact
    fine-structure isotopologue distribution (every combination of which
    specific atoms are a heavy isotope) up to a total-probability
    coverage threshold, which is fast even at intact-protein scale
    (milliseconds for a >10 kDa protein) -- but that fine structure is
    finer than any real spectrum ever shows: an instrument can't tell a
    +1 Da shift caused by one heavy carbon apart from one caused by one
    heavy nitrogen, so those fine-structure peaks are indistinguishable
    in practice and must be summed into one binned nominal-mass peak
    before comparing to (or plotting like) an observed spectrum. The bin
    width defaults to 1.0 Da, matching real isotope spacing.

    @param sequence bare amino acid sequence (standard 20 AA letters, no
        PTMs -- a PTM's mass delta is applied by the R caller as a flat
        shift to every returned peak's mass afterwards, consistent with
        how PTMs are represented everywhere else in this codebase: a
        mass delta only, not a full elemental-formula change)
    @param prob_to_cover total probability mass to enumerate before
        stopping (the excluded tail is genuinely negligible at the
        default 99.99%)
    @param bin_width nominal-mass bin width in Da (1.0 = real isotope
        spacing; matches every adjacent isotope peak)
    @return dict(masses=[...], probs=[...]) ascending by mass, probs
        renormalized to sum to 1 over the enumerated coverage
    """
    if IsoSpecPy is None:
        raise ImportError(
            "IsoSpecPy is not installed -- run: "
            "Rscript scripts/setup_python_env.R"
        )
    comp = mass.Composition(sequence=sequence)
    formula = "".join(f"{el}{int(n)}" for el, n in comp.items())
    iso = IsoSpecPy.IsoTotalProb(formula=formula, prob_to_cover=prob_to_cover)
    binned = iso.binned(bin_width)
    pairs = sorted(zip(binned.masses, binned.probs))
    total = sum(p for _, p in pairs)
    return {
        "masses": [m for m, _ in pairs],
        "probs": [p / total for _, p in pairs],
    }
