# from Bio import
import math

NO_DEC_PLACES = 3

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    # ignore sites with indels and Ns
    valid_pairs = [(a, b) for a, b in zip(reference_sequence, distant_sequence) if a != "-" and b != "-" and a != 'N' and b != 'N']

    ham_dist = sum(a != b for a, b in valid_pairs)
    no_valid_sites = len(valid_pairs)

    # p = no. mismatches (aka Hamming distance) / total valid alignment length
    p_val = ham_dist / no_valid_sites

    if p_val >= 0.75:
        return 0.0

    # correction - 3/4 * ln(1-4/3*p) [per site]
    jc_correction = -3/4 * math.log(1 - (4/3) * p_val) # mutation rate per nucleotide per day
    corrected_dist = jc_correction * no_valid_sites

    return round(corrected_dist, NO_DEC_PLACES)


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()

if __name__ == "__main__":
    # (e.g. AGGT-GA)
    # distant_sequence: str
    #     A string of nucleotides in a sequence after the alignment
    #     with a reference (e.g. AGC-AGA)
    print(jukes_cantor("AGGT-GA", "AGC-AGA"))