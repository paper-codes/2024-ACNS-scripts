from typing import List

from measures.common import Code, SecurityLevels

# Levels at https://csrc.nist.gov/CSRC/media/Projects/post-quantum-cryptography/documents/round-3/seminars/oct-2020-gaj-kris-presentation.pdf, pag. 13
# n, k, t, name, level
CODES: List[Code] = [
    Code(
        "McEliece", SecurityLevels(1), 4096, 3556, 45, False
    ),  # Should be of category 1
    # FROM https://classic.mceliece.org/nist/mceliece-20171129.pdf, pag.11
    # FROM https://classic.mceliece.org/nist/mceliece-20201010.pdf, chap.3
    # k = n - mt, t; n = 2^m
    Code("McEliece", SecurityLevels(1), 3488, 2720, 64, False),  # Category 1, m = 12
    Code("McEliece", SecurityLevels(3), 4608, 3360, 96, False),  # Category 3, m = 13
    Code("McEliece", SecurityLevels(5), 6688, 5024, 128, False),  # Category 5, m = 13
    Code("McEliece", SecurityLevels(5), 6960, 5413, 119, False),  # Category 5, m = 13
    Code("McEliece", SecurityLevels(5), 8192, 6528, 128, False),  # Category 5, m = 13
    # from BIKE (message)
    # n = 2r, k = r
    Code("BIKE (message)", SecurityLevels(1), 24646, 12323, 134, True),  # Level 1
    Code("BIKE (message)", SecurityLevels(3), 49318, 24659, 199, True),  # Level 3
    Code("BIKE (message)", SecurityLevels(5), 81946, 40973, 264, True),  # Level 5
    # from BIKE (key)
    Code("BIKE (key)", SecurityLevels(1), 24646, 12323, 142, True),  # Level 1
    Code("BIKE (key)", SecurityLevels(3), 49318, 24659, 206, True),  # Level 3
    Code("BIKE (key)", SecurityLevels(5), 81946, 40973, 274, True),  # Level 5
    # from Esser_Bellin Syndrome Decoding Estimator
    Code("HQC", SecurityLevels(1), 35338, 17669, 132, True),  # Level 1
    Code("HQC", SecurityLevels(3), 71702, 35851, 200, True),  # Level 3
    Code("HQC", SecurityLevels(5), 115274, 57637, 262, True),  # Level 5
    # from https://csrc.nist.gov/CSRC/media/Projects/post-quantum-cryptography/documents/round-3/seminars/oct-2020-gaj-kris-presentation.pdf, pag.79;
    # each 1st row of each 2nd block
    # n = 2p, k = p, w = t
    # Code("leda", SecurityLevels(1), 28277 * 2, 28277, 129),  # Level 1
    # Code("leda", SecurityLevels(3), 52667 * 2, 52667, 195),  # Level 3
    # Code("leda", SecurityLevels(5), 83579 * 2, 83579, 260),  # Level 5
]
