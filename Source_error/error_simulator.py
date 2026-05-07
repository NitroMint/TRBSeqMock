import numpy as np
try:
    from numba import jit, uint8, int32 # type: ignore
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False


BASE_ENCODE = {"A": 0, "T": 1, "C": 2, "G": 3}
BASE_DECODE = ["A", "T", "C", "G"]
REPLACE_MAP = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]], dtype=np.uint8)


def _batch_add_seq_error_py(seq_encoded: np.ndarray, replace_mask: np.ndarray) -> np.ndarray:
    """Pure-Python fallback for when numba is unavailable."""
    import random
    seq_len = len(seq_encoded)
    for i in range(seq_len):
        if replace_mask[i]:
            current_base = int(seq_encoded[i])
            replace_idx = random.randint(0, 2)
            seq_encoded[i] = REPLACE_MAP[current_base, replace_idx]
    return seq_encoded


if _HAS_NUMBA:
    @jit(nopython=True, cache=True)
    def _batch_add_seq_error_jit(seq_encoded: np.ndarray, replace_mask: np.ndarray) -> np.ndarray:
        seq_len = len(seq_encoded)
        for i in range(seq_len):
            if replace_mask[i]:
                current_base = seq_encoded[i]
                replace_idx = np.random.randint(0, 3)
                seq_encoded[i] = REPLACE_MAP[current_base, replace_idx]
        return seq_encoded

    _batch_add_seq_error = _batch_add_seq_error_jit
else:
    _batch_add_seq_error = _batch_add_seq_error_py


def import_error(seq, p: float = 0.001) -> str:
    encoded = np.array([BASE_ENCODE[b] for b in seq], dtype=np.uint8)
    replace_mask = np.random.random(len(encoded)) < p
    error_encoded = _batch_add_seq_error(encoded.copy(), replace_mask)
    return "".join([BASE_DECODE[b] for b in error_encoded])