import numba as nb
import numpy as np


@nb.njit(cache=True, fastmath=True)
def compute_Pi_L_numba(
    out: np.ndarray,
    aa_idx_seq: np.ndarray,
    Pixy: np.ndarray,
    Pins: np.ndarray,
    fnb: np.ndarray,
    znb: np.ndarray,
    S: np.ndarray,
    D: np.ndarray,
    T: np.ndarray,
    lD: np.ndarray,
    lT: np.ndarray,
    max_align: int
) -> None:
    
    max_insertions = Pins.size - 1
    M_codons = max_insertions // 3

    v0 = np.empty(4)
    cur = np.empty(4)

    _fill_from_left(
        out, aa_idx_seq, Pixy, Pins, fnb, znb,
        S, D, T, lD, lT, max_align, v0, cur,
        M_codons
    )


@nb.njit(cache=True, fastmath=True)
def compute_Pi_V_insVJ_given_J_numba(
    out: np.ndarray,
    aa_idx_seq: np.ndarray,
    Pixy: np.ndarray,
    Pins: np.ndarray,
    fnb: np.ndarray,
    znb: np.ndarray,
    S: np.ndarray,
    D: np.ndarray,
    T: np.ndarray,
    lD: np.ndarray,
    lT: np.ndarray,
    max_align: int
) -> None:

    numD = Pixy.shape[0]
    max_insertions = Pins.size - 1
    M_codons = max_insertions // 3

    v0 = np.empty(4)
    cur = np.empty(4)

    for d in range(numD):
        O = out[d]
        Pij = Pixy[d]

        _fill_from_left(
            O, aa_idx_seq, Pij, Pins, fnb, znb, 
            S, D, T, lD, lT, max_align, v0, cur, 
            M_codons
        )


@nb.njit(fastmath=True)
def _fill_from_left(
    O: np.ndarray,
    aa_idx_seq: np.ndarray,
    Pij: np.ndarray,
    Pins: np.ndarray,
    fnb: np.ndarray,
    znb: np.ndarray,
    S: np.ndarray,
    D: np.ndarray,
    T: np.ndarray,
    lD: np.ndarray,
    lT: np.ndarray,
    max_align: int,
    v0: np.ndarray,
    cur: np.ndarray,
    M_codons: int
) -> None:

    # Frame 1
    j0 = 0
    for init_pos in range(0, max_align, 3):
        _load4(Pij, v0, init_pos)
        _axpy4(Pins[0], v0, O, init_pos)

        aa0 = aa_idx_seq[j0]
        _dot4_axpy4(Pins[1], lD[aa0], v0, O, init_pos + 1)

        _dot4(lT[aa0], v0, cur)
        s = _sum4(cur)
        O[0, init_pos + 2] += Pins[2] * s

        base_ins = 2
        _walk_from_left(
            base_ins, j0, M_codons, aa_idx_seq,
            O, init_pos, Pins, S, D, T, cur
        )
        j0 += 1

    # Frame 2
    j0 = 0
    for init_pos in range(1, max_align, 3):
        _load4(Pij, v0, init_pos)
        _axpy4(Pins[0], v0, O, init_pos)
    
        for k in range(4):
            cur[k] = v0[k] * fnb[k]

        s = _sum4(cur)
        O[0, init_pos + 1] += Pins[1] * s

        base_ins = 1
        _walk_from_left(
            base_ins, j0, M_codons, aa_idx_seq,
            O, init_pos, Pins, S, D, T, cur
        )
        j0 += 1

    # Frame 3
    j0 = 0
    for init_pos in range(2, max_align, 3):
        a0 = Pij[0, init_pos]
        O[0, init_pos] += Pins[0] * a0

        for k in range(4):
            cur[k] = znb[k] * a0

        base_ins = 0
        _walk_from_left(
            base_ins, j0, M_codons, aa_idx_seq,
            O, init_pos, Pins, S, D, T, cur
        )
        j0 += 1


@nb.njit(cache=True, fastmath=True)
def compute_Pi_JinsDJ_given_D_numba(
    out: np.ndarray,
    aa_idx_seq: np.ndarray,
    Pixy: np.ndarray,
    Pins: np.ndarray,
    fnb: np.ndarray,
    znb: np.ndarray,
    S: np.ndarray,
    D: np.ndarray,
    T: np.ndarray,
    rD: np.ndarray,
    rT: np.ndarray,
    max_align: int
) -> None:

    L = aa_idx_seq.shape[0]
    numD = Pixy.shape[0]
    max_insertions = Pins.size - 1
    M_codons = max_insertions // 3

    v0 = np.empty(4)
    cur = np.empty(4)

    for d in range(numD):
        O = out[d]
        Pij = Pixy[d]

        # Frame 1
        j0 = L - 1
        for init_pos in range(-1, -(max_align + 1), -3):
            _load4(Pij, v0, init_pos)
            _axpy4(Pins[0], v0, O, init_pos)

            aa0 = aa_idx_seq[j0]
            _dot4_axpy4(Pins[1], rD[aa0], v0, O, init_pos - 1)
            _dot4(rT[aa0], v0, cur)

            s = _sum4(cur)
            O[0, init_pos - 2] += Pins[2] * s
            base_ins = 2
            _walk_from_right(base_ins, j0, M_codons, aa_idx_seq, O, init_pos, Pins, S, D, T, cur)
            j0 -= 1

        # Frame 2
        j0 = L - 1
        for init_pos in range(-2, -(max_align + 1), -3):
            _load4(Pij, v0, init_pos)
            _axpy4(Pins[0], v0, O, init_pos)

            for k in range(4):
                cur[k] = v0[k] * fnb[k]

            s = _sum4(cur)
            O[0, init_pos - 1] += Pins[1] * s
            base_ins = 1
            _walk_from_right(base_ins, j0, M_codons, aa_idx_seq, O, init_pos, Pins, S, D, T, cur)
            j0 -= 1

        # Frame 3
        j0 = L - 1
        for init_pos in range(-3, -(max_align + 1), -3):
            a0 = Pij[0, init_pos]
            O[0, init_pos] += Pins[0] * a0

            for k in range(4):
                cur[k] = znb[k] * a0

            base_ins = 0
            _walk_from_right(base_ins, j0, M_codons, aa_idx_seq, O, init_pos, Pins, S, D, T, cur)
            j0 -= 1


@nb.njit(fastmath=True)
def _walk_from_right(
    base_ins: int,
    j0: int,
    M_codons: int,
    aa_idx_seq: np.ndarray,
    O: np.ndarray,
    init_pos: int,
    Pins: np.ndarray,
    S: np.ndarray,
    D: np.ndarray,
    T: np.ndarray,
    cur: np.ndarray
) -> None:
    start_j = j0 - 1
    end_j_inc = max(0, j0 - M_codons + 1)
    for j in range(start_j, end_j_inc - 1, -1):
        aidx = aa_idx_seq[j]

        p0 = base_ins + 1
        p = init_pos - base_ins - 1
        _dot4_axpy4(Pins[p0], S[aidx], cur, O, p)

        p0 += 1
        p = init_pos - base_ins - 2
        _dot4_axpy4(Pins[p0], D[aidx], cur, O, p)

        _dot4_inplace(T[aidx], cur)

        p = init_pos - base_ins - 3
        s = _sum4(cur)
        O[0, p] += Pins[base_ins + 3] * s
        base_ins += 3


@nb.njit(fastmath=True)
def _walk_from_left(
    base_ins: int,
    j0: int,
    M_codons: int,
    aa_idx_seq: np.ndarray,
    O: np.ndarray,
    init_pos: int,
    Pins: np.ndarray,
    S: np.ndarray,
    D: np.ndarray,
    T: np.ndarray,
    cur: np.ndarray
) -> None:
    L = aa_idx_seq.shape[0]
    start_j = j0 + 1
    end_j = min(L, j0 + M_codons)

    for j in range(start_j, end_j):
        aidx = aa_idx_seq[j]

        p0 = base_ins + 1
        p = init_pos + p0
        _dot4_axpy4(Pins[p0], S[aidx], cur, O, p)

        p0 += 1
        p = init_pos + p0
        _dot4_axpy4(Pins[p0], D[aidx], cur, O, p)

        _dot4_inplace(T[aidx], cur)
        p = init_pos + base_ins + 3
        s = _sum4(cur)
        O[0, p] += Pins[base_ins + 3] * s
        base_ins += 3


@nb.njit(cache=True, fastmath=True)
def compute_Pi_R_one_numba(
    Pi_R: np.ndarray,
    L3: int,
    min_pos: int,
    dseq: np.ndarray,
    Pij: np.ndarray,
    Pdel_D: np.ndarray,
    max_del_D: np.ndarray,
    min_del_D: np.ndarray,
    PD_nt: np.ndarray,
    PD_2nd_stack: np.ndarray,
    zeroD: float,
    aa_idx_seq: np.ndarray,
    allow_lsf_stack,
    allow1: np.ndarray,
    allow2: np.ndarray,
    allow3: np.ndarray
) -> None:
    Ld = dseq.shape[0]
    num_delr = Pdel_D.shape[1]

    base_pair_prob = np.empty((4, 4))
    v0 = np.empty(4)

    # Frame 1
    for init_pos in range(-1, -L3 - 1, -3):
        _load4(Pij, v0, init_pos)
        for k in range(4):
            Pi_R[k, init_pos] += zeroD * v0[k]

        aa0 = aa_idx_seq[(init_pos + L3) // 3]
        lsf = allow_lsf_stack[aa0]

        for l in range(4):
            for s in range(4):
                acc = 0.0
                for t in range(4):
                    if lsf[l, s, t]:
                        acc += v0[t]
                base_pair_prob[l, s] = acc

        for delDr in range(num_delr):
            if min_del_D[delDr] == -1:
                continue
            
            dpos1 = Ld - delDr - 1
            if dpos1 <= max_del_D[delDr]:
                second_idx = dseq[dpos1]
                w = Pdel_D[dpos1, delDr]
                Pi_R[0, init_pos - 1] += w * base_pair_prob[0, second_idx]
                Pi_R[1, init_pos - 1] += w * base_pair_prob[1, second_idx]
                Pi_R[2, init_pos - 1] += w * base_pair_prob[2, second_idx]
                Pi_R[3, init_pos - 1] += w * base_pair_prob[3, second_idx]

            dpos2 = Ld - delDr - 2
            
            last_idx = dseq[dpos2]
            second_idx = dseq[dpos2 + 1]
            base_prob = base_pair_prob[last_idx, second_idx]
            if base_prob == 0.0:
                continue

            if dpos2 <= max_del_D[delDr]:
                Pi_R[0, init_pos - 2] += Pdel_D[dpos2, delDr] * base_prob

            stop = max(init_pos - Ld + delDr, min_pos)
            for pos in range(init_pos - 3, stop - 1, -1):
                D_pos = Ld - delDr - 1 - ((init_pos - 1) - pos)
                ok = _update_Pi_R(
                    D_pos, max_del_D, delDr, Pdel_D, aa_idx_seq,
                    pos, dseq, allow1, allow2, allow3, Pi_R, PD_nt,
                    PD_2nd_stack, base_prob
                )
                if not ok:
                    break
  
    # Frame 2
    for init_pos in range(-2, -L3 - 1, -3):
        _load4(Pij, v0, init_pos)
        for k in range(4):
            Pi_R[k, init_pos] += zeroD * v0[k]

        p_last = v0

        for delDr in range(num_delr):
            if min_del_D[delDr] == -1:
                continue

            dpos1 = Ld - delDr - 1
            last_idx = dseq[dpos1]
            base_prob = p_last[last_idx]
            if base_prob == 0.0:
                continue

            if dpos1 <= max_del_D[delDr]:
                Pi_R[0, init_pos - 1] += Pdel_D[dpos1, delDr] * base_prob

            stop = max(init_pos - Ld + delDr, min_pos)
            for pos in range(init_pos - 2, stop - 1, -1):
                D_pos = dpos1 - ((init_pos - 1) - pos)  
                ok = _update_Pi_R(
                    D_pos, max_del_D, delDr, Pdel_D, aa_idx_seq,
                    pos, dseq, allow1, allow2, allow3, Pi_R, PD_nt,
                    PD_2nd_stack, base_prob
                )
                if not ok:
                    break

    # Frame 3
    for init_pos in range(-3, -L3 - 1, -3):
        pA = Pij[0, init_pos]
        Pi_R[0, init_pos] += zeroD * pA
 
        base_prob = pA
        for delDr in range(num_delr):
            if min_del_D[delDr] == -1:
                continue

            stop = max(init_pos - Ld + delDr, min_pos)
            for pos in range(init_pos - 1, stop - 1, -1):
                D_pos = Ld - delDr - 1 - ((init_pos - 1) - pos)
                ok = _update_Pi_R(
                    D_pos, max_del_D, delDr, Pdel_D, aa_idx_seq,
                    pos, dseq, allow1, allow2, allow3, Pi_R, PD_nt,
                    PD_2nd_stack, base_prob
                )
                if not ok:
                    break


@nb.njit(fastmath=True)
def _update_Pi_R(
    D_pos: int,
    max_del_D: np.ndarray,
    delDr: int,
    Pdel_D: np.ndarray,
    aa_idx_seq: np.ndarray,
    pos: int,
    dseq: np.ndarray,
    allow1: np.ndarray,
    allow2: np.ndarray,
    allow3: np.ndarray,
    Pi_R: np.ndarray,
    PD_nt: np.ndarray,
    PD_2nd_stack: np.ndarray,
    base_prob: float
) -> bool:

    if D_pos > max_del_D[delDr]:
        current_PdelDldelDr = 0
    else:
        current_PdelDldelDr = Pdel_D[D_pos, delDr]

    aa = aa_idx_seq[pos // 3]
    r = pos % 3
    if r == 2:
        single = dseq[D_pos]
        if allow1[aa, single] != 0:
            w = current_PdelDldelDr * base_prob
            Pi_R[0, pos] += w * PD_nt[0, D_pos]
            Pi_R[1, pos] += w * PD_nt[1, D_pos]
            Pi_R[2, pos] += w * PD_nt[2, D_pos]
            Pi_R[3, pos] += w * PD_nt[3, D_pos]
            return True
    elif r == 1:
        pair = 4*dseq[D_pos] + dseq[D_pos + 1]
        if allow2[aa, pair] != 0:
            w = current_PdelDldelDr * base_prob
            Pi_R[0, pos] += w * PD_2nd_stack[aa, 0, D_pos]
            Pi_R[1, pos] += w * PD_2nd_stack[aa, 1, D_pos]
            Pi_R[2, pos] += w * PD_2nd_stack[aa, 2, D_pos]
            Pi_R[3, pos] += w * PD_2nd_stack[aa, 3, D_pos]
            return True
    else:
        triplet = 16*dseq[D_pos] + 4*dseq[D_pos + 1] + dseq[D_pos + 2]
        if allow3[aa, triplet] != 0:
            Pi_R[0, pos] += current_PdelDldelDr * base_prob
            return True

    return False


@nb.njit(inline='always')
def _load4(src: np.ndarray, dst: np.ndarray, init_pos: int) -> None:
    """ Load a integer into a 4-length array"""
    dst[0] = src[0, init_pos]
    dst[1] = src[1, init_pos]
    dst[2] = src[2, init_pos]
    dst[3] = src[3, init_pos]


@nb.njit(inline='always')
def _dot4(M: np.ndarray, x: np.ndarray, out: np.ndarray) -> None:
    """ Dot product of two 4-length arrays to avoid temp array
    allocation. Saves to existing third array `out`.
    """

    out[0] = M[0,0]*x[0] + M[0,1]*x[1] + M[0,2]*x[2] + M[0,3]*x[3]
    out[1] = M[1,0]*x[0] + M[1,1]*x[1] + M[1,2]*x[2] + M[1,3]*x[3]
    out[2] = M[2,0]*x[0] + M[2,1]*x[1] + M[2,2]*x[2] + M[2,3]*x[3]
    out[3] = M[3,0]*x[0] + M[3,1]*x[1] + M[3,2]*x[2] + M[3,3]*x[3]


@nb.njit(inline='always')
def _dot4_inplace(M: np.ndarray, x: np.ndarray) -> None:
    """ Dot product of two 4-length arrays to avoid temp array
    allocation. Rewrites to `x`.
    """

    a = M[0,0]*x[0] + M[0,1]*x[1] + M[0,2]*x[2] + M[0,3]*x[3]
    b = M[1,0]*x[0] + M[1,1]*x[1] + M[1,2]*x[2] + M[1,3]*x[3]
    c = M[2,0]*x[0] + M[2,1]*x[1] + M[2,2]*x[2] + M[2,3]*x[3]
    d = M[3,0]*x[0] + M[3,1]*x[1] + M[3,2]*x[2] + M[3,3]*x[3]
    x[0] = a; x[1] = b; x[2] = c; x[3] = d


@nb.njit(inline='always')
def _axpy4(scale: float, src: np.ndarray, dst: np.ndarray, col: int) -> None:
    """ For 4-length array, dst[:, col] += scale * src """

    dst[0, col] += scale * src[0]
    dst[1, col] += scale * src[1]
    dst[2, col] += scale * src[2]
    dst[3, col] += scale * src[3]


@nb.njit(inline='always')
def _dot4_axpy4(scale: float, M: np.ndarray, x: np.ndarray, dst: np.ndarray, col: int) -> None:
    """ dot4 -> axpy4 without temp allocation """

    a = M[0,0]*x[0] + M[0,1]*x[1] + M[0,2]*x[2] + M[0,3]*x[3]
    b = M[1,0]*x[0] + M[1,1]*x[1] + M[1,2]*x[2] + M[1,3]*x[3]
    c = M[2,0]*x[0] + M[2,1]*x[1] + M[2,2]*x[2] + M[2,3]*x[3]
    d = M[3,0]*x[0] + M[3,1]*x[1] + M[3,2]*x[2] + M[3,3]*x[3]

    dst[0, col] += scale * a
    dst[1, col] += scale * b
    dst[2, col] += scale * c
    dst[3, col] += scale * d


@nb.njit(inline='always')
def _sum4(x: np.ndarray) -> float:
    """ Sum a 4-length array """
    return x[0] + x[1] + x[2] + x[3]
