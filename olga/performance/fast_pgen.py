from typing import List, Union
from types import MethodType
import numpy as np

from olga.generation_probability import GenerationProbabilityVJ, GenerationProbabilityVDJ
from olga.performance.kernels import (
    compute_Pi_R_one_numba,
    compute_Pi_L_numba,
    compute_Pi_JinsDJ_given_D_numba,
    compute_Pi_V_insVJ_given_J_numba
)

try:
    import numba as nb
    NUMBA_INSTALLED = True
except ImportError:
    NUMBA_INSTALLED = False


class FastPgen:
    """ An wrapper around OLGA's `GenerationProbabilityVJ` or `GenerationProbabilityVDJ`
    object for calculating Pgen with a faster implementation of `compute_CDR3_pgen` thanks 
    to numba.

    Call any compute method such as `compute_aa_CDR3_pgen` or `compute_nt_CDR3_pgen`
    as usual, and the faster `compute_CDR3_pgen` will be used under the hood.
 
    Parameters
    ----------
    impl : Union[GenerationProbabilityVJ, GenerationProbabilityVDJ]
        An existing OLGA generation probability model.

    Examples
    --------
    >>> from olga import GenerationProbabilityVJ
    >>> base_model = GenerationProbabilityVJ(...)
    >>> fast_model = FastPgen(base_model)
    >>> p = pgen_model.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30', 'TRBJ1-2')
    """

    def __init__(self, impl: Union[GenerationProbabilityVJ, GenerationProbabilityVDJ]) -> None:
        self._impl = impl

        if not NUMBA_INSTALLED:
            raise RuntimeError('FastPgen requires numba. ')

        self._init_common()
        if isinstance(self._impl, GenerationProbabilityVJ):
            self._init_vj()
        elif isinstance(self._impl, GenerationProbabilityVDJ):
            self._init_vdj()
        else:
            raise ValueError('Invalid model passed to FastPgen')

        # Monkey patch to override the default compute_CDR3_pgen
        def _proxy(impl_self, *args, **kwargs):
            return self.compute_CDR3_pgen(*args, **kwargs)
        impl.compute_CDR3_pgen = MethodType(_proxy, impl)
    

    def _init_common(self) -> None:
        self._nt_lut = np.full(256, -1, dtype=np.int8)
        for i, c in enumerate(b'ACGT'):
            self._nt_lut[c] = i

        self.i2aa = list(self.codons_dict.keys())
        self.aa2i = {aa:i for i, aa in enumerate(self.i2aa)}

        self._aa_lut = np.full(256, -1, dtype=np.int16)
        for aa, i in self.aa2i.items():
            self._aa_lut[ord(aa)] = i


    def _init_vj(self) -> None:
        # Need these to be pure numpy arrays to use in numba
        self.Svj_stack = np.stack([self.Svj[a] for a in self.i2aa])
        self.Dvj_stack = np.stack([self.Dvj[a] for a in self.i2aa])
        self.Tvj_stack = np.stack([self.Tvj[a] for a in self.i2aa])
        self.lDvj_stack = np.stack([self.lDvj[a] for a in self.i2aa])
        self.lTvj_stack = np.stack([self.lTvj[a] for a in self.i2aa])


    def _init_vdj(self) -> None:
        # Need these to be pure numpy arrays to use in numba
        self.Svd_stack = np.stack([self.Svd[a] for a in self.i2aa])
        self.Dvd_stack = np.stack([self.Dvd[a] for a in self.i2aa])
        self.Tvd_stack = np.stack([self.Tvd[a] for a in self.i2aa])
        self.lDvd_stack = np.stack([self.lDvd[a] for a in self.i2aa])
        self.lTvd_stack = np.stack([self.lTvd[a] for a in self.i2aa])

        self.Sdj_stack = np.stack([self.Sdj[a] for a in self.i2aa])
        self.Ddj_stack = np.stack([self.Ddj[a] for a in self.i2aa])
        self.Tdj_stack = np.stack([self.Tdj[a] for a in self.i2aa])
        self.rDdj_stack = np.stack([self.rDdj[a] for a in self.i2aa])
        self.rTdj_stack = np.stack([self.rTdj[a] for a in self.i2aa])

        # Preparing the required variables for compute_Pi_R (the slowest step)
        self.D_pre = []
        for D_in, cutD in enumerate(self.cutD_genomic_CDR3_segs):
            dseq = self._nt_lut[np.frombuffer(cutD.encode('ascii'), dtype=np.uint8)]
            PD_2nd_all = self.PD_2nd_nt_pos_per_aa_vec[D_in]
            PD_2nd_stack = np.empty((len(self.i2aa), 4, dseq.shape[0]))
            for aidx, aa in enumerate(self.i2aa):
                PD_2nd_stack[aidx] = PD_2nd_all[aa]

            self.D_pre.append({
                'dseq': dseq,
                'Pdel_D': self.PdelDldelDr_given_D[:, :, D_in],
                'max_del_D': np.asarray(self.max_delDl_given_DdelDr[D_in]),
                'min_del_D': np.asarray(self.min_delDl_given_DdelDr[D_in]),
                'PD_nt': self.PD_nt_pos_vec[D_in],
                'PD_2nd_stack': PD_2nd_stack,
                'zeroD': self.zeroD_given_D[D_in],
            })

        # Converting the codon lookup tables into arrays so we can access
        # in the numba functions
        num_aa = len(self.i2aa)
        self.allow1 = np.zeros((num_aa,  4), dtype=np.uint8)
        self.allow2 = np.zeros((num_aa, 16), dtype=np.uint8)
        self.allow3 = np.zeros((num_aa, 64), dtype=np.uint8)
        self.allow_lsf_stack = np.zeros((num_aa,  4, 4, 4), dtype=np.uint8)

        self.nt2i = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        for aa, aidx in self.aa2i.items():
            for sfx in self.sub_codons_right[aa]:
                if len(sfx) == 1:
                    l = self.nt2i[sfx[0]]
                    self.allow1[aidx, l] = 1
                elif len(sfx) == 2:
                    s = self.nt2i[sfx[0]]
                    l = self.nt2i[sfx[1]]
                    self.allow2[aidx, s*4 + l] = 1

        for aa, aidx in self.aa2i.items():
            for c in self.codons_dict[aa]:
                l = self.nt2i[c[0]]
                s = self.nt2i[c[1]]
                f = self.nt2i[c[2]]
                self.allow3[aidx, l*16 + s*4 + f] = 1
                self.allow_lsf_stack[aidx, l, s, f] = 1


    def __getattr__(self, name):
        return getattr(self._impl, name)


    def compute_CDR3_pgen(
        self,
        CDR3_seq: str,
        V_usage_mask: List[int],
        J_usage_mask: List[int]
    ) -> float:

        if isinstance(self._impl, GenerationProbabilityVDJ):
            Pi_V, max_V_align = self._impl.compute_Pi_V(CDR3_seq, V_usage_mask)
            Pi_L = self.compute_Pi_L(CDR3_seq, Pi_V, max_V_align)

            Pi_J_given_D, max_J_align = self._impl.compute_Pi_J_given_D(CDR3_seq, J_usage_mask)
            Pi_JinsDJ_given_D = self.compute_Pi_JinsDJ_given_D(CDR3_seq, Pi_J_given_D, max_J_align)
            Pi_R = self.compute_Pi_R(CDR3_seq, Pi_JinsDJ_given_D)
   
            pgen = np.sum(Pi_L[:, :-1] * Pi_R[:, 1:])
        else:
            Pi_J, r_J_usage_mask = self._impl.compute_Pi_J(CDR3_seq, J_usage_mask)
            if len(Pi_J) == 0:
                return 0.0
            
            Pi_V_given_J, max_V_align = self._impl.compute_Pi_V_given_J(CDR3_seq, V_usage_mask, r_J_usage_mask)
            if len(Pi_V_given_J) == 0:
                return 0.0
            
            Pi_V_insVJ_given_J = self.compute_Pi_V_insVJ_given_J(CDR3_seq, Pi_V_given_J, max_V_align)    
            Pi_J = np.stack(Pi_J)
            pgen = np.dot(Pi_V_insVJ_given_J[:, :, :-1].ravel(), Pi_J[:, :, 1:].ravel())

        return pgen


    def compute_Pi_V_insVJ_given_J(
        self,
        CDR3_seq: str,
        Pi_V_given_J: List[np.ndarray],
        max_V_align: int
    ) -> np.ndarray:
        L3 = 3 * len(CDR3_seq)
        numD = len(Pi_V_given_J)

        out = np.zeros((numD, 4, L3))
        Pi_V_given_J = np.stack(Pi_V_given_J, axis=0)
        aa_idx_seq = self._aa_lut[np.frombuffer(CDR3_seq.encode('ascii'), np.uint8)]
        
        compute_Pi_V_insVJ_given_J_numba(
            out,
            aa_idx_seq,
            Pi_V_given_J,
            self.PinsVJ,
            self.first_nt_bias_insVJ,
            self.zero_nt_bias_insVJ,
            self.Svj_stack, self.Dvj_stack, self.Tvj_stack,
            self.lDvj_stack, self.lTvj_stack,
            max_V_align
        )
        return out


    def compute_Pi_JinsDJ_given_D(
        self,
        CDR3_seq: str,
        Pi_J_given_D: List[np.ndarray],
        max_J_align: int
    ) -> np.ndarray:
        L3 = 3 * len(CDR3_seq)
        numD = len(Pi_J_given_D)

        out = np.zeros((numD, 4, L3))
        Pi_J_given_D = np.stack(Pi_J_given_D, axis=0)
        aa_idx_seq = self._aa_lut[np.frombuffer(CDR3_seq.encode('ascii'), np.uint8)]

        compute_Pi_JinsDJ_given_D_numba(
            out,
            aa_idx_seq,
            Pi_J_given_D,
            self.PinsDJ,
            self.first_nt_bias_insDJ,
            self.zero_nt_bias_insDJ,
            self.Sdj_stack, self.Ddj_stack, self.Tdj_stack,
            self.rDdj_stack, self.rTdj_stack,
            max_J_align
        )
        return out
 

    def compute_Pi_R(self, CDR3_seq: str, Pi_JinsDJ_given_D: np.ndarray) -> np.ndarray:
        L3 = len(CDR3_seq) * 3
        Pi_R = np.zeros((4, L3))
        aa_idx_seq = self._aa_lut[np.frombuffer(CDR3_seq.encode('ascii'), np.uint8)]

        allow_lsf_stack = self.allow_lsf_stack
        allow1 = self.allow1
        allow2 = self.allow2
        allow3 = self.allow3

        for D_in, pre in enumerate(self.D_pre):
            compute_Pi_R_one_numba(
                Pi_R, L3, -L3,
                pre['dseq'],
                Pi_JinsDJ_given_D[D_in],
                pre['Pdel_D'],
                pre['max_del_D'],
                pre['min_del_D'],
                pre['PD_nt'], pre['PD_2nd_stack'], pre['zeroD'],
                aa_idx_seq, allow_lsf_stack, allow1, allow2, allow3
            )

        return Pi_R
    

    def compute_Pi_L(self, CDR3_seq: str, Pi_V: np.ndarray, max_V_align: int) -> np.ndarray:
        L3 = len(CDR3_seq) * 3
        Pi_L = np.zeros((4, L3))
        aa_idx_seq = self._aa_lut[np.frombuffer(CDR3_seq.encode('ascii'), np.uint8)]

        compute_Pi_L_numba(
            Pi_L, aa_idx_seq, Pi_V,
            self.PinsVD,
            self.first_nt_bias_insVD,
            self.zero_nt_bias_insVD,
            self.Svd_stack, self.Dvd_stack, self.Tvd_stack,
            self.lDvd_stack, self.lTvd_stack,
            max_V_align
        )
        return Pi_L
