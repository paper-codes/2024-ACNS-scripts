"""Lee-Brickell classic. Measures taken from ??
Can be used also to derive Prange's ISD complexity by setting p=0.
"""

from sage.all import sqrt, var
from sage.functions.all import binomial, log
from utils.datamanipulation import replace_exp_valuedict

from measures.common import CodeExtendedSage, ValueDicts, ValueNormalDictInt


class LeeBrickellClassic(CodeExtendedSage):
    def _init_additional(self, **kwords):
        p_o = var("p_o", domain="integer")
        self.var_by_names["p_o"] = p_o
        self.input_params["p_o"] = p_o
        self.var_subs[p_o] = kwords["p"]

    def go(self, **kwords) -> ValueDicts:
        n = self.input_params["n_o"]
        k = self.input_params["k_o"]
        r = self.input_params["r_o"]
        t = self.input_params["t_o"]
        p = self.input_params["p_o"]
        M = self._get_m()

        ret_dic = ValueDicts(code_ext=self.code_ext)
        if self.code_ext:
            ret_dic.code_ext = self.code_ext
        ret_dic.in_params = self.get_dict_input_params(subs=False)

        ret_dic.normal = ValueNormalDictInt(
            grover_iters=None,
            qubits=None,
            gates_iteration=None,
            depth_iteration=None,
        )

        classic_iters = binomial(n, t) / binomial(r, t - p) / binomial(k, p)
        # DOOM
        classic_iters /= sqrt(M)
        ret_dic.normal.classic_iters = classic_iters
        ret_dic.normal.classic_space = r * n + r * M
        ret_dic.normal.classic_gates_iter = self._classic_gates_iter()
        ret_dic.normal.classic_gates_total = (
            ret_dic.normal.classic_gates_iter * ret_dic.normal.classic_iters
        )
        ret_dic.normal.classic_space_time = (
            ret_dic.normal.classic_gates_total * ret_dic.normal.classic_space
        )

        ####
        # Return dictionary composed only of parameterized expressions
        ####
        ret_dic.normal2 = ret_dic.normal.to_log2(
            skip_grover_iters=True, skip_classic_iters=True
        )

        classic_iters_l2 = 0
        classic_iters_l2 += log(binomial(n, t), 2)
        classic_iters_l2 -= log(binomial(k, p), 2)
        classic_iters_l2 -= log(binomial(r, t - p), 2)
        classic_iters_l2 -= log(0.288, 2)

        ret_dic.normal2.classic_iters = classic_iters_l2
        ret_dic.normal2.set_all_computable()

        if not kwords["subs"]:
            return ret_dic

        ret_dic.in_params = self.get_dict_input_params(subs=True)
        ret_dic.normal = replace_exp_valuedict(
            ret_dic.normal, self.var_subs, numerical=True
        )
        ret_dic.normal2 = replace_exp_valuedict(
            ret_dic.normal2, self.var_subs, numerical=True
        )
        return ret_dic

    def _get_m(self):
        return self.input_params["r_o"] if self.is_cyclic() else 1

    def _get_gje1(self, n, r):
        # I computed this, but I don't remember how. Anyway, the results are
        # quite similar to Baldi et al. 2019
        return (3 / 4) * (n * r**2) + n * r / 4 - n / 2 + 3 / 4 * r**2 - r / 2

    def _get_gje2(self, n, r):
        # From Baldi et al. 2019, A Finite regime analysis ...
        return (1 / 2) * (n * r**2) + n * r / 2 - r**3 / 6 + r**2 + r / 6 - 1


    def _classic_gates_iter(self):
        """I think I computed these on my own"""
        n = self.input_params["n_o"]
        k = self.input_params["k_o"]
        r = self.input_params["r_o"]
        p = self.input_params["p_o"]
        rl2 = log(r, 2)
        M = self._get_m()

        classic_gates = 0
        # GJE
        #
        # n + M since the M syndromes are appended to H in the GJE algorithm.
        #
        # classic_gates += self._get_gje1(n + M, r)
        classic_gates += self._get_gje2(n + M, r)

        # weight check. I think is given by the fast population count circuit
        # measures as 1 layer of half adders + r layer of full adders??
        classic_gates_exhaustive = 7 * r * M + rl2 * M
        # summing ncr(k;p) columns of length r and weight p to M syndromes
        classic_gates_exhaustive += p * r * M
        classic_gates += (classic_gates_exhaustive * binomial(k, p))

        return classic_gates
