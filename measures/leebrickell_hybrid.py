"""Direct measures from phdthesis. Adaptation strategy 4.
"""
import copy
from typing import Tuple

from sage.all import sqrt, var
from sage.functions.all import binomial, factorial, log
from sage.symbolic.all import Expression
from utils.datamanipulation import replace_exp_valuedict, replace_expression

from measures.common import (
    CodeExtendedSage,
    ValueDicts,
    ValueNormalDictInt,
)
from measures.decomposition import DecompositionType, decompose_mcx


class LeeBrickellHybrid(CodeExtendedSage):
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

        grov_iters = self.get_grover_iterations(binomial(k, p), M)
        # anticipate computation of grover iters, if 0 it's useless to continue
        grover_iters = replace_expression(grov_iters, self.var_subs, numerical=True)
        if grover_iters <= 1:
            quantum_off = True
            ret_dic.normal = ValueNormalDictInt(
                grover_iters=None,
                qubits=None,
                gates_iteration=None,
                depth_iteration=None,
            )
        else:
            quantum_off = False
            ret_dic.normal = ValueNormalDictInt(grover_iters=grov_iters, qubits=0)
            self._single_iter(ret_dic.normal)
            # too much, better to set the computable only in log2 form
            # ret_dic.normal.set_all_computable()
            ret_dic.normal.set_gates_sum_iter()

        # binomial(n_o, t_o) / ( binomial(r_o, t_o - p_o) * binomial(k_o, p_o))
        classic_iters = (
            (
                factorial(n)
                * factorial(t - p)
                * factorial(r - t + p)
                * factorial(p)
                * factorial(k - p)
            )
            / factorial(t)
            / factorial(n - t)
            / factorial(r)
            / factorial(k)
            # already included in the GJE cost in Baldi et al. 2019
            # / 0.288
        )
        if quantum_off:
            classic_iters /= sqrt(M)

        ret_dic.normal.classic_iters = classic_iters
        ret_dic.normal.classic_space = r * n + r * M

        ret_dic.normal.classic_gates_iter = self._classic_gates_iter(quantum_off)
        ret_dic.normal.classic_gates_total = (
            ret_dic.normal.classic_gates_iter * ret_dic.normal.classic_iters
        )
        ret_dic.normal.classic_space_time = (
            ret_dic.normal.classic_gates_total * ret_dic.normal.classic_space
        )
        # too much, better to set the computable only in log2 form
        # ret_dic.normal.set_all_computable()

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
        if quantum_off:
            classic_iters_l2 -= log(0.288, 2)
            grov_iter2 = 0
        else:
            grov_iter2 = self.get_grover_iterations_log2(
                (factorial(k),),
                (
                    factorial(p),
                    factorial(k - p),
                    self._get_m(),
                ),
            )

        ret_dic.normal2.classic_iters = classic_iters_l2
        ret_dic.normal2.grover_iters = grov_iter2
        ret_dic.normal2.set_all_computable()

        # Expanded
        ret_dic.normal_expanded = copy.deepcopy(ret_dic.normal)
        if not quantum_off:
            decompose_mcx(
                DecompositionType("ours"),
                ret_dic.normal_expanded,
                *self._get_mcgate_params(),
                self._get_m()
            )
        # too much, better to set the computable only in log2 form
        # ret_dic.normal_expanded.set_all_computable()
        ret_dic.normal_expanded2 = ret_dic.normal_expanded.to_log2(
            skip_grover_iters=True,
            skip_classic_iters=True,
        )
        ret_dic.normal_expanded2.grover_iters = grov_iter2
        ret_dic.normal_expanded2.classic_iters = classic_iters_l2
        ret_dic.normal_expanded2.set_all_computable()

        if not kwords["subs"]:
            return ret_dic

        ret_dic.in_params = self.get_dict_input_params(subs=True)
        ret_dic.normal = replace_exp_valuedict(
            ret_dic.normal, self.var_subs, numerical=True
        )
        ret_dic.normal2 = replace_exp_valuedict(
            ret_dic.normal2, self.var_subs, numerical=True
        )
        ret_dic.normal_expanded = replace_exp_valuedict(
            ret_dic.normal_expanded, self.var_subs, numerical=True
        )
        ret_dic.normal_expanded2 = replace_exp_valuedict(
            ret_dic.normal_expanded2, self.var_subs, numerical=True
        )
        return ret_dic

    def _get_mcgate_params(self) -> Tuple[Expression, Expression, Expression]:
        mcz_diff_ctrls = self.var_by_names["k_o"] - 1
        mcz_ora_ctrls = log(self.var_by_names["r_o"], 2) + 1
        return mcz_diff_ctrls, 0, mcz_ora_ctrls

    def _get_m(self):
        return self.input_params['r_o'] if self.is_cyclic() else 1

    # Ignored
    # def _input_prep(self):
    #     nl2 = log(n, 2)
    #     rl2 = log(r, 2)
    #     M = self._get_m()
    #     xc, cnc, ccnc, cswc, ryc = 0, 0, 0, 0, 0

    #     # Dicke
    #     xc += p
    #     cnc += 5*k*r - 5 * (p**2) - 2*k
    #     ryc += 4 * k * p - 4 * p**2 - 2 * k + 1
    #     depth = k

    def _classic_gates_iter(self, quantum_off: bool):
        n = self.input_params['n_o']
        k = self.input_params['k_o']
        r = self.input_params['r_o']
        p = self.input_params['p_o']
        rl2 = log(r, 2)
        M = self._get_m()

        classic_gates = 0
        # GJE
        classic_gates += (
            (3 / 4) * (n * r**2) + n * r / 4 - n / 2 + 3 / 4 * r**2 - r / 2
        )
        # compute on s
        classic_gates += r**2
        # weight check
        classic_gates += 7 * r + rl2
        if quantum_off:
            # summing ncr(k;p) columns of length r and weight p to M syndromes
            classic_gates_exhaustive = (p * r * M) * binomial(k, p)
            classic_gates += classic_gates_exhaustive

        return classic_gates

    def _single_iter(self, normal: ValueNormalDictInt):
        p = self.input_params['p_o']
        k = self.input_params['k_o']
        r = self.input_params['r_o']
        t = self.input_params['t_o']
        rl2 = log(r, 2)
        M = self._get_m()

        #########################
        xc, cnc, ccnc, ryc = 0, 0, 0, 0
        #########################
        xc += r * M / 2
        xc_hwcom = 4 * r - 2 * rl2 - 4
        xc_hwchk = rl2 + 1 - log(t - p, 2)
        xc += (xc_hwcom + xc_hwchk) * M
        xc *= 2
        # Diffusion
        xc += r + 2 * p
        #########################
        cnc += r * k * M / 2
        cnc_hwcom = 9 / 2 * r - 5 * rl2 - 11
        cnc += M * cnc_hwcom
        cnc += 10 * k * p - 10 * p**2 - 4 * k
        cnc *= 2
        #########################
        ccnc_hwcom = 3 * r - 2 * rl2 - 3
        ccnc += M * ccnc_hwcom
        ccnc *= 2
        #########################
        ryc_dicke = 4 * k * p - 4 * p**2 - 2 * k + 1
        ryc += 2 * ryc_dicke
        #########################
        depth = 1
        # Sum
        depth += k
        # HWCC
        depth_hwcom = rl2**2 + 7 * rl2 - 4
        depth += depth_hwcom
        # diff
        depth += k

        depth *= 2
        # flip
        depth += M
        #########################
        qubits = k
        # syndrome(s)
        qubits += r * M
        qubits += M * (5 * r / 4 - 1)

        #########################
        normal.qubits = qubits
        normal.gates_iteration = {
            "X": xc,
            "CNOT": cnc,
            "CCNOT": ccnc,
            "RY": ryc,
        }
        normal.depth_iteration = depth
