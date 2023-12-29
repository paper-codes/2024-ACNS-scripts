"""Direct measures from phdthesis. Adaptation strategy 2.
"""
import copy
from typing import Tuple

from sage.all import sqrt, var
from sage.functions.all import binomial, factorial, log
from sage.symbolic.all import Expression
from utils.datamanipulation import replace_exp_valuedict, replace_expression

from measures.conversion import TGateConversion, convert_gate
from measures.common import (
    CodeExtendedSage,
    StandardGate,
    ValueCliffordTDict,
    ValueDicts,
    ValueNormalDictInt,
    ValueNormalDictLog2,
    ValueParallelized,
    get_qaes_nist,
)
from measures.decomposition import DecompositionType, decompose_mcx


class LeeBrickellQuantum(CodeExtendedSage):
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

        grov_iters = self.get_grover_iterations(
            binomial(n, t), 0.2888 * binomial(k, p) * binomial(r, t - p) * M
        )

        ret_dic.normal = ValueNormalDictInt(grover_iters=grov_iters, qubits=0)
        self._single_iter(ret_dic.normal)
        # too much, better to set the computable only in log2 form
        # ret_dic.normal.set_all_computable()
        ret_dic.normal.set_gates_sum_iter()

        # binomial(n_o, t_o) / ( binomial(r_o, t_o - p_o) * binomial(k_o, p_o))
        classic_iters = 1

        ret_dic.normal.classic_iters = classic_iters
        ret_dic.normal.classic_space = r * n + r * M

        # too much, better to set the computable only in log2 form
        # ret_dic.normal.set_all_computable()

        ####
        # Return dictionary composed only of parameterized expressions
        ####
        ret_dic.normal2 = ret_dic.normal.to_log2(
            skip_grover_iters=True, skip_classic_iters=True
        )

        classic_iters_l2 = 0
        grov_iter2 = self.get_grover_iterations_log2(
            (binomial(n, t),), (0.2888, binomial(k, p), binomial(r, t - p), M)
        )

        ret_dic.normal2.classic_iters = classic_iters_l2
        ret_dic.normal2.grover_iters = grov_iter2
        ret_dic.normal2.set_all_computable()

        # Expanded
        ret_dic.normal_expanded = copy.deepcopy(ret_dic.normal)
        # TODO double check
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

        ret_dic.tmeas, ret_dic.tmeas2 = self._t_measures(
            ret_dic.normal_expanded, ret_dic.normal_expanded2
        )

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
        ret_dic.tmeas = replace_exp_valuedict(
            ret_dic.tmeas, self.var_subs, numerical=True
        )
        ret_dic.tmeas2 = replace_exp_valuedict(
            ret_dic.tmeas2, self.var_subs, numerical=True
        )
        return ret_dic

    def _get_mcgate_params(self) -> Tuple[Expression, Expression, Expression]:
        n = self.var_by_names["n_o"]
        r = self.var_by_names["r_o"]
        #########################
        # mcz_orac = 1 * grov_iter
        # mcz_diffc = 1 * grov_iter

        mcz_diff_ctrls = n - 1
        if self.is_cyclic():
            # Oracle CZ
            mcx_ora_ctrls = r  # for identity

            # for each of the M controls, we have as control qubits the one
            # storing the result of HWC, +1 containing the result of identity
            mcz_ora_ctrls = log(r, 2) + 1
            # Diffusion CZ
        else:
            mcz_ora_ctrls = log(r, 2) + r
            mcx_ora_ctrls = 0

        return mcz_diff_ctrls, mcx_ora_ctrls, mcz_ora_ctrls

    def _get_m(self):
        return self.input_params["r_o"] if self.is_cyclic() else 1


    def _single_iter(self, normal: ValueNormalDictInt):
        n = self.input_params["n_o"]
        k = self.input_params["k_o"]
        r = self.input_params["r_o"]
        t = self.input_params["t_o"]
        p = self.input_params["p_o"]
        nl2 = log(n, 2)
        rl2 = log(r, 2)
        M = self._get_m()

        xc, cnc, ccnc, cswc, ryc = 0, 0, 0, 0, 0

        #########################
        # MEASURES FOR 1 round. We add the number of rounds later. We also skip
        # the initial dicke state, the difference is not so much
        #########################
        # H
        xc += n * r / 2
        # syndrome(s)
        # <5.2
        xc += (r * M) / 2  # r instances of syndromes
        # />
        # pack columns
        xc += 2 * (n - 1) * nl2 * (nl2 - 1)
        # GJE
        xc += 2 * (r - 2)  # we skip also last check pivot
        # HWC
        xc_hwcom = 4 * r - 2 * rl2 - 4
        xc_hwchk = rl2 + 1 - log(t, 2)
        xc += (xc_hwcom + xc_hwchk) * M
        # because of uncompute
        xc *= 2
        # Diffusion
        xc += n + 2 * r
        # Dicke
        xc += r
        #########################
        # GJE checkers (pivot and element under pivot)
        # Impr.: C vector not required
        # cnc_gauss = 3 / 2 * r * (r - 1)
        cnc_gauss = 1 / 2 * r * (r - 1)
        # Impr.: C vector not required
        # 1 CNOT per column addition
        # cnc_gauss += r * (r - 1)
        cnc += cnc_gauss
        # cnc += 9 / 2 * r - 5 * rl2 - 11
        # Pack columns
        cnc += (n - 1) * nl2 * (nl2 - 1)
        cnc *= 2
        # diffusione
        cnc += 10 * n * r - 10 * r**2 - 4 * n
        # L-B combination generation
        # le p somme iniziali
        cnc += p
        # Le successive ncr(k;p) - 1 combinazioni in cui tolgo una colonna e
        # aggiungo una colonna
        cnc += 2 * (binomial(k, p) - 1)
        # Per ogni combinazione, un circuito di somma
        # HWC
        cnc_hwcom = binomial(k, p) * (17 / 2 * r - 5 * rl2 - 11)
        cnc += M * cnc_hwcom
        #########################
        # Pack columns
        # ccnc += (n - 1) * nl2 * (nl2 - 1)
        # GJE
        # <5.2
        # fake swap ccnots
        ccnc_gauss = 1 / 6 * r * (r - 1) * (2 * r + 3 * M + 2)
        # column add ccnots
        ccnc_gauss += 1 / 2 * r * (r - 1) * (r + 2 * M - 1)
        # />
        ccnc += ccnc_gauss
        # HWC
        # ccnc += 3 * r - 2 * rl2 - 3
        ccnc_hwcom = 3 * r - 2 * rl2 - 3
        ccnc += M * ccnc_hwcom
        ccnc *= 2
        #########################
        cswc += (n - 1) * nl2 * (nl2 - 1) * (r + 1)
        cswc *= 2
        #########################
        ryc += 4 * n * r - 4 * r**2 - 2 * n + 1
        ryc *= 2
        ryc += 2 * (4 * n * r - 4 * r**2 - 2 * n + 1)
        #########################

        # depth = (3 * (9 * n * r - 4 * n - 9 * r**2 + 1)) / (r - 2)
        # Dicke
        depth = n
        # Data
        depth += 1
        # Pack columns
        depth += nl2 * (nl2 + 1) + r
        # Gauss
        depth_gauss = 3 / 2 * r * (r - 1) + r
        if self.is_cyclic():
            depth_gauss += r
        depth += depth_gauss
        # HWC
        # Depth is equal for HW com/chk since all the adders can run in
        # parallel
        # Cuccaro
        # depth_hwcom = rl2**2 + 7 * rl2 - 4
        # TakahashiKK10
        depth_hwcom = 2.5 * rl2**2 - 0.5 * rl2
        depth += depth_hwcom
        depth += 1
        depth *= 2
        depth += (3 * (9 * n * r - 4 * n - 9 * r**2 + 1)) / (r - 2)
        # L-B generation + hwcom. No UNDO.
        depth += binomial(k, p) * (depth_hwcom + 1)
        #########################
        # dicke
        qubits = n
        # H
        qubits += r * n
        # syndrome(s)
        qubits += r * M
        # permutation
        qubits += (n - 1) * nl2 * (nl2 - 1)
        # GJE ancillae
        # qubits += 3 / 2 * r * (r - 1)
        qubits += r * (r + 1) / 2
        # HWC
        # Cuccaro
        # qubits += M * (r + r / 4 - 1)
        # TakahashiKK10
        qubits += M * r
        # for the oracle flip, adding the identity
        if self.is_cyclic():
            qubits += 1
        #########################

        # we do not report the multi-controlled CZs here since they have
        # different number of controls and therefore it has no meaning
        normal.qubits = qubits
        normal.gates_iteration={
            "X": xc,
            "CNOT": cnc,
            "CCNOT": ccnc,
            "RY": ryc,
            "CSWAP": cswc,
            }
        normal.depth_iteration = depth

    def _t_measures(
        self, normal_expanded: ValueNormalDictInt, normal_expanded2: ValueNormalDictLog2
    ):
        n = self.var_by_names["n_o"]
        k = self.var_by_names["k_o"]
        r = self.var_by_names["r_o"]
        p = self.var_by_names["p_o"]
        nl2 = log(n, 2)
        rl2 = log(r, 2)
        # ##############
        # # T measures
        # ##############
        t_dict: ValueCliffordTDict
        # relevant_gates = ('X', 'CCNOT', 'CCSWAP', 'CZ', 'Z')
        gates_total_iter = {}
        # add_qubits is useless. It gives the number of additional qubits
        # required in each conversion. However, such qubits are reset at the
        # end of the conversion, so they can be reused by the other conversions.
        add_qubits = 0
        for gate_name, gate_value in normal_expanded.gates_iteration.items():
            gates: TGateConversion = convert_gate(StandardGate(gate_name, gate_value))
            for key, val in gates._asdict().items():
                if key == "add_qubits":
                    add_qubits += val
                else:
                    # v1 = v if v is not None else 0
                    gates_total_iter[key] = gates_total_iter.get(key, 0) + val
        gates_total = {
            key: val * normal_expanded.grover_iters
            for key, val in gates_total_iter.items()
        }
        gates_total2 = {
            key: log(val, 2) + normal_expanded2.grover_iters
            for key, val in gates_total_iter.items()
        }

        # For the distinct depths, check the notes of 2021_IQW_ISD: expres.tex
        depth_hwcom = rl2**2 + 7 * rl2 - 4
        depth_gauss = 3 / 2 * r * (r - 1) + r
        if self.is_cyclic():
            depth_gauss += r

        # Dicke
        tdepth = n
        # Pack columns.
        #
        # .5 is because of roughly half times depth, since AND gates have depth
        # 1 and AND^\dagger depth 0 
        tdepth += .5 * nl2 * (nl2 + 1)
        # GJE, same as before
        tdepth += 0.5 * depth_gauss
        # HW check, same as before
        tdepth += 0.5 * binomial(k, p) * depth_hwcom

        # MCZ ~\cite{DBLP:journals/corr/abs-1210-0974} showed a decomposition
        # requiring $8m$ \qgate{T} gates and a \qgate{T}-depth of $2\lbp{m} +
        # 1$, using the same amount of additional qubits of our previous
        # decomposition. tdepth -= (cz1ctrls + cz2ctrls)
        mcz_diff_cs, mcx_ora_cs, mcz_ora_cs = self._get_mcgate_params()

        tdepth_mcz = log(mcz_diff_cs, 2) + 1
        if self.is_cyclic():
            # + 1 because mcnot set an additional qubit acting as control

            # MCX seems to have  log() decomposition
            tdepth_mcz += log(mcx_ora_cs, 2)
            # rough estimate, based on the fact that they can all be run in
            # parallel by interleaving controls
            tdepth_mcz += 2 * log(mcz_ora_cs, 2)
        else:
            tdepth_mcz += 2 * log(mcz_ora_cs, 2)
        tdepth += tdepth_mcz

        tdepth_norm = tdepth * normal_expanded.grover_iters
        tdepth_log2 = log(tdepth, 2) + normal_expanded2.grover_iters
        del tdepth

        # add_qubits
        add_qubits = r**2
        # last dicke
        # tdepth += n
        # tdepth_c = tdepth * classic_iters
        t_dict = ValueCliffordTDict(
            qubits=normal_expanded.qubits + add_qubits,
            # qubits=normal_expanded.qubits,
            grover_iters=normal_expanded.grover_iters,
            gates_total=gates_total,
            t_count=gates_total["T"],
            t_depth=tdepth_norm,
            tq=normal_expanded.qubits * tdepth_norm,
        )
        qubits2 = log(normal_expanded.qubits + add_qubits, 2)
        t_dict2 = ValueCliffordTDict(
            qubits=qubits2,
            # qubits=normal_expanded2.qubits,
            grover_iters=normal_expanded2.grover_iters,
            gates_total=gates_total2,
            t_count=gates_total2["T"],
            t_depth=tdepth_log2,
            tq=qubits2 + tdepth_log2,
        )
        return t_dict, t_dict2





