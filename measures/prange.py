"""Direct measures t from paper 'A complete quantum circuit to solve the
information set decoding problem' and 'TQuantum' paper.

In this improvement, we reduce the number of columns from n to an

"""
import copy
from typing import Optional, Tuple

from sage.all import assume, var
from sage.functions.all import binomial, ceil, factorial, log, max_symbolic
from sage.symbolic.all import Expression
from utils.datamanipulation import replace_exp_valuedict

from measures.common import (
    MAX_DEPTHS,
    CodeExtendedSage,
    StandardGate,
    ValueCliffordTDict,
    ValueDicts,
    ValueNormalDictInt,
    ValueNormalDictLog2,
    ValueParallelized,
    get_qaes_nist,
)
from measures.conversion import TGateConversion, convert_gate
from measures.decomposition import DecompositionType, decompose_mcx


class Prange(CodeExtendedSage):
    # optimization, alpha, delta
    def _compute_params(
        self,
        optimization: int,
        alphaa: Optional[float] = None,
        deltad: Optional[float] = None,
    ):
        alpha = var("alpha", domain="real")
        self.var_by_names["alpha"] = alpha
        assume(alpha >= 0)

        if optimization == 1:
            if not alphaa or alphaa <= 0:
                if not deltad or deltad <= 0:
                    raise Exception("Either alpha or delta must be positive")
                # Esser, Bellin
                # alpha = (1 - namespace.delta) * R
                # Our computation
                alphaa = 1 - deltad
            tmp = ceil(alpha * self.var_by_names["n_o"])
            n_new = self.var_by_names["n_o"] - tmp
            k_new = self.var_by_names["k_o"] - tmp
            r_new = self.var_by_names["n_o"] - self.var_by_names["k_o"]
            t_new = self.var_by_names["t_o"]
            classic_iters = binomial(
                self.var_by_names["n_o"], self.var_by_names["t_o"]
            ) / binomial(n_new, self.var_by_names["t_o"])
        elif optimization == 2:
            if not alphaa or alphaa <= 0:
                raise AttributeError("Alpha must be not None")
            # it's a percentage of the total weight, we fix alpha*t bit
            # positions to 1
            tmp = ceil(alpha * self.var_by_names["t_o"])
            n_new = self.var_by_names["n_o"] - tmp
            k_new = self.var_by_names["k_o"] - tmp
            r_new = self.var_by_names["n_o"] - self.var_by_names["k_o"]
            t_new = self.var_by_names["t_o"] - tmp
            # classic_iters = binomial(n, tmp)
            classic_iters = binomial(
                self.var_by_names["n_o"], self.var_by_names["t_o"]
            ) / binomial(n_new, t_new)
        else:
            # no optimization
            n_new = self.var_by_names["n_o"]
            k_new = self.var_by_names["k_o"]
            r_new = self.var_by_names["n_o"] - self.var_by_names["k_o"]
            t_new = self.var_by_names["t_o"]
            classic_iters = 1
            alphaa = 0.0

        # p stands for prime
        self.var_by_names["n_p"] = n_new
        self.var_by_names["k_p"] = k_new
        self.var_by_names["r_p"] = r_new
        self.var_by_names["t_p"] = t_new
        self.var_by_names["classic_iters"] = classic_iters
        assume(alpha > k_new / n_new)
        if self.code_ext.init_real_values:
            self.var_subs[alpha] = self.REAL_FIELD(alphaa)

    def _check_params(
        self,
        optimization: int,
        alphaa: Optional[float] = None,
        deltad: Optional[float] = None,
    ) -> bool:
        # For now it's the only check we have
        if alphaa:
            return alphaa < self.code_ext.kk / self.code_ext.nn
        return True

    def _get_mcgate_params(self) -> Tuple[Expression, Expression, Expression]:
        n = self.var_by_names["n_p"]
        r = self.var_by_names["r_p"]
        #########################
        # mcz_orac = 1 * grov_iter
        # mcz_diffc = 1 * grov_iter

        mcz_diff_ctrls = n - 1
        # if namespace.name.lower() in self.CYCLIC_CODE_NAMES:
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

    def _single_iter(self):
        n = self.var_by_names["n_p"]
        k = self.var_by_names["k_p"]
        r = self.var_by_names["r_p"]
        t = self.var_by_names["t_p"]
        nl2 = log(n, 2)
        rl2 = log(r, 2)

        M = self._get_m()
        grov_iter = self.get_grover_iterations(
            binomial(n, r), (0.2888 * binomial(n - t, k)) * M
        )

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
        cnc += (n - 1) * nl2 * (nl2 - 1)
        # cnc += 9 / 2 * r - 5 * rl2 - 11
        # HWC
        cnc_hwcom = 17 / 2 * r - 5 * rl2 - 11
        cnc += M * cnc_hwcom
        cnc *= 2
        cnc += 10 * n * r - 10 * r**2 - 4 * n
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
        normal = ValueNormalDictInt(
            grover_iters=grov_iter,
            qubits=qubits,
            gates_iteration={
                "X": xc,
                "CNOT": cnc,
                "CCNOT": ccnc,
                "RY": ryc,
            },
            depth_iteration=depth,
        )
        # normal.set_all_computable()
        return normal

    def _t_measures(
        self, normal_expanded: ValueNormalDictInt, normal_expanded2: ValueNormalDictLog2
    ):
        n = self.var_by_names["n_p"]
        k = self.var_by_names["k_p"]
        r = self.var_by_names["r_p"]
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
        # Pack columns
        tdepth += 1 / 2 * nl2 * (nl2 + 1)
        # GJE, roughly half times depth, since AND gates have depth 1 and AND^\dagger depth 0
        tdepth += 0.5 * depth_gauss
        # HW check, same as before
        tdepth += 0.5 * depth_hwcom
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

    def go(self, **kwords) -> ValueDicts:
        # namespace = kwords['namespace']
        optimization = kwords["optimization"]
        decomposition = kwords["expansion"]
        alphaa = kwords["alphaa"]
        deltad = kwords["deltad"]
        max_depthh2 = kwords.get("max_depthh2", None)
        qaess2 = kwords.get("qaess2", None)
        # if we have to substitute values or not
        subs_values = kwords["subs"]

        M = self._get_m()
        max_depth2 = self.var_by_names["max_depth2"]

        ret_dic = ValueDicts(code_ext=self.code_ext)
        if self.code_ext:
            ret_dic.code_ext = self.code_ext

        self._compute_params(optimization, alphaa, deltad)
        n = self.var_by_names["n_p"]
        k = self.var_by_names["k_p"]
        t = self.var_by_names["t_p"]
        r = self.var_by_names["r_p"]

        ####
        # Return dictionary composed only of parameterized expressions
        ####
        ret_dic.in_params = self.get_dict_input_params(subs=False)
        ret_dic.normal = self._single_iter()
        ret_dic.normal.classic_iters = self.var_by_names["classic_iters"]
        # Avoiding, too expensive
        # ret_dic.normal.set_all_computable()
        ret_dic.normal.set_gates_sum_iter()
        ret_dic.normal.set_gates_sum_total()

        ret_dic.normal2 = ret_dic.normal.to_log2()
        grov_iter2 = self.get_grover_iterations_log2(
            (factorial(n), factorial(r - t)),
            (factorial(r), factorial(n - t), M, 0.2888),
        )
        ret_dic.normal2.grover_iters = grov_iter2
        ret_dic.normal2.set_all_computable()
        # print(ret_dic.normal2.grover_iters)
        # input("n2")

        # Expanded
        ret_dic.normal_expanded = copy.deepcopy(ret_dic.normal)
        decompose_mcx(
            DecompositionType(decomposition),
            ret_dic.normal_expanded,
            *self._get_mcgate_params(),
            M
        )
        # ret_dic.normal_expanded.set_all_computable()
        ret_dic.normal_expanded2 = ret_dic.normal_expanded.to_log2()
        # WIP This is necessary to simplify some factorials
        ret_dic.normal_expanded2.grover_iters = grov_iter2
        ret_dic.normal_expanded2.set_all_computable()

        # print(ret_dic.normal_expanded2.grover_iters)
        # input("n2_exp")

        ret_dic.tmeas, ret_dic.tmeas2 = self._t_measures(
            ret_dic.normal_expanded, ret_dic.normal_expanded2
        )

        #### Parallelized version
        prange_depth2 = ret_dic.normal_expanded2.depth_total
        aes_gates = self.var_by_names["qaes2"] - self.var_by_names["max_depth2"]

        n_instances2 = max_symbolic(
            0,
            2 * ret_dic.normal_expanded2.depth_total
            # 2 * ret_dic.normal_expanded2.depth_iteration  # OK
            # - (-2 * ret_dic.normal_expanded2.grover_iters) # probs
            - max_depth2 * 2,  # OK
        )
        overall_gates = ret_dic.normal_expanded2.gates_sum_total + 0.5 * n_instances2

        # TODO convert to if else expression
        # overall_depth = (
        #     max_depth + n_instances
        #     if n_instances > 0
        #     else ret_dic.normal_expanded2.depth_total
        # )
        overall_depth = max_depth2 + n_instances2

        parallel2 = ValueParallelized(
            max_depth=max_depth2,
            qaes=self.var_by_names["qaes2"],
            aes_gates=aes_gates,
            # D^2 / MAXD
            prange_depth=2 * prange_depth2 - max_depth2,
            # D^2/QAES
            security_margin=2 * prange_depth2 - self.var_by_names["qaes2"],
            #
            n_instances=n_instances2,
            depth=overall_depth,
            gates=overall_gates,
            gates_margin=overall_gates - aes_gates,
        )
        ret_dic.parallel_expanded2_maxd = parallel2

        # print(ret_dic.tmeas2.grover_iters)
        # input("T2")
        ####

        if not subs_values:
            return ret_dic

        #####
        # Return dictionary composed of values
        #####
        is_feasible = self._check_params(optimization, alphaa, deltad)
        if not is_feasible:
            # print(f"alpha {alpha} is greater than or equal to R {kk/nn}")
            return ret_dic
        # I add also equality check since I don't want to analyse limit cases
        # also here for n, k, ....
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
        if max_depthh2:
            self.var_subs[self.var_by_names["max_depth2"]] = max_depthh2
            if qaess2:
                self.var_subs[self.var_by_names["qaes2"]] = qaess2
            else:
                self.var_subs[self.var_by_names["qaes2"]] = get_qaes_nist(
                    self.code_ext.level
                )
            parallel_expanded2_maxd_given = replace_exp_valuedict(
                ret_dic.parallel_expanded2_maxd, self.var_subs, numerical=True
            )
            ret_dic.parallel_expanded2_maxds[
                self.code_ext.max_depth_log2
            ] = parallel_expanded2_maxd_given
        else:
            # We check for all the depth given by NIST
            for max_depth_nist in MAX_DEPTHS:
                # self.var_subs[alpha] = self.REAL_FIELD(alphaa)
                running = copy.deepcopy(ret_dic.parallel_expanded2_maxd)
                self.var_subs[self.var_by_names["max_depth2"]] = self.REAL_FIELD(
                    max_depth_nist
                )
                if qaess2:
                    self.var_subs[self.var_by_names["qaes2"]] = qaess2
                else:
                    self.var_subs[self.var_by_names["qaes2"]] = get_qaes_nist(
                        self.code_ext.level
                    )

                parallel_expanded2_maxd_given = replace_exp_valuedict(
                    running, self.var_subs, numerical=True
                )
                ret_dic.parallel_expanded2_maxds[
                    max_depth_nist
                ] = parallel_expanded2_maxd_given
        ret_dic.parallel_expanded2_maxd = None
        return ret_dic
