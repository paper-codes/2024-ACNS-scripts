from enum import Enum
from typing import Dict, List, Tuple

from sage.symbolic.all import Expression

from measures.common import ValueNormalDictInt


class DecompositionType(Enum):
    OURS = "ours"
    BARENCO = "barenco"
    SAEEDI = "saeedi"


def decompose_mcx(
    decomposition: DecompositionType,
    dic_to_expand: ValueNormalDictInt,
    mcz_diff_ctrls: Expression,
    mcx_ora_ctrls: Expression,
    mcz_ora_ctrls: Expression,
    M: Expression,
):
    """Returns additional depth and qubits"""
    # (gates_add: {gate_name: value}, qubits_add: exp, depth_add: exp)
    add_list: List[Tuple[Dict[str, Expression], Expression, Expression]]
    if decomposition == DecompositionType.OURS:
        diff_dec = _decompose_mcx1(mcz_diff_ctrls, True)
        add_list = [diff_dec]
        if M > 1:
            ora_add = _decompose_mcx1(mcx_ora_ctrls, False)
            add_list.append(ora_add)
            ora_add = _decompose_mcx1(mcz_ora_ctrls, True)
            for k, v in ora_add[0].items():
                ora_add[0][k] = M * v
            # Depth is not multiplied by M since all the MCZ can be interleaved
            add_list.append((ora_add[0], ora_add[1] * M, ora_add[2]))
        else:
            ora_add = _decompose_mcx1(mcz_ora_ctrls, True)
    elif decomposition == DecompositionType.BARENCO:
        raise NotImplementedError("Barenco decomposition not implemented yet")
    elif decomposition == DecompositionType.SAEEDI:
        # Diffusion expansion
        diff_add = _decompose_mcx3(mcz_diff_ctrls, True)
        add_list = [diff_add]
        # Oracle expansion
        if M > 1:
            ora_add = _decompose_mcx3(mcx_ora_ctrls, False)
            add_list.append(ora_add)
            ora_add = _decompose_mcx3(mcz_ora_ctrls, True)
            for k, v in ora_add[0].items():
                ora_add[0][k] = M * v
            # Depth is not multiplied by M since all the MCZ can be interleaved
            add_list.append((ora_add[0], ora_add[1] * M, ora_add[2]))

        else:
            ora_add = _decompose_mcx3(mcz_ora_ctrls + mcz_ora_ctrls, True)
            add_list.append(ora_add)

    for dic_gates, qubits, depth in add_list:
        for k, v in dic_gates.items():
            dic_to_expand.gates_iteration[k] = (
                dic_to_expand.gates_iteration.get(k, 0) + v
            )
        dic_to_expand.qubits += qubits
        dic_to_expand.depth_iteration += depth


def _decompose_mcx1(
    n_ctrls: Expression, final_cz: bool
) -> Tuple[Dict[str, Expression], Expression, Expression]:
    # For mCNOT, instead, 2(m-2) CCNOTs, (m-2) additional qubits and depth
    # log(m) + 1.
    # dic_to_expand['CCNOT'] = dic_to_expand['CCNOT'] + 2 * (n_ctrls - 1)
    # dic_to_expand['qubits'] = dic_to_expand['qubits'] + (n_ctrls - 1)
    # dic_to_expand['depth'] = dic_to_expand['depth'] + log(n_ctrls, 2)
    ccnot_add = 2 * (n_ctrls - 1)
    qubits_add = n_ctrls - 1
    # depth_add = log(n_ctrls, 2)
    depth_add = 2 * n_ctrls
    gates_add = {"CCNOT": ccnot_add}
    if final_cz:
        gates_add["CZ"] = int(1)
    return gates_add, qubits_add, depth_add


def _decompose_mcx2(
    n_ctrls: Expression,
) -> Tuple[Dict[str, Expression], Expression, Expression]:
    # # Barenco et al. 1995
    # ccnc_diff = (mcz_diff_ctrls - 1)**2
    # # cyclic
    # if M > 1:
    #     ccnc_ora = (mcx_ora_ctrls - 2)**2
    #     ccnc_ora += (M * (mcz_ora_ctrls - 1))**2
    # else:
    #     ccnc_ora = 2 * (mcz_ora_ctrls - 1)
    # qubits_exp = dic_to_expand['qubits']
    # dic_to_expand['CCNOT'] = dic_to_expand['CCNOT'] + ccnc_diff + ccnc_ora
    # # 1 MCZ for the diffusion stage and M for the oracle at each iteration
    # dic_to_expand['CZ'] = dic_to_expand.get('CZ', 0) + 1
    # dic_to_expand['CZ'] = dic_to_expand.get('CZ', 0) + M
    # depth_exp = dic_to_expand['depth'] + log(mcz_diff_ctrls, 2)
    # # cyclic
    # if M > 1:
    #     # Depth is not multiplied by M since all the MCZ can be interleaved
    #     # (they have only 1 qubit in common ouof t of (M-1))
    #     depth_exp += log(mcx_ora_ctrls, 2) + log(mcz_ora_ctrls, 2)
    # else:
    #     depth_exp += log(mcz_ora_ctrls, 2)
    pass


def _decompose_mcx3(
    n_ctrls: Expression, final_z: bool
) -> Tuple[Dict[str, Expression], Expression, Expression]:
    # x is the number of C-Rx required
    gates_add = {"CRX": 2 * (n_ctrls) ** 2 - 6 * (n_ctrls) + 5}
    if final_z:
        gates_add["Z"] = int(1)
    qubits_add = 1
    depth_add = 8 * n_ctrls - 20

    return gates_add, qubits_add, depth_add

    # # Each mCZ is applied once at each grover iterations. An
    # # m-controlled CZ is expanded in 2(m-1) CCNOTs, 1 CZ, (m-1)
    # # additional qubits and depth log(m) + 2. But we can consider
    # # log(m) since the last 2 layers consist of only 1 gate each that
    # # can be interleaved w/ other parts of the circuit.

    # # Note that the M>1 mCZ have have log(r) + 1 ctrls, where log(r)
    # # comes from the HWCC circuit and 1 from the output of the identity
    # # check. If M = 1, instead, we can use a single mCZ having log(r) +
    # # r control qubits.

    # # For mCNOT, instead, 2(m-2) CCNOTs, (m-2) additional qubits and depth
    # # log(m) + 1.
    # ccnc_diff = 2 * (mcz_diff_ctrls - 1)
    # qubits_exp = dic_to_expand['qubits'] + (mcz_ora_ctrls +
    #                                         mcz_diff_ctrls - 2)
    # if M > 1:
    #     ccnc_ora = 2 * (mcx_ora_ctrls - 2)  # iden check
    #     ccnc_ora += 2 * M * (mcz_ora_ctrls - 1)  # M HWCCs
    #     qubits_exp += mcx_ora_ctrls - 2
    # else:
    #     ccnc_ora = 2 * (mcz_ora_ctrls - 1)
    # dic_to_expand['CCNOT'] = dic_to_expand['CCNOT'] + ccnc_diff + ccnc_ora
    # # 1 MCZ for the diffusion stage and M for the oracle at each iteration
    # dic_to_expand['CZ'] = dic_to_expand.get('CZ', 0) + 1
    # dic_to_expand['CZ'] = dic_to_expand.get('CZ', 0) + M
    # #
    # depth_exp = dic_to_expand['depth'] + log(mcz_diff_ctrls, 2)
    # if M > 1:
    #     # Depth is not multiplied by M since all the MCZ can be interleaved
    #     # (they have only 1 qubit in common ouof t of (M-1))
    #     depth_exp += log(mcx_ora_ctrls, 2) + log(mcz_ora_ctrls, 2)
    # else:
    #     depth_exp += log(mcz_ora_ctrls, 2)
