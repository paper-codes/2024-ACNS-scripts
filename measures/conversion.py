""" Convert to Clifford + T gate set"""

from typing import NamedTuple

from measures.common import StandardGate


class TGateConversion(NamedTuple):
    S: int = 0
    H: int = 0
    CNOT: int = 0
    T: int = 0
    add_qubits: int = 0


def convert_gate(gate: StandardGate) -> TGateConversion:
    match gate.name:
        case "X":
            res = _convert_gate_x()
        case "CCNOT":
            res = _convert_gate_cCNOT2()
        case "CSWAP":
            res = _convert_gate_cswap2()
        case "CZ":
            res = _convert_gate_cz()
        case "RY":
            res = _convert_gate_ry()
        case "H":
            res = TGateConversion(H=1)
        case "CNOT":
            res = TGateConversion(CNOT=1)
        case "T":
            res = TGateConversion(T=1)
        case "S":
            res = TGateConversion(S=1)
        case "Z":
            res = TGateConversion(S=2)
        case "CRX":
            res = _convert_gate_crx()
        case _:
            raise Exception(f"Unhandled gate {gate.name}")
    dic = {k: v * gate.count if v else 0 for k, v in res._asdict().items()}
    return TGateConversion(**dic)


def _convert_gate_x():
    return TGateConversion(H=2, T=4)


def _convert_gate_cCNOT1():
    return TGateConversion(H=2, T=7, CNOT=7)

def _convert_gate_cCNOT2():
    return TGateConversion(H=2, T=4, CNOT=8, add_qubits=1)


def _convert_gate_cswap1():
    return TGateConversion(H=2, T=7, CNOT=8)

def _convert_gate_cswap2():
    return TGateConversion(H=2, T=4, CNOT=10, add_qubits=1)


def _convert_gate_cz():
    return TGateConversion(H=2, CNOT=1)


def _convert_gate_crx():
    # - CRY(a) reg = CNOT reg, RY(-a/2) reg[1], CNOT reg, RY(a/2) reg[1]
    #   +  = CNOT reg, S, H, RZ(-a/2), H, S.dag() reg[1], CNOT reg, S, H, RZ(a/2), H, S.dag() reg[1]
    # - CRX(a) reg = S.dag() reg[1], CRY(-a), S
    # 2 RZ, 4 S, 2 CNOT, 4 H
    dic = _convert_gate_rz()._asdict()
    dic = {k: 2 * v for k, v in dic.items()}
    dic["S"] += 4
    dic["CNOT"] += 2
    dic["H"] += 4
    return TGateConversion(**dic)


def _convert_gate_ry():
    # 2 H, 2 S + RZ decomposition
    dic = _convert_gate_rz()._asdict()
    dic["H"] += 2
    dic["S"] += 2
    return TGateConversion(**dic)


def _convert_gate_rz():
    return TGateConversion(T=150, H=150, S=80)
