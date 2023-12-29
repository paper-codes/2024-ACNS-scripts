from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field, fields
from enum import Enum, IntEnum
import functools
from operator import add, mul
from typing import Callable, Dict, NamedTuple, Optional, Sequence

from sage.all import var
from sage.functions.all import log, sqrt, min_symbolic, max_symbolic
from sage.rings.all import RealField
from sage.symbolic.all import Expression
from sage.symbolic.constants import NaN

# NIST requirement parameters for depth
MAX_DEPTHS = (40, 64, 96)

import utils.datamanipulation as uda  # Must be after sage imports, that's why it's here


class SecurityLevels(IntEnum):
    AES_128 = 1
    AES_192 = 3
    AES_256 = 5


def get_qaes_nist(security_level: SecurityLevels):
    # The world is not ready for 3.10
    match security_level:
        case SecurityLevels.AES_128:
            return 170
        case SecurityLevels.AES_192:
            return 233
        case SecurityLevels.AES_256:
            return 298
        case _:
            raise AttributeError("Wrong security level %s" % str(security_level))


def get_qaes_zou(security_level: SecurityLevels):
    match security_level:
        case SecurityLevels.AES_128:
            return 159
        case SecurityLevels.AES_192:
            return 223
        case SecurityLevels.AES_256:
            return 288
        case _:
            raise AttributeError("Wrong security level %s" % str(security_level))


def get_qaes_langenberg(security_level: SecurityLevels):
    match security_level:
        case SecurityLevels.AES_128:
            # lb(1.55*2^86+1.19*2.86) + lb(1.16*2^81)
            return 167
        case SecurityLevels.AES_192:
            # lb(1.81*2^118 + 1.17*2^119) + lb(1.33*2^113)
            return 233
        case SecurityLevels.AES_256:
            # lb(1.41*2^151 + 1.83*2^151) + lb(1.57*2^145)
            return 298
        case _:
            raise AttributeError("Wrong security level %s" % str(security_level))


def get_qaes_jang22(security_level: SecurityLevels):
    match security_level:
        case SecurityLevels.AES_128:
            # lb(81312+800+12240) + lb(978) + 128
            return 154
        case SecurityLevels.AES_192:
            # lb(1.55*2^86+1.19*2.86) + lb(1.16*2^81)
            return 219
        case SecurityLevels.AES_256:
            # lb(113744 + 1103 + 17408) + lb(1377) + 256
            return 283
        case _:
            raise AttributeError("Wrong security level %s" % str(security_level))


def get_qaes_gates(
    security_level: SecurityLevels, maxdepth: int, version: int = 1
) -> int:
    match version:
        case 1:
            return get_qaes_nist(security_level) - maxdepth
        case 2:
            return get_qaes_zou(security_level) - maxdepth
        case 3:
            return get_qaes_langenberg(security_level) - maxdepth
        case _:
            raise AttributeError("Version %d not known" % version)


def get_n_parallel_instances(
    depth_iter: Expression, max_depth_log2: Expression, success_probability: Expression
):
    """Parallelize Grover's algorithm with a no. of instances dictated by max_depth, as per EB22
    depth_iter: depth of a single grover iteration
    max_depth: max depth accepted for a single quantum circuit
    success_probability: of the algorithm; usually no. of Grover's iterations is evaluated as sqrt(reciprocal(success_probability))
    """
    n_instances = (
        (depth_iter**2) * (1 / success_probability) / (2 ** (2 * max_depth_log2))
    )
    if n_instances < 1:
        return 1
    return n_instances


def get_n_parallel_instances_log2(
    depth_iter_log2: Expression,
    max_depth_log2: Expression,
    success_probability_log2: Expression,
):
    """Parallelize Grover's algorithm with a no. of instances dictated by max_depth, as per EB22
    depth_iter: depth of a single grover iteration
    max_depth: max depth accepted for a single quantum circuit
    success_probability: of the algorithm; usually no. of Grover's iterations is evaluated as sqrt(reciprocal(success_probability))
    """
    n_instances = 2 * depth_iter_log2 - success_probability_log2 - max_depth_log2 * 2
    if n_instances < 0:
        return 0
    return n_instances


class LinearCorrectingCodeNames(Enum):
    HQC = "HQC"
    BIKE = "BIKE"
    MCE = "McEliece"


# class CodeExpression(NamedTuple):
#     is_cyclic: bool


class Code(NamedTuple):
    codename: str
    level: SecurityLevels
    nn: int
    kk: int
    tt: int
    is_cyclic: bool
    # True for a message attack, False for a key recovery attack
    is_key: bool = False
    # max_depth_log2: int = infinity
    # set to true if you want also the condsider the real values
    # TODO change it, there should be a better way to design all of these classes
    init_real_values: bool = True


class CodeExtended(Code):
    @property
    def rr(self):
        return self.nn - self.kk


class CodeAES(NamedTuple):
    codename: str
    level: SecurityLevels
    kk: int


class StandardGate(NamedTuple):
    name: str
    count: Expression


@dataclass  # type: ignore[misc]
class ValueDict:
    grover_iters: Expression
    qubits: Expression

    def as_dict(self):
        return asdict(self)


@dataclass
class ValueCliffordTDict(ValueDict):
    """Measures are total, meaning it's not for single iteration"""

    t_count: Expression
    t_depth: Expression
    tq: Expression
    gates_total: Optional[Dict[str, Expression]] = None
    # depth considering both clifford and t at the same level
    overall_depth: Expression = None
    clifford_count: Expression = None


@dataclass
class ValueParallelized:
    # NIST gives 3 max_depths: 40, 64, 96
    max_depth: Optional[Expression]
    # value depending on the security level: 1-> 170, 3->233, 5-> 298
    qaes: Optional[Expression]
    # no. gates for aes is given by qaes/max_depth
    aes_gates: Optional[Expression]
    # depth^2/max_depth, w/ depth being the depth of the whole prange, sequential, circuit
    prange_depth: Optional[Expression]
    # depth^2/qaes
    security_margin: Optional[Expression]
    # All the values below depend on n_instances as well
    n_instances: Optional[Expression]
    depth: Optional[Expression] = NaN
    gates: Optional[Expression] = NaN
    gates_margin: Optional[Expression] = NaN

    def as_dict(self):
        return asdict(self)


@dataclass  # type: ignore[misc]
class ValueNormalDict(ValueDict):
    """_iter: single grover iteration for quantum stuff, or classical iter for
    classical stuff

    _total: total across all grover iterations

    _sum_total: sum of gates across all grover iterations

    _hybrid: total across all grover and classical iterations

    _sum_hybrid: sum of gates across all grover and classical iterations

    _mixed: mixing classical and quantum measures in arbitrary ways

    """

    n_grover_instances: Optional[Expression] = 1
    # gates per each iteration
    gates_iteration: Optional[Dict[str, Expression]] = None
    # sum of all the gates, considering an equal weight for all of them
    gates_sum_iter: Optional[Expression] = None
    # depth of single iteration
    depth_iteration: Optional[Expression] = None
    # list of gates multiplied by n iterations
    gates_total: Optional[Dict[str, Expression]] = None
    # total of gates multiplied by n iterations
    gates_sum_total: Optional[Expression] = None
    # depth iteration multiplied by n iterations
    depth_total: Optional[Expression] = None
    # depth_total * qubits
    dq: Optional[Expression] = None
    # No. of classical iterations
    classic_iters: Optional[Expression] = None
    # No. of classic gates per each classical iteration
    classic_gates_iter: Optional[Expression] = None
    # No. of classic gates across all classical iteration
    classic_gates_total: Optional[Expression] = None
    # amount of space require
    classic_space: Optional[Expression] = None
    # space * classic_gates_total
    classic_space_time: Optional[Expression] = None
    # quantum measures, taking also into account classic no. of iterations
    gates_total_hybrid: Optional[Dict[str, Expression]] = None
    gates_sum_hybrid: Optional[Expression] = None
    depth_total_hybrid: Optional[Expression] = None
    dq_hybrid: Optional[Expression] = None
    # total quantum + classical gates, considering having same weight
    gates_sum_hybrid_mixed: Optional[Expression] = None
    # quantum dq + classic space_time
    dq_hybrid_mixed: Optional[Expression] = None

    def _get_mult(self):
        raise Exception("Should be implemented in subclass")

    def _get_add(self):
        raise Exception("Should be implemented in subclass")

    def _get_n_iterations(self):
        if self.classic_iters:
            return self._get_mult()(self.classic_iters, self.grover_iters)
        return self.grover_iters

    def set_classic_gates_total(self):
        self.classic_gates_total = self._get_mult()(
            self.classic_gates_iter, self.classic_iters
        )

    def set_classic_space_time(self):
        self.classic_space_time = self._get_mult()(
            self.classic_gates_total, self.classic_space
        )

    def set_gates_sum_hybrid(self):
        self.gates_sum_hybrid = self._get_mult()(
            self.classic_iters, self.gates_sum_total
        )

    def set_gates_sum_hybrid_mixed(self):
        self.gates_sum_hybrid_mixed = self._get_add()(
            self.classic_gates_total, self.gates_sum_hybrid
        )

    # def set_dq_hybrid_mixed(self):

    def set_gates_total(self):
        if not self.gates_iteration:
            self.gates_total = None
            return
        gates1 = {}
        toadd1 = self.grover_iters
        gates2 = {}
        toadd2 = self._get_n_iterations()
        for k, v in self.gates_iteration.items():
            gates1[k] = self._get_mult()(v, toadd1)
            gates2[k] = self._get_mult()(v, toadd2)

        self.gates_total = gates1
        self.gates_total_hybrid = gates2

    def set_depth_total(self):
        if not self.depth_iteration:
            self.depth_total = 0
            return
        self.depth_total = self._get_mult()(self.depth_iteration, self.grover_iters)
        self.depth_total_hybrid = self._get_mult()(
            self.depth_iteration, self._get_n_iterations()
        )

    def set_dq_total(self):
        if not self.depth_total or not self.qubits:
            self.dq = 0
            return
        self.dq = self._get_mult()(self.depth_total, self.qubits)
        self.dq_hybrid = self._get_mult()(self.depth_total_hybrid, self.qubits)

    def set_all_computable(self):
        if self.grover_iters and self.grover_iters != float("inf"):
            if not self.gates_total:
                self.set_gates_total()
            if not self.depth_total or not self.depth_total_hybrid:
                self.set_depth_total()
            if not self.gates_sum_total:
                if self.gates_sum_iter and self.grover_iters:
                    self.gates_sum_total = self._get_mult()(self.gates_sum_iter, self.grover_iters)
            if not self.dq or not self.dq_hybrid:
                self.set_dq_total()
                # self.dq = self._get_mult()(self.depth_total, self.qubits)
                # self.dq_hybrid = self._get_mult()(self.depth_total_hybrid, self.qubits)
        else:
            self.grover_iters = 0
            self.n_grover_instances = 0
            self.gates_total = 0
            self.gates_total_hybrid = 0
            self.depth_total = 0
            self.depth_total_hybrid = 0
            self.gates_sum_iter = 0
            self.gates_sum_total = 0
            self.dq = 0
            self.dq_hybrid = 0

        if self.classic_iters and self.classic_gates_iter:
            if not self.classic_gates_total:
                self.set_classic_gates_total()
            if not self.classic_space_time:
                self.set_classic_space_time()
        else:
            # at least one !?
            self.classic_iters = 1
            self.classic_gates_iter = 0
            self.classic_gates_total = 0
            self.classic_space_time = 0

        if not self.gates_sum_hybrid:
            self.set_gates_sum_hybrid()
        if not self.gates_sum_hybrid_mixed:
            self.set_gates_sum_hybrid_mixed()
        if not self.dq_hybrid_mixed:
            self.dq_hybrid_mixed = self._get_add()(self.dq_hybrid, self.classic_space_time)


@dataclass
class ValueNormalDictLog2(ValueNormalDict):
    def _get_mult(self) -> Callable[[Expression, Expression], Expression]:
        return add

    def _get_add(self) -> Callable[[Expression, Expression], Expression]:
        # approximation of a sum between two log2 values. That is
        # a_2 + b_2 = max(a_2, b_2) + 1 - min(|a_2 - b_2|, 1)
        return functools.partial(
            lambda a, b: max_symbolic(a, b) + 1 - min_symbolic(abs(a - b), 1)
        )


@dataclass
class ValueNormalDictInt(ValueNormalDict):
    def _get_mult(self) -> Callable[[Expression, Expression], Expression]:
        return mul

    def _get_add(self) -> Callable[[Expression, Expression], Expression]:
        return add

    def set_all_computable(self):
        # These two make only sense with normal; with log2 not really
        self.set_gates_sum_iter()
        self.set_gates_sum_total()
        super().set_all_computable()

    def set_gates_sum_iter(self):
        if self.gates_iteration:
            self._set_gates_sum(False)

    def set_gates_sum_total(self):
        if self.gates_iteration:
            self._set_gates_sum(True)

    def _set_gates_sum(self, total):
        total_gates = sum(self.gates_iteration.values())
        if total:
            total_gates *= self.grover_iters
            self.gates_sum_total = total_gates
        else:
            self.gates_sum_iter = total_gates

    def to_log2(
        self, skip_grover_iters: bool = False, skip_classic_iters: bool = False
    ) -> ValueNormalDictLog2:
        """Transform each value in its log2
        # skip_grov_iters is sometimes needed since sage has problems with log2 of binomials
        keywords:
        """
        qubits = log(self.qubits, 2) if self.qubits else None
        if self.grover_iters and not skip_grover_iters:
            grov_iters = log(self.grover_iters, 2)
        else:
            grov_iters = None

        vdict = ValueNormalDictLog2(grover_iters=grov_iters, qubits=qubits)

        for f in fields(self):
            f_name = f.name
            f_value = getattr(self, f_name)
            if f_value is None:
                continue
            # already set before
            if f_name in ("qubits", "grover_iters"):
                continue
            if skip_classic_iters and f_name == "classic_iters":
                continue
            if f_value == 0:
                continue

            if f_name in ("gates_iteration", "gates_total", "gates_total_hybrid"):
                gates = {}
                # field value is a dic
                for k, v in f_value.items():
                    gates[k] = log(v, 2)
                setattr(vdict, f_name, gates)
            else:
                setattr(vdict, f_name, log(f_value, 2))
        return vdict


@dataclass
class ValueDicts:
    code_ext: Optional[CodeExtended | CodeAES]
    in_params: Optional[Dict[str, Expression]] = None
    normal: Optional[ValueNormalDictInt] = None
    normal2: Optional[ValueNormalDictLog2] = None
    normal_expanded: Optional[ValueNormalDictInt] = None
    normal_expanded2: Optional[ValueNormalDictLog2] = None
    tmeas: Optional[ValueCliffordTDict] = None
    tmeas2: Optional[ValueCliffordTDict] = None
    # Generic, not depending on the security level
    # parallel_expanded2_depth: Optional[float] = None
    # This is generic
    parallel_expanded2_maxd: Optional[ValueParallelized] = None
    # This is supposed to hold actual values
    parallel_expanded2_maxds: Dict[int, ValueParallelized] = field(default_factory=dict)

    def as_dict(self):
        data = asdict(self)
        return {key: value for key, value in data.items() if value is not None}


@dataclass
class CodeExtendedSage(ABC):
    REAL_PREC = 1000
    REAL_FIELD = RealField(prec=REAL_PREC)

    def __init__(self, code_ext: CodeExtended, **kwords) -> None:
        # self.code = CodeExtended(codename, level, nn, kk, tt)
        self.code_ext = code_ext
        # o stands for the original values
        n_o, k_o, t_o, r_o = var("n_o, k_o, t_o, r_o", domain="integer")
        # sec_level = var("sec_level", domain="integer")
        max_depth2 = var("max_depth2", domain="integer")
        # # quantum gates estimate for aes, in log2
        qaes2 = var("qaes2", domain="integer")
        # Dict associting to each string a corresponding name
        self.var_by_names = {str(v): v for v in (n_o, k_o, t_o, r_o, max_depth2, qaes2)}
        self.input_params = {str(v): v for v in (n_o, k_o, t_o, r_o)}

        # Dictionary of variables with their actual corresponding values, used
        # when we want to replace the expressions with the actual value
        if code_ext.init_real_values:
            self.var_subs = {
                self.var_by_names["n_o"]: self.REAL_FIELD(code_ext.nn),
                self.var_by_names["k_o"]: self.REAL_FIELD(code_ext.kk),
                self.var_by_names["t_o"]: self.REAL_FIELD(code_ext.tt),
                self.var_by_names["r_o"]: self.REAL_FIELD(code_ext.rr),
                self.var_by_names["r_o"]: self.REAL_FIELD(code_ext.rr),
            }

        self._init_additional(**kwords)

    def _init_additional(self, **kwords):
        pass

    @abstractmethod
    def go(self, **kwords) -> ValueDicts:
        pass

    def is_cyclic(self):
        return self.code_ext.is_cyclic

    def get_dict_input_params(self, subs: bool) -> Dict[str, Expression]:
        if not subs:
            return self.var_by_names
        if hasattr(self, "input_params"):
            to_check = self.input_params
        else:
            to_check = self.var_by_names
        return uda.replace_exp_dict(to_check, self.var_subs, numerical=True)

    @classmethod
    def get_grover_iterations(
        cls, domain_size: Expression, n_solutions: Expression, red_factor=0.58278
    ) -> Expression:
        return red_factor * sqrt(domain_size / n_solutions)

    @classmethod
    def get_grover_iterations_log2(
        cls,
        numerators: Sequence[Expression],
        denominators: Sequence[Expression],
        red_factor: Expression = 0.58278,
    ) -> Expression:
        iters = log(red_factor, 2)
        for exp in numerators:
            iters += log(exp, 2) / 2
        for exp in denominators:
            iters -= log(exp, 2) / 2
        return iters
