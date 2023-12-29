"""Actual measures from Esser+21
"""

from sage.all import sqrt, var
from sage.functions.all import binomial, log
from utils.datamanipulation import replace_exp_valuedict

from measures.common import CodeExtendedSage, ValueNormalDictInt, ValueDicts

class Esser(CodeExtendedSage):
    def _init_additional(self, **kwords):
        p_o = var("p_o", domain="integer")

        self.var_by_names["p_o"] = p_o
        self.input_params["p_o"] = p_o
        self.var_subs[p_o] = kwords["p"]

    def _dicke(self, running_dic: ValueNormalDictInt):
        n = self.var_by_names['n_o']
        k = self.var_by_names['k_o']
        r = n - k
        # running_dic.qubits += n + 2 * log(r + 1, 2)
        running_dic.depth_iteration += n * r * log(r, 2) * log(log(r, 2), 2)

    def _permutation(self, running_dic: ValueNormalDictInt):
        n = self.var_by_names['n_o']
        k = self.var_by_names['k_o']
        running_dic.depth_iteration += (n - k) * (n**2)

    def _gje(self, running_dic: ValueNormalDictInt):
        n = self.var_by_names['n_o']
        k = self.var_by_names['k_o']
        t = self.var_by_names['t_o']
        # running_dic.qubits += n^2 + n - 2
        running_dic.depth_iteration += (n**3) * log(n, 2)

    def _exhaustive_enumeration(self, running_dic: ValueNormalDictInt):
        n = self.var_by_names['n_o']
        k = self.var_by_names['k_o']
        r = n - k
        t = self.var_by_names['t_o']
        p = self.var_by_names['p_o']
        # running_dic.qubits += n^2 + n - 2
        running_dic.depth_iteration += binomial(k, p) * (r * 2 + r * 2)

    def _get_m(self):
        return self.var_by_names['r_o'] if self.is_cyclic() else 1

    def go(self, **kwords) -> 'ValueDicts':
        n = self.var_by_names['n_o']
        k = self.var_by_names['k_o']
        r = n - k
        t = self.var_by_names['t_o']

        ret_dic = ValueDicts(code_ext=self.code_ext)
        ret_dic.in_params = self.get_dict_input_params(subs=True)
        M = self._get_m()
        grover_iters = self.get_grover_iterations(binomial(
            n, r), (0.2888 * binomial(n - t, k)) * M)

        qubits = (n - k + 1) * (n + 2) - 4 + M
        dic = ValueNormalDictInt(
            grover_iters=grover_iters,
            qubits=qubits,
            depth_iteration=0,
            gates_iteration=None,
            gates_sum_iter=0,
            gates_sum_total=0,
        )
        self._dicke(dic)
        self._permutation(dic)
        self._gje(dic)
        self._exhaustive_enumeration(dic)

        # dic.set_gates_sum_iter()
        # dic.set_gates_sum_total()
        dic.depth_iteration = dic.depth_iteration * 2

        ret_dic.normal_expanded = replace_exp_valuedict(dic,
                                                        self.var_subs,
                                                        numerical=True)

        dic2 = dic.to_log2()
        dic2.set_depth_total()
        dic2.set_dq_total()
        ret_dic.normal_expanded2 = replace_exp_valuedict(dic2,
                                                         self.var_subs,
                                                         numerical=True)
        return ret_dic
