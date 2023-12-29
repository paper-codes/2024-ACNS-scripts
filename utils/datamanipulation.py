from typing import TYPE_CHECKING, Any, Dict, Iterator, Tuple, Union

from sage.functions.all import log
from sage.symbolic.all import Expression

from utils.typing import TYPE_DICT_ITERATOR, TYPE_DICT_VALUES

if TYPE_CHECKING:
    from measures.common import ValueDict


def get_log2_dict(dic: TYPE_DICT_VALUES):
    dicl2 = {}
    for k, v in dic.items():
        if isinstance(v, (Expression, float, int)):
            dicl2[k] = log(v, 2)
            # elif isinstance(v, dict):
            #     dicl2[k] = get_log2_dict(v)
            # elif isinstance(v, str):
            #     dicl2[k] = v
        else:
            raise Exception(f"Undefined value {type(v)}")
    return dicl2


def replace_exp_valuedict(
    value_dict: "ValueDict", subs_dic: Dict[Expression, Expression], numerical=False
) -> "ValueDict":
    dic = replace_exp_dict(value_dict.as_dict(), subs_dic, numerical)
    return value_dict.__class__(**dic)


def replace_exp_dict(
    dic: TYPE_DICT_VALUES, subs_dic: Dict[Expression, Expression], numerical=False
) -> dict:
    if subs_dic is None or len(subs_dic) == 0:
        print("None")
        return dic

    return replace_exp_dictitemiter(iter(dic.items()), subs_dic, numerical=numerical)


def replace_exp_dictitemiter(
    result_exp_iter: TYPE_DICT_ITERATOR,
    subs_dic: Dict[Expression, Expression],
    numerical=False,
):
    return _replace_exp_dictitemiter_support(
        result_exp_iter, subs_dic, numerical=numerical
    )


def _replace_exp_dictitemiter_support(
    result_exp_iter: Iterator[Tuple[str, Expression]],
    subs_dic: Dict[Expression, Expression],
    numerical=False,
):
    result_sub: Dict[str, Union[Expression, str, Dict[str, Any]]] = {}
    for k, v in result_exp_iter:
        _replace_dictitemiter_elem(k, v, subs_dic, numerical, result_sub)
    return result_sub


def _replace_dictitemiter_elem(k, v, subs_dic, numerical, result_sub, keep_negative=True):
    if isinstance(v, dict):
        result_sub[k] = _replace_exp_dictitemiter_support(
            iter(v.items()), subs_dic, numerical
        )
    elif isinstance(v, Expression):
        exp = replace_expression(v, subs_dic, numerical)
        if exp > 0 or (exp < 0 and keep_negative):
            result_sub[k] = exp
        else:
            # TODO check, we are replacing negative values with 0
            result_sub[k] = 0
    else:
        result_sub[k] = str(v)


def replace_expression(
    v: Expression, subs_dic: Dict[Expression, Expression], numerical: bool
):
    exp = v.subs(subs_dic)
    interset = set(exp.variables()).intersection(set(subs_dic))
    while len(interset) > 0:
        exp = exp.subs(subs_dic)
        interset = set(exp.variables()).intersection(set(subs_dic))
    # assert len(interset) == 0, f"{k}: {v}, interset {interset}"
    if numerical:
        exp = exp.unhold().n()
    return exp


def merge_dicts_by_keys(*dics: TYPE_DICT_VALUES, exclude_keys=()) -> dict:
    # Can handle also symbolic computation
    return merge_dicts_by_keys_mul_factor(1, *dics, exclude_keys=exclude_keys)


def merge_dicts_by_keys_mul_factor(
    mul_factor: Expression, *dics: TYPE_DICT_VALUES, exclude_keys=()
):
    # Can handle also symbolic computation
    all_ks = set().union(*(d.keys() for d in dics))
    all_ks_filtered = filter(lambda x: x not in exclude_keys, all_ks)
    return {k: sum(t.get(k, 0) * mul_factor for t in dics) for k in all_ks_filtered}


def aggregate_gate_counts(mul_factor: Expression, *dics: TYPE_DICT_VALUES):
    # Assuming all gates count as 1
    overall = merge_dicts_by_keys_mul_factor(mul_factor, *dics)
    # return sum(overall.values()), overall
    return overall
