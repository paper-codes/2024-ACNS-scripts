import os
import re
import pickle
from typing import Any, Callable, Dict, Iterator, TYPE_CHECKING, Union

if TYPE_CHECKING:
    from measures.common import ValueDicts  # , Code, SecurityLevels

# Lee-Brickell
OUT_FILES_LB_PART_FMT = (
    "{code_name:*^20}_L{code_level}_n{n:06}_k{k:06}_t{t:03}_p{p:03}.{ext}"
)
OUT_FILES_LB_DIR: str = os.path.join(".", "out", "{variant}", "{out_type}")
OUT_FILES_LB_FMT: str = os.path.join(OUT_FILES_LB_DIR, OUT_FILES_LB_PART_FMT)

# Prange
OUT_FILES_PRANGE_DIR: str = os.path.join(
    ".", "out", "{variant}", "{expansion}", "{out_type}"
)
OUT_FILES_DJB_DIR: str = os.path.join(".", "out", "{variant}", "{out_type}")
OUT_FILES_ESSER_DIR: str = OUT_FILES_DJB_DIR
OUT_FILES_EB22_DIR: str = OUT_FILES_DJB_DIR
OUT_FILES_EB22_ESTIMATOR_DIR: str = OUT_FILES_DJB_DIR
OUT_FILES_PRANGE_ESTIMATOR_DIR: str = os.path.join(
    ".", "out", "{variant}", "{out_type}"
)
OUT_FILES_K16_DIR: str = OUT_FILES_DJB_DIR
OUT_FILES_AES_DIR: str = os.path.join(".", "out", "{variant}", "{out_type}")
#
OUT_FILES_PART_FMT: str = "{code_name:*^20}_L{code_level}_n{n:06}_k{k:06}_t{t:03}.{ext}"
OUT_FILES_PRANGE_FMT: str = os.path.join(OUT_FILES_PRANGE_DIR, OUT_FILES_PART_FMT)
OUT_FILES_DJB_FMT: str = os.path.join(OUT_FILES_DJB_DIR, OUT_FILES_PART_FMT)
OUT_FILES_ESSER_FMT: str = os.path.join(OUT_FILES_ESSER_DIR, OUT_FILES_PART_FMT)
OUT_FILES_EB22_FMT: str = os.path.join(OUT_FILES_EB22_DIR, OUT_FILES_PART_FMT)
OUT_FILES_K16_FMT: str = os.path.join(OUT_FILES_K16_DIR, OUT_FILES_PART_FMT)
# Expression file
OUT_FILES_EXP_PRANGE_FMT: str = os.path.join(OUT_FILES_PRANGE_DIR, "0_{cyclic}.{ext}")
OUT_FILES_EXP_EB22_FMT: str = os.path.join(OUT_FILES_EB22_DIR, "0_{cyclic}.{ext}")
# Estimator file (variable n)
OUT_FILES_PART_ESTIMATOR_FMT: str = (
    "{rate:.2f}_n{n:10}_k{k:10}_t{weight_fun_name}_{t:10}.{ext}"
)
OUT_FILES_PRANGE_ESTIMATOR_FMT: str = os.path.join(
    OUT_FILES_PRANGE_ESTIMATOR_DIR, OUT_FILES_PART_ESTIMATOR_FMT
)
OUT_FILES_EB22_ESTIMATOR_FMT: str = os.path.join(
    OUT_FILES_EB22_ESTIMATOR_DIR, OUT_FILES_PART_ESTIMATOR_FMT
)
#
OUT_FILES_PART_AES_FMT: str = "{code_name:*^10}_L{code_level}.{ext}"
OUT_FILES_AES_FMT: str = os.path.join(OUT_FILES_AES_DIR, OUT_FILES_PART_AES_FMT)
####
# Figures
FIG_DIR: str = os.path.join(".", "figs")
####


def save_dict_to_txt_file(
    filename: str, dictionary: Union["ValueDicts", Dict[str, Any]], mode: str
):
    with open(filename, mode, encoding="utf-8") as fp:
        dic = dictionary if isinstance(dictionary, dict) else dictionary.as_dict()
        for k, v in dic.items():
            print(">", file=fp)
            print(k, file=fp)
            print(v, file=fp)
        print("*" * 50, file=fp)


def save_obj_to_obj_file(filename: str, obj: Any):
    with open(filename, "wb") as fp:
        pickle.dump(obj, fp)


def save_vdict_to_obj_file(filename: str, dictionary: "ValueDicts"):
    save_obj_to_obj_file(filename, dictionary)


def load_obj_from_obj_file(filename: str) -> Any:
    with open(filename, "rb") as fp:
        return pickle.load(fp)


def load_vdict_from_obj_file(filename: str) -> "ValueDicts":
    return load_obj_from_obj_file(filename)


def get_filtered_files_of_dir(
    path: str,
    filter_path: Callable[[str], bool] = lambda _: True,
    keysort: Callable[[str], str] = lambda x: re.sub(r"[^A-Za-z0-9]+", "", x),
) -> Iterator[str]:
    allfiles = filter(
        filter_path,
        [f"{path}/{f}" for f in sorted(os.listdir(path), key=keysort)],
    )
    return allfiles


# def get_code_by_fname(fname: str) -> Code:
#     """returns code, level, n, k, t"""
#     pos = fname.rfind('/')
#     end_pos = fname.find('_')
#     code = fname[pos + 1:end_pos]
#     code = code.replace('*', '')
#     pos = fname.find('L', end_pos)
#     end_pos = fname.find('_', pos)
#     level = fname[pos + 1:end_pos]
#     pos = fname.find('n', end_pos)
#     end_pos = pos + 7
#     ns = fname[pos + 1:end_pos]
#     n = int(ns)
#     pos = fname.find('k', end_pos + 1)
#     end_pos = pos + 7
#     ks = fname[pos + 1:end_pos]
#     k = int(ks)
#     pos = fname.find('t', end_pos)
#     end_pos = pos + 4
#     ts = fname[pos + 1:end_pos]
#     t = int(ts)
#     return Code(code, SecurityLevels(int(level)), n, k, t)
