import argparse
import os
from typing import Type, Optional

from statics.codes import CODES
from utils.export.file_mgmt import (
    OUT_FILES_DJB_DIR,
    OUT_FILES_DJB_FMT,
    OUT_FILES_ESSER_DIR,
    OUT_FILES_ESSER_FMT,
    OUT_FILES_EXP_PRANGE_FMT,
    OUT_FILES_PRANGE_DIR,
    OUT_FILES_PRANGE_FMT,
    save_dict_to_txt_file,
    save_vdict_to_obj_file,
)

from measures.common import CodeExtended, CodeExtendedSage
from measures.esser import Esser
from measures.prange import Prange

from math import inf

# from measures.zou import Zou


def parse_arguments():
    parser = argparse.ArgumentParser("Counts the zeros of the RREF reversible function")
    parser.add_argument(
        "launcher",
        choices=["prange", "esser"],
        default="prange",
    )
    parser.add_argument(
        "--optimization",
        type=int,
        choices=(0, 1, 2, 3, 4),
        default=0,
        # Check A. Esser, S. Ramos-Calderer, E. Bellini, J. I. Latorre, and M.
        # Manzano, “An Optimized Quantum Implementation of ISD on Scalable
        # Quantum Resources,” Cryptology ePrint Archive, 2021.
        help=("Optimizations: 0 no opt, 1 hybrid, 2 puncturing, 3 both, 4 column"),
    )
    parser.add_argument("--alpha", type=float, help="required for 1 and 3 optimization")
    parser.add_argument("--beta", type=float, help="required for 2 and 3 optimization")
    parser.add_argument("--delta-var", action="store_true")
    parser.add_argument(
        "--delta",
        type=float,
        help="required for 1, 2 and 3 optimization instead of alpha or delta",
    )
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--out-format", choices=["txt", "bin"], default="txt")
    parser.add_argument("--expansion", choices=["ours", "saeedi"], default="ours")
    parser.add_argument("--aes-diff", action="store_true")
    parser.add_argument("--skip-code-system", nargs="+", default=["leda"])
    parser.add_argument(
        "--formula-only",
        action="store_true",
        help="Just store the formula without computing the actual values",
    )
    return parser


def main(raw_args: Optional[list[str]] = None):
    print("#" * 80)
    parser = parse_arguments()
    if raw_args and len(raw_args) != 0:
        namespace = parser.parse_args(raw_args)
    else:
        namespace = parser.parse_args()
    print(namespace)

    kwords = {}
    launch: Type[CodeExtendedSage]
    if namespace.launcher == "prange":
        launch = Prange
        kwords["alphaa"] = namespace.alpha
        kwords["deltad"] = namespace.delta
        kwords["optimization"] = namespace.optimization
        kwords["expansion"] = namespace.expansion
        kwords["subs"] = not namespace.formula_only

        out_dirs = OUT_FILES_PRANGE_DIR.format(
            variant=namespace.launcher,
            out_type=namespace.out_format,
            expansion=namespace.expansion,
        )
        if namespace.formula_only:
            out_file = OUT_FILES_EXP_PRANGE_FMT
        else:
            out_file = OUT_FILES_PRANGE_FMT
    elif namespace.launcher == "esser":
        launch = Esser
        out_dirs = OUT_FILES_ESSER_DIR.format(
            variant=namespace.launcher, out_type=namespace.out_format
        )
        out_file = OUT_FILES_ESSER_FMT
    else:
        raise Exception(f"wrong launcher {namespace.launcher}")

    if namespace.out_format == "txt":
        ext = "txt"
    elif namespace.out_format == "bin":
        ext = "pkl"
    else:
        raise Exception("Unrecognized out format")

    os.makedirs(out_dirs, exist_ok=True)

    # To get only the parametric expressions
    if namespace.launcher == "prange" and namespace.formula_only:
        for is_cyclic in True, False:
            filename = out_file.format(
                variant=namespace.launcher,
                out_type=namespace.out_format,
                expansion=namespace.expansion,
                ext=ext,
                cyclic=is_cyclic,
            )
            meas = launch(
                CodeExtended(
                    "", "", inf, inf, inf, is_cyclic, init_real_values=False
                )
            )
            res = meas.go(**kwords)
            if namespace.out_format == "txt":
                save_dict_to_txt_file(filename, res, "w")
            elif namespace.out_format == "bin":
                save_vdict_to_obj_file(filename, res)
        return

    codes2 = []
    for code in CODES:
        found = False
        for skipcode in namespace.skip_code_system:
            if code.codename.lower().find(skipcode.lower()) >= 0:
                found = True
                break
        if not found:
            codes2.append(code)
    for code in codes2:

        filename = out_file.format(
            variant=namespace.launcher,
            out_type=namespace.out_format,
            expansion=namespace.expansion,
            code_name=code.codename,
            code_level=code.level,
            n=code.nn,
            k=code.kk,
            t=code.tt,
            ext=ext,
        )
        print(filename)
        if namespace.skip_existing and os.path.isfile(filename):
            print("skipping")
            continue

        # if namespace.out_format == 'txt':
        #     save_dict_to_txt_file(filename, {'namespace': namespace}, 'w')

        tm_min = float("inf")
        delta_min = float("inf")
        if namespace.delta_var and namespace.optimization > 0:
            # delta varying from 1 to 0 in .05 steps
            for delta in range(100, -1, -5):
                namespace.delta = delta / 100
                meas = launch(
                    CodeExtended(
                        code.codename,
                        code.level,
                        code.nn,
                        code.kk,
                        code.tt,
                        code.is_cyclic,
                    )
                )
                res = meas.go(**kwords)
                # print(f"alpha {alpha} is greater than or equal to R {kk/nn}")
                if len(res) == 0:
                    # print(f"alpha {alpha} >= R {kk/nn}")
                    break
                tm = res["T*M"]
                if tm < tm_min:
                    tm_min = tm
                    delta_min = namespace.delta
                if namespace.out_format == "txt":
                    save_dict_to_txt_file(filename, res, "w")
            if namespace.out_format == "txt":
                res = {}
                res["tm min"] = tm_min
                res["delta min"] = delta_min
                save_dict_to_txt_file(filename, res, "a")
        else:
            meas = launch(
                CodeExtended(
                    code.codename, code.level, code.nn, code.kk, code.tt, code.is_cyclic
                )
            )
            res = meas.go(**kwords)
            if namespace.out_format == "txt":
                save_dict_to_txt_file(filename, res, "w")
            elif namespace.out_format == "bin":
                save_vdict_to_obj_file(filename, res)

    print("#" * 80)


if __name__ == "__main__":
    main()
