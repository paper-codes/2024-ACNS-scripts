import argparse
import os
from typing import Optional, Type

from statics.codes import CODES
from utils.export.file_mgmt import (
    OUT_FILES_LB_DIR,
    OUT_FILES_LB_FMT,
    save_dict_to_txt_file,
    save_vdict_to_obj_file,
)

from measures.common import CodeExtended, CodeExtendedSage
from measures.leebrickell_classic import LeeBrickellClassic
from measures.leebrickell_hybrid import LeeBrickellHybrid
from measures.leebrickell_quantum import LeeBrickellQuantum
from measures.esser_lb import Esser


def parse_arguments():
    parser = argparse.ArgumentParser("Launch Lee-Brickell")
    parser.add_argument(
        "launcher",
        choices=["lee-brickell-2", "lee-brickell-4", "lee-brickell-classic", "esser21-lb"],
        default="lee-brickell-4",
        help=(
            "Variants 1, 2, ecc. refer to the variants numbering in the TQC23 paper."
        ),
    )
    parser.add_argument("-p", type=int, required=True)
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--out-format", choices=["txt", "bin"], default="txt")
    parser.add_argument("--skip-code-system", nargs="+", default=[])
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
    kwords["subs"] = not namespace.formula_only
    launch: Type[CodeExtendedSage]
    if namespace.launcher == "lee-brickell-4":
        launch = LeeBrickellHybrid
        out_dirs = OUT_FILES_LB_DIR.format(
            variant=namespace.launcher,
            out_type=namespace.out_format,
        )
        if namespace.formula_only:
            raise NotImplementedError("still not implemented")
        else:
            out_file = OUT_FILES_LB_FMT
    elif namespace.launcher == "lee-brickell-2":
        launch = LeeBrickellQuantum
        out_dirs = OUT_FILES_LB_DIR.format(
            variant=namespace.launcher,
            out_type=namespace.out_format,
        )
        if namespace.formula_only:
            raise NotImplementedError("still not implemented")
        else:
            out_file = OUT_FILES_LB_FMT
    elif namespace.launcher == "lee-brickell-classic":
        launch = LeeBrickellClassic
        out_dirs = OUT_FILES_LB_DIR.format(
            variant=namespace.launcher,
            out_type=namespace.out_format,
        )
        if namespace.formula_only:
            raise NotImplementedError("still not implemented")
        else:
            out_file = OUT_FILES_LB_FMT
    elif namespace.launcher == "esser21-lb":
        launch = Esser
        out_dirs = OUT_FILES_LB_DIR.format(
            variant=namespace.launcher,
            out_type=namespace.out_format,
        )
        if namespace.formula_only:
            raise NotImplementedError("still not implemented")
        else:
            out_file = OUT_FILES_LB_FMT
    else:
        raise ValueError(f"wrong launcher {namespace.launcher}")

    if namespace.out_format == "txt":
        ext = "txt"
    elif namespace.out_format == "bin":
        ext = "pkl"
    else:
        raise ValueError("Unrecognized out format")

    os.makedirs(out_dirs, exist_ok=True)

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
            code_name=code.codename,
            code_level=code.level,
            n=code.nn,
            k=code.kk,
            t=code.tt,
            p=namespace.p,
            ext=ext,
        )
        print(filename)
        if namespace.skip_existing and os.path.isfile(filename):
            print("skipping")
            continue

        # if namespace.out_format == 'txt':
        #     save_dict_to_txt_file(filename, {'namespace': namespace}, 'w')

        meas = launch(
            CodeExtended(
                code.codename, code.level, code.nn, code.kk, code.tt, code.is_cyclic
            ),
            p=namespace.p,
        )
        res = meas.go(**kwords)
        if namespace.out_format == "txt":
            save_dict_to_txt_file(filename, res, "w")
        elif namespace.out_format == "bin":
            save_vdict_to_obj_file(filename, res)

    print("#" * 80)


if __name__ == "__main__":
    main()
