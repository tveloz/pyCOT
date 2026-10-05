"""
run_all_eso_eo.py -- run the ESO/EO-aware compare_endosymbiosis on all four
endosymbiosis cases (toy, mitochondrial-type, chloroplast-type dark/light),
saving each log to outputs/. Consolidates what were previously four
separate one-off invocations into a single reproducible script.
"""
from __future__ import annotations

import os
import sys
import contextlib
import io

_here = os.path.dirname(os.path.abspath(__file__))
_root = os.path.join(_here, '..')
_out_dir = os.path.join(_root, 'outputs')
os.makedirs(_out_dir, exist_ok=True)
sys.path.insert(0, _here)

from endosymbiosis_analysis import compare_endosymbiosis


def run_and_save(label, log_name, **kwargs):
    print(f"\n########## {label} ##########", flush=True)
    buf = io.StringIO()

    class Tee(io.StringIO):
        def write(self_inner, s):
            sys.__stdout__.write(s)
            sys.__stdout__.flush()
            buf.write(s)
            return len(s)

    tee = Tee()
    with contextlib.redirect_stdout(tee):
        compare_endosymbiosis(**kwargs)
    out_path = os.path.join(_out_dir, log_name)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(buf.getvalue())
    print(f"[saved] {out_path}", flush=True)


def main():
    toy = os.path.join(_root, 'toy_model')
    run_and_save(
        "Toy model", "toy_eso_eo_analysis.log",
        host_path=os.path.join(toy, 'host_alone.txt'),
        symbiont_path=os.path.join(toy, 'symbiont_alone.txt'),
        merged_path=os.path.join(toy, 'merged.txt'),
        host_suffix="_h", symbiont_suffix="__endo",
    )

    rd = os.path.join(_root, 'real_data')
    run_and_save(
        "Mitochondrial-type", "mito_eso_eo_analysis.log",
        host_path=os.path.join(rd, 'mito_host_fermentative.txt'),
        symbiont_path=os.path.join(rd, 'mito_symbiont_aerobic.txt'),
        merged_path=os.path.join(rd, 'mito_merged.txt'),
        symbiont_suffix="__endo",
    )

    run_and_save(
        "Chloroplast-type, dark", "chloro_dark_eso_eo_analysis.log",
        host_path=os.path.join(rd, 'chloro_host_yeast.txt'),
        symbiont_path=os.path.join(rd, 'chloro_symbiont_synecho_dark.txt'),
        merged_path=os.path.join(rd, 'chloro_merged_dark.txt'),
        symbiont_suffix="__endo",
    )

    run_and_save(
        "Chloroplast-type, +light", "chloro_light_eso_eo_analysis.log",
        host_path=os.path.join(rd, 'chloro_host_yeast.txt'),
        symbiont_path=os.path.join(rd, 'chloro_symbiont_synecho_light.txt'),
        merged_path=os.path.join(rd, 'chloro_merged_light.txt'),
        symbiont_suffix="__endo",
    )


if __name__ == "__main__":
    main()
