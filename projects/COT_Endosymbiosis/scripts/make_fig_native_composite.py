"""
make_fig_native_composite.py -- combine the three native-pyCOT-rendered
network PNGs (host, symbiont, merged) into one composite figure for the
PNAS paper: host and symbiont side by side on top, merged network below
(cropped/scaled to fit), with panel labels (a)/(b)/(c).
"""
from __future__ import annotations

import os
from PIL import Image, ImageDraw, ImageFont

_here = os.path.dirname(os.path.abspath(__file__))
_fig_dir = os.path.join(_here, '..', 'figures')

host = Image.open(os.path.join(_fig_dir, 'native_host.png')).convert("RGB")
symb = Image.open(os.path.join(_fig_dir, 'native_symbiont.png')).convert("RGB")
merged = Image.open(os.path.join(_fig_dir, 'native_merged.png')).convert("RGB")

PAD = 30
LABEL_H = 44


def add_label(img, text):
    w, h = img.size
    canvas = Image.new("RGB", (w, h + LABEL_H), "white")
    canvas.paste(img, (0, LABEL_H))
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("arialbd.ttf", 32)
    except Exception:
        font = ImageFont.load_default()
    draw.text((6, 4), text, fill="black", font=font)
    return canvas


host_l = add_label(host, "(a) host alone")
symb_l = add_label(symb, "(b) symbiont alone")
merged_l = add_label(merged, "(c) merged (interface highlighted)")

# Top row: host + symbiont, matched height.
top_h = max(host_l.height, symb_l.height)


def resize_to_h(img, h):
    w = int(img.width * h / img.height)
    return img.resize((w, h), Image.LANCZOS)


host_r = resize_to_h(host_l, top_h)
symb_r = resize_to_h(symb_l, top_h)

top_w = host_r.width + PAD + symb_r.width
top_row = Image.new("RGB", (top_w, top_h), "white")
top_row.paste(host_r, (0, 0))
top_row.paste(symb_r, (host_r.width + PAD, 0))

# Bottom row: merged, scaled so its width matches the top row's width.
merged_r = merged_l.resize((top_w, int(merged_l.height * top_w / merged_l.width)), Image.LANCZOS)

total_h = top_h + PAD + merged_r.height
canvas = Image.new("RGB", (top_w, total_h), "white")
canvas.paste(top_row, (0, 0))
canvas.paste(merged_r, (0, top_h + PAD))

out_path = os.path.join(_fig_dir, "fig_native_networks_composite.png")
canvas.save(out_path, dpi=(200, 200))
print(f"wrote {out_path}  size={canvas.size}")

# Also a PDF version for LaTeX inclusion.
out_pdf = os.path.join(_fig_dir, "fig_native_networks_composite.pdf")
canvas.save(out_pdf, "PDF", resolution=200)
print(f"wrote {out_pdf}")
