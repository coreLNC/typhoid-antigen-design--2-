"""Render the ranked epitope segments on the three structures (Figure 5) with 3Dmol.js.

usage (from the repository root):
    npm install 3dmol                      # provides node_modules/3dmol/build/3Dmol-min.js
    pip install playwright pillow matplotlib && playwright install chromium
    python3 -I scripts/render_structures.py node_modules/3dmol/build/3Dmol-min.js

Ribbon colour: ranked segments in the protein colour, HlyE 175-188 (flagged literature segment)
in pink, everything else grey; the two Typhi-selective residues (HlyE N213, CdtB G228) are red spheres.
HlyE = 1QOY (E. coli HlyE, numbering identical to Typhi at these positions); CdtB and YncE = AlphaFold models.
Writes figures/fig5_structures.png.
"""
import csv, json, os, shutil, sys, tempfile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageChops
from playwright.sync_api import sync_playwright

JS = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else sys.exit(__doc__)
COL = {"hlyE": "#2a78d6", "cdtB": "#eb6834", "yncE": "#1baf7a"}
FILES = {"hlyE": "structures/1QOY.pdb", "cdtB": "structures/cdtB_mature_23-269.pdb", "yncE": "structures/yncE_mature_31-353.pdb"}
TYPHI_SEL = {"hlyE": 213, "cdtB": 228}
segs = {p: [] for p in COL}
for r in csv.DictReader(open("results/epitopes/epitope_ranking_v3.csv")):
    segs[r["protein"]].append((int(r["start"]), int(r["end"])))

def html(p, rot):
    js = f'let v=$3Dmol.createViewer("c",{{backgroundColor:"white"}});v.addModel({json.dumps(open(FILES[p]).read())},"pdb");'
    js += 'v.setStyle({},{cartoon:{color:"#d9d8d3"}});'
    for s, e in segs[p]:
        js += f'v.setStyle({{resi:"{s}-{e}"}},{{cartoon:{{color:"{COL[p]}"}}}});'
    if p == "hlyE":
        js += 'v.setStyle({resi:"175-188"},{cartoon:{color:"#c2185b"}});'
    if p in TYPHI_SEL:
        js += f'v.addStyle({{resi:{TYPHI_SEL[p]}}},{{stick:{{color:"#e34948",radius:0.35}},sphere:{{color:"#e34948",radius:0.9}}}});'
    js += f'v.zoomTo();v.rotate({rot},"y");v.render();window.done=1;'
    return f'<html><body style="margin:0"><div id="c" style="width:420px;height:560px;position:relative"></div><script src="3Dmol-min.js"></script><script>{js}</script></body></html>'

tmp = tempfile.mkdtemp(); shutil.copy(JS, os.path.join(tmp, "3Dmol-min.js")); shots = {}
with sync_playwright() as pw:
    b = pw.chromium.launch(args=["--no-sandbox", "--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"])
    for p in COL:
        for rot in (0, 90):
            page_html = os.path.join(tmp, "t.html"); open(page_html, "w").write(html(p, rot))
            pg = b.new_page(viewport={"width": 420, "height": 560}, device_scale_factor=2)
            pg.goto("file://" + page_html); pg.wait_for_function("window.done==1", timeout=30000); pg.wait_for_timeout(500)
            out = os.path.join(tmp, f"{p}_{rot}.png"); pg.screenshot(path=out); pg.close(); shots[(p, rot)] = out
    b.close()

def crop(path):
    im = Image.open(path).convert("RGB"); W, H = im.size
    bb = ImageChops.difference(im, Image.new("RGB", im.size, "white")).getbbox()
    return im.crop((max(bb[0] - 15, 0), max(bb[1] - 15, 0), min(bb[2] + 15, W), min(bb[3] + 15, H)))

titles = {"hlyE": "HlyE (1QOY, 1–303)", "cdtB": "CdtB (AlphaFold, 23–269)", "yncE": "YncE (AlphaFold, 31–353)"}
plt.rcParams.update({"font.family": "serif", "font.size": 9})
fig, axs = plt.subplots(2, 3, figsize=(6.4, 5.6), gridspec_kw={"hspace": .02, "wspace": .02})
for j, p in enumerate(COL):
    for i, rot in enumerate((0, 90)):
        ax = axs[i, j]; ax.imshow(crop(shots[(p, rot)])); ax.axis("off")
        if i == 0: ax.set_title(titles[p], fontsize=9)
        ax.text(.02, .02, "front" if rot == 0 else "rotated 90°", transform=ax.transAxes, fontsize=7, color="#52514e")
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/fig5_structures.png", dpi=220, bbox_inches="tight"); print("wrote figures/fig5_structures.png")
