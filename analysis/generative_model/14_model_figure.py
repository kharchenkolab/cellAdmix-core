"""The model figure of docs/generative.md: the overlapping-cells
illustration above the model's formula, its components colored as in the
illustration and its named quantities pointed out.

Needs pdflatex (with the standalone, sansmath and tikz packages) and
pdftoppm. Writes docs/figures/generative_model.png.
"""
import subprocess
import tempfile
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "docs" / "figures" / "generative_model.png"
WIDTH = 1400  # pixels of the composed figure
CELLS_WIDTH = WIDTH // 3  # the illustration, centered above the formula

with tempfile.TemporaryDirectory() as tmp:
    for _ in range(2):  # the second pass places the overlay labels
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "-output-directory", tmp,
                        str(HERE / "model_formula.tex")], check=True, capture_output=True)
    subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile",
                    str(Path(tmp) / "model_formula.pdf"), str(Path(tmp) / "formula")], check=True)
    formula = Image.open(Path(tmp) / "formula.png").convert("RGB")

cells = Image.open(HERE / "figures" / "overlapping_cells.png").convert("RGB")
# trim the illustration's blank margins (its background is near-white)
box = Image.eval(cells.convert("L"), lambda v: 255 if v < 235 else 0).getbbox()
cells = cells.crop((0, max(box[1] - 30, 0), cells.width, min(box[3] + 30, cells.height)))
cells = cells.resize((CELLS_WIDTH, round(cells.height * CELLS_WIDTH / cells.width)), Image.LANCZOS)
fw = round(WIDTH * 0.92)
formula = formula.resize((fw, round(formula.height * fw / formula.width)), Image.LANCZOS)

gap = 20
fig = Image.new("RGB", (WIDTH, cells.height + gap + formula.height + 30), "white")
fig.paste(cells, ((WIDTH - CELLS_WIDTH) // 2, 0))
fig.paste(formula, ((WIDTH - fw) // 2, cells.height + gap))
fig.quantize(colors=256, method=Image.Quantize.MEDIANCUT).save(OUT, optimize=True)
print(f"wrote {OUT} ({fig.width} x {fig.height})")
