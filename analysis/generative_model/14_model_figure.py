"""The model figure of docs/generative.md: the overlapping-cells
illustration beside the model's formula, its components colored as in the
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
FORMULA_WIDTH = 1100  # pixels of the formula; the illustration sits to its left
CELLS_SCALE = 1.25  # the illustration's height relative to the formula's
GAP = 40


def trim(img, pad):
    """Crop an image's near-white margins, keeping `pad` pixels around the content."""
    box = Image.eval(img.convert("L"), lambda v: 255 if v < 235 else 0).getbbox()
    return img.crop((max(box[0] - pad, 0), max(box[1] - pad, 0),
                     min(box[2] + pad, img.width), min(box[3] + pad, img.height)))


with tempfile.TemporaryDirectory() as tmp:
    for _ in range(2):  # the second pass places the overlay labels
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "-output-directory", tmp,
                        str(HERE / "model_formula.tex")], check=True, capture_output=True)
    subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile",
                    str(Path(tmp) / "model_formula.pdf"), str(Path(tmp) / "formula")], check=True)
    formula = trim(Image.open(Path(tmp) / "formula.png").convert("RGB"), 10)
formula = formula.resize((FORMULA_WIDTH, round(formula.height * FORMULA_WIDTH / formula.width)),
                         Image.LANCZOS)

cells = trim(Image.open(HERE / "figures" / "overlapping_cells.png").convert("RGB"), 10)
ch = round(formula.height * CELLS_SCALE)
cells = cells.resize((round(cells.width * ch / cells.height), ch), Image.LANCZOS)

height = max(cells.height, formula.height) + 20
fig = Image.new("RGB", (cells.width + GAP + formula.width + 20, height), "white")
fig.paste(cells, (10, (height - cells.height) // 2))
fig.paste(formula, (10 + cells.width + GAP, (height - formula.height) // 2))
fig.quantize(colors=256, method=Image.Quantize.MEDIANCUT).save(OUT, optimize=True)
print(f"wrote {OUT} ({fig.width} x {fig.height})")
