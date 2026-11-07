"""The figure of the neighbor-association test in docs/generative.md: the
regression that follows the model's fit, with the components colored as in the
model's formula and its named quantities pointed out.

Needs pdflatex (with the standalone, sansmath and tikz packages) and
pdftoppm. Writes docs/figures/generative_fig4.png.
"""
import subprocess
import tempfile
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "docs" / "figures" / "generative_fig4.png"
WIDTH = 1500


def trim(img, pad):
    """Crop an image's near-white margins, keeping `pad` pixels around the content."""
    box = Image.eval(img.convert("L"), lambda v: 255 if v < 235 else 0).getbbox()
    return img.crop((max(box[0] - pad, 0), max(box[1] - pad, 0),
                     min(box[2] + pad, img.width), min(box[3] + pad, img.height)))


with tempfile.TemporaryDirectory() as tmp:
    for _ in range(2):  # the second pass places the overlay labels
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "-output-directory", tmp,
                        str(HERE / "test_formula.tex")], check=True, capture_output=True)
    subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile",
                    str(Path(tmp) / "test_formula.pdf"), str(Path(tmp) / "formula")], check=True)
    formula = trim(Image.open(Path(tmp) / "formula.png").convert("RGB"), 10)
fig = formula.resize((WIDTH, round(formula.height * WIDTH / formula.width)), Image.LANCZOS)
fig.quantize(colors=256, method=Image.Quantize.MEDIANCUT).save(OUT, optimize=True)
print(f"wrote {OUT} ({fig.width} x {fig.height})")
