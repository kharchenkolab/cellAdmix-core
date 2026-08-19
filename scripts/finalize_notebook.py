#!/usr/bin/env python3
"""Post-process a quarto-rendered ipynb into its publishable form.

Quarto's md->ipynb writer keeps only the first code segment of any chunk
whose outputs interleave with its code (all quarto versions tested, 1.0-1.6),
and it stamps Python notebook metadata. This script makes the rendered
notebook whole:

1. restores every code cell's full source from its Rmd chunk,
2. restores every markdown cell's text from the Rmd prose (undoing quarto's
   smart-quote styling, so notebook text equals the Rmd source verbatim),
3. sets the R kernelspec / language_info metadata,
4. sets retina display widths on figure outputs (png pixel width / 2),
5. assigns stable position-based cell ids (quarto randomizes them per
   render, which churns git diffs).

Usage: finalize_notebook.py path/to/notebook   (no extension; expects
       notebook.Rmd and notebook.ipynb side by side)

Fails loudly on any chunk/cell count mismatch.
"""
import base64
import json
import re
import struct
import sys

R_METADATA = {
    "kernelspec": {"display_name": "R", "language": "R", "name": "ir"},
    "language_info": {
        "codemirror_mode": "text",
        "file_extension": ".r",
        "mimetype": "text/plain",
        "name": "text",
        "pygments_lexer": "text",
        "version": "4.4.3",
    },
}


def parse_rmd(rmd_path):
    """Return (title, visible R chunk bodies, non-R fenced bodies, prose blocks)."""
    lines = open(rmd_path).read().splitlines()
    title = None
    i = 0
    if lines and lines[0].strip() == '---':
        i = 1
        while i < len(lines) and lines[i].strip() != '---':
            m = re.match(r'^title:\s*["\']?(.*?)["\']?\s*$', lines[i])
            if m:
                title = m.group(1)
            i += 1
        i += 1
    r_re = re.compile(r'^```\{r([^}]*)\}\s*$')
    other_re = re.compile(r'^```(\w+)\s*$|^```\{(?!r[,}\s])(\w+)[^}]*\}\s*$')
    chunks, fences, prose, cur = [], [], [], []
    while i < len(lines):
        m = r_re.match(lines[i])
        o = other_re.match(lines[i]) if not m else None
        if m or o:
            if any(l.strip() for l in cur):
                prose.append('\n'.join(cur).strip('\n'))
            cur = []
            body = []
            i += 1
            while i < len(lines) and lines[i].strip() != '```':
                body.append(lines[i])
                i += 1
            if m:
                header = m.group(1).replace(' ', '')
                if 'include=FALSE' not in header and 'echo=FALSE' not in header:
                    chunks.append('\n'.join(body))
            else:
                fences.append('\n'.join(body))
        else:
            cur.append(lines[i])
        i += 1
    if any(l.strip() for l in cur):
        prose.append('\n'.join(cur).strip('\n'))
    return title, chunks, fences, prose


def as_source(text):
    lines = [l + '\n' for l in text.rstrip('\n').splitlines()]
    if lines:
        lines[-1] = lines[-1].rstrip('\n')
    return lines


def finalize(base):
    title, chunks, fences, prose = parse_rmd(base + '.Rmd')
    fence_firsts = {f.splitlines()[0] for f in fences if f.splitlines()}
    nb = json.load(open(base + '.ipynb'))

    md_cells = [c for c in nb['cells'] if c['cell_type'] == 'markdown']
    if len(md_cells) == len(prose):
        targets = ["# {}\n\n{}".format(title, prose[0])] + prose[1:]
    elif len(md_cells) == len(prose) + 1:
        targets = ["# {}".format(title)] + prose
    else:
        sys.exit(f"ERROR {base}: {len(md_cells)} markdown cells "
                 f"vs {len(prose)} prose blocks")
    md_resynced = 0
    for cell, text in zip(md_cells, targets):
        if ''.join(cell['source']).rstrip('\n') != text:
            cell['source'] = as_source(text)
            md_resynced += 1

    ci = resynced = figures = 0
    for cell in nb['cells']:
        if cell['cell_type'] != 'code':
            continue
        src = ''.join(cell['source'])
        first = src.splitlines()[0] if src.splitlines() else ''
        if first in fence_firsts:
            continue
        if ci >= len(chunks):
            sys.exit(f"ERROR {base}: more code cells than visible chunks")
        if src.rstrip('\n') != chunks[ci].rstrip('\n'):
            cell['source'] = as_source(chunks[ci])
            resynced += 1
        for out in cell.get('outputs', []):
            png = out.get('data', {}).get('image/png')
            if png:
                width = struct.unpack('>I', base64.b64decode(png)[16:20])[0]
                out.setdefault('metadata', {})['image/png'] = {'width': width // 2}
                figures += 1
        ci += 1
    if ci != len(chunks):
        sys.exit(f"ERROR {base}: only {ci}/{len(chunks)} chunks placed")

    for idx, cell in enumerate(nb['cells']):
        cell['id'] = 'cell-%03d' % idx

    nb['metadata'] = R_METADATA
    json.dump(nb, open(base + '.ipynb', 'w'), indent=1, ensure_ascii=False)
    print(f"{base}: {resynced} code + {md_resynced} markdown cells resynced, "
          f"{figures} figure widths set, {ci} chunks verified")


if __name__ == '__main__':
    for arg in sys.argv[1:]:
        for suffix in ('.Rmd', '.ipynb'):
            if arg.endswith(suffix):
                arg = arg[:-len(suffix)]
        finalize(arg)
