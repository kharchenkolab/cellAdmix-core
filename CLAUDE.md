# Project guidance

## Commits

- Write short commit messages: a single imperative subject line, no body
  unless truly essential.
- Do not add Co-Authored-By lines or any other signatures/trailers.

## Documentation

- Keep the documentation coherently current: maintain the user-facing
  narrative as an integral whole that reflects the current state of the code.
  When behavior, defaults, or APIs change, rework the affected pages
  (README.md, docs/*.md, roxygen and docstrings) as part of the same change so
  each page still reads as one unified description.
- Never keep a journal of updates in the docs: no "new in this version",
  "previously/now", or appended change notes - integrate changes into the
  existing narrative as if the code had always worked this way.
