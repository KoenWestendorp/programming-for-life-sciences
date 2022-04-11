
## Usage 

### Package

For use of the package, see the [docs](file:./doc/aligner.html).

### CLI

```bash
$ python -m aligner 
```

## Installation

### Pip install

When in the aligner package directory, run

```bash
$ pip install .
```

For development, it is useful to run:

```bash
$ pip install -e .
```

The `-e` flag makes the installation editable, which means that the source
files can be edited and on import after saving that edit, the package will be
updated.

### Generate docs

To generate docs, run:

```bash
$ pdoc ./aligner -d numpy -o doc 
```

The docs can then be found in the _doc_ directory, such as
[`./doc/aligner.html`](file:./doc/aligner.html).

### Build the report

The source file for the report is `report.md`.

It can be converted to a number of different output formats using
[pandoc](https://pandoc.org/).

```bash
$ pandoc --pdf-engine=xelatex -i report.md -o report.pdf
$ pandoc --mathml --standalone --self-contained -i report.md -o report.html
```

## Dependencies

### Build dependecies

- [matplotlib](https://matplotlib.org/)
- [numpy](https://numpy.org/)

### Notebook dependencies

- [MDAnalysis](https://www.mdanalysis.org/)
- [nglviewer](https://nglviewer.org/)

## Testing

Tests are located in the _test_ directory.

Tests can be run using PyTest:

```bash
$ pytest
```
