# genbankx

A lightweight GenBank parser with ordered CDS IntervalList and translation capabilities.

## Features

- Parse standard GenBank files
- Extract CDS features with ordered intervals
- Handle join() and complement() location syntax
- Translate DNA sequences to protein sequences
- Export sequences in FASTA format

## Installation
bash
pip install genbankx

## Usage
### Python API
python
from genbankx import parse_genbank, translate

Parse a GenBank file
with open("example.gb", "r") as f:
records = list(parse_genbank(f))
Extract and translate CDS sequences
for record in records:
for cds in record.cds:
dna_seq = "".join(record.sequence[iv.start-1:iv.end] for iv in cds)
if cds.strand.value == -1:
dna_seq = reverse_complement(dna_seq)
protein_seq = translate(dna_seq)
print(f"{cds.gene}: {protein_seq}")

### Command Line Interface
bash
Extract CDS DNA sequences
genbankx input.gb -o cds.fasta
Extract protein sequences
genbankx input.gb --proteins -o proteins.fasta

## License
MIT

