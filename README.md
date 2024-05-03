

## Repository Summary

Python implementation of non-homologous stellar collapse simulation, based on IDL implementation by S. Cranmer.


## Installation
This package will run with the standard python library.

## Usage

The model can be run using the following syntax:\
`python non_homologous.py --epsilon [EPSILON] --outdir [OUTDIR]`

The arguments are:\
`--epsilon`: (float) .\
`--outdir`: (string) .

## Examples
To produce the plots in the `plots` directory, see below.

$\epsilon = 0:$
`python non_homologous.py --epsilon 0.0 --outdir ./plots`

$\epsilon = 0.1:$
`python non_homologous.py --epsilon 0.1 --outdir ./plots`

$\epsilon = 0.9:$
`python non_homologous.py --epsilon 0.9 --outdir ./plots`

