# Mendel's First Law Solver

Solution to the [Rosalind IPRB](https://rosalind.info/problems/iprb/) problem.

Given a population of `k` homozygous dominant (AA), `m` heterozygous (Aa),
and `n` homozygous recessive (aa) organisms, calculates the probability that
two randomly selected mating organisms produce offspring with the dominant phenotype.

## Usage

```bash
# Interactive mode
python mendels_first_law.py

# Direct input
python mendels_first_law.py 2 2 2
# Output: 0.78333
```

## Features

- Class-based design (`Population`, `MendelSolver`, `InputHandler`)
- Exact rational arithmetic using Python's `Fraction` class
- Full input validation with helpful error messages
- Detailed cross-by-cross breakdown output
- Works as a standalone script or importable module
