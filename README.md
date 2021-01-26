# Stepik exercises - Coursera Bioinformatics II - UC Sandiego

Python codes and tests for interactive excercises on Stepik
The course introduces algorithms to:
1. Assemble a genome from reads and contigs
2. Reconstruct antibiotic sequences

To run all the tests

```
pytest
```

## Week1
In week1, we learn how to construct a Debruijn graph from segments of DNA (reads)
I use dictionary to store a Debruijn graph. The keys and values of the dictionary represent the nodes of the graph. The mapping of keys and values are edges that connect the nodes.

## Week2
In week2, we learn how to construct Eulerian cycles from a Debruijn graph.
An Eulerian cycle show a way to cross all edges of a Debruijn graph exactly one time
We then build longer reads from short pairs of reads, and construct a Debruijn graph from the read pairs. This makes the graph less tangle and reduces the edges.

## Week3
In week3, we learn how to construct a peptide from a spectral data obtained from an experiment
First we generate theoretical spectra from single amino acid
Then we compare each theoretical spectrum with the experimental spectrum and keep the peptide whose spectrum completely match the spectrum from the experiment

## Week4
In week4, we learn how to construct a peptide from a spectral data obtained from an experiment, but allowing for mismatches as results from experiment are often messy.