# Profile Hidden Markov Models

I implemented a Profile Hidden Markov Model (PHMM) in Python, which is used for classification of protein or gene sequences to families.

[![Build Status](https://travis-ci.org/ikar1234/PHMM.svg?branch=master)](https://travis-ci.org/ikar1234/PHMM)

## Input

The model takes a either a multiple sequence alignment as an input or a set of unaligned sequences.


## Training
- Known parameters (multiple alignment was given)
 The model uses the Viterbi algorithm for training.

- Unknown parameters
 When the parameters are unknown, the model first estimates them using the test data with the Baum-Welch algorithm. Then, it aligns the sequences to a multiple alignment.
## Evaluation

The model can estimate the observation probability of a sequence given a multiple alignment.
