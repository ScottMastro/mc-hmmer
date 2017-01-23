#!/usr/bin/env python3

import random
import math
from moves import mcmc_move
from data import *
from viterbi import *
import operator
import pickle
import csv
import sys

class Model:
    """A class representing an HMM.

    e is the matrix of emission probabilities
    a is the matrix of transition probabilities
    labels is a list of the labels for the nodes (ie. alpha, beta, coil, none)
    """
    def __init__(self, a, e, labels):
        self.a = a
        self.e = e
        self.labels = labels

def consensus(seqs):
    """Return the plurality consensus of a list of sequences"""
    cons = []
    for i in range(len(seqs[0])):
        counts = {}
        for s in seqs:
            try:
                counts[s[i]] += 1
            except KeyError:
                counts[s[i]] = 1
        cons.append(sorted(counts.items(), key=lambda x: -x[1])[0][0])
    return cons

def log_likelihood(m, seqs):
    loglik = 0
    for s in seqs:
        loglik += forwardAlgorithm(m.e, m.a, s)[2]
    
    return loglik

def get_path(m, seq):
    path = run_viterbi(m.e, m.a, seq, True)
    return [m.labels[x-1] for x in path]

def main():

    # how many iterations to do total, and how often to sample models
    n_iter = 10000
    sample_every = 100

    # don't start sampling until after burnin steps
    burnin = 1000

    # load the data
    print("Loading training data")
    data, states = load_train_data()

    # harvest amino acid frequencies
    print("Measuring amino acid frequencies")
    e = harvest_e(data, states)

    # initialize the model from the assignment
    a = [[1.0/3, 1.0/3, 1.0/3],
         [0.9, 0, 0.1],
         [0, 0.89, 0.11],
         [0.1, 0.11, 0.79]]
    labels = ["A", "B", "O"]
    m = Model(a, e, labels)

    log_likelihood_m = log_likelihood(m, data)

    print("Loading test data")
    test_data, test_states = load_test_data()
    counts = []

    # to log the results
    header = ["loglik", "move", "lik.ratio", "accept"]
    writer = csv.DictWriter(sys.stderr, fieldnames=header, delimiter="\t")
    writer.writeheader()

    # run MCMC
    sampled_models = []
    for i in range(n_iter):
                
        print("------------------------------")
        print("Iteration", i, "-", len(m.labels), "states")

        row = dict.fromkeys(header)

        # generate a new model
        row["move"], m2 = mcmc_move(m)

        # find its log likelihood 
        log_likelihood_m2 = log_likelihood(m2, data)

        # accept the new model with probability equal to the likelihood ratio
        ratio = math.exp(log_likelihood_m2 - log_likelihood_m)
        print("Likelyhood ratio:", ratio)
        row["lik.ratio"] = ratio
        row["accept"] = "FALSE"
        if (len(m2.labels) <= 10 and
            (log_likelihood_m2 > log_likelihood_m or
             random.random() < ratio)):
            m = m2
            log_likelihood_m = log_likelihood_m2
            print("Switch to new model")
            row["accept"] = "TRUE"

        row["loglik"] = log_likelihood_m
        writer.writerow(row)

        # keep some subset of the models
        if i % sample_every == 0 and i > burnin:
            print("Sampled iteration", i)
            sampled_models.append(m)

    with open("models.pkl", "wb") as f:
        pickle.dump(sampled_models, f)

    data, true_states = load_test_data()

    results = [0, 0] # number of correct/incorrect bases
    for i in range(len(data)):
        states = consensus([get_path(m, data[i]) for m in sampled_models])
        for j in range(len(true_states[i])):
            if states[j] == true_states[i][j]:
                results[0] += 1
            else:
                results[1] += 1

    print ("{} correct bases, {} incorrect bases".format(*results))

    return 0

if __name__ == "__main__":
    main()
