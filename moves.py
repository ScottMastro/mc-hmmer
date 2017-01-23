#!/usr/bin/env python3

import random
import copy
import math
import sys

def check_model(m):
    """Make sure a model is coherent"""
    assert(all([math.fsum(row) - 1 < 1e-10 for row in m.a]))
    for j in range(len(m.e[0])):
        assert(math.fsum([m.e[i][j] for i in range(len(m.e))]) - 1 < 1e-10)
    assert(len(m.a) == len(m.labels)+1)
    assert(all([len(x) == len(m.labels) for x in m.a]))
    assert(len(m.e) == 20)
    assert(all([len(x) == len(m.labels) for x in m.e]))

def divsum(v):
    """Divide a vector by its sum"""
    sum = math.fsum(v)
    return [x/sum for x in v]

def is_strongly_connected(a):
    """Check if an adjacency matrix represents a strongly connected digraph"""
    # first DFS
    stack = [0]
    found = set([])
    while stack != []:
        v = stack.pop()
        if not v in found:
            found.add(v)
            children = [w for w in range(len(a)) if a[v][w] != 0 and w != v]
            for w in children:
                stack.append(w)
    if len(found) != len(a):
        return False

    # second DFS, with reversed edges
    stack = [0]
    found = set([])
    while stack != []:
        v = stack.pop()
        if not v in found:
            found.add(v)
            children = [w for w in range(len(a)) if a[w][v] != 0 and w != v]
            for w in children:
                stack.append(w)
    return len(found) == len(a)

def split(m):
    """Split one node in m into two nodes, preserving edges"""

    # make a new model
    m = copy.deepcopy(m)
    nnodes = len(m.labels)

    # choose a node to split
    node = random.randint(0, nnodes-1)

    # add the node to the labels
    m.labels.append(m.labels[node])

    # split edges into the new node
    m.a[0] = [1.0/(nnodes+1) for i in range(nnodes+1)]
    for i in range(nnodes):
        m.a[i+1][node] /= 2
        m.a[i+1].append(m.a[i+1][node])

    # add edges out of the new node
    m.a.append(m.a[node+1][:])

    # give the new node the same emission probabilities
    for i in range(len(m.e)):
        m.e[i].append(m.e[i][node])

    # make sure everything worked ok
    #check_model(m)

    return m

def join(m):
    """Join two nodes in m into one"""

    # if there are no nodes to join, just propose the same model
    if len(m.labels) == 3:
        return m
    
    # make a new model
    m = copy.deepcopy(m)
    nnodes = len(m.labels)

    # choose nodes to join
    node1_choices = [i for i in range(nnodes) if m.labels.count(m.labels[i]) > 1]
    node1 = random.choice(node1_choices)
    node2_choices = [i for i in range(nnodes) if m.labels[i] == m.labels[node1] and i != node1]
    node2 = random.choice(node2_choices)
    if (node1 > node2):
        node1, node2 = node2, node1

    # combine edges going into and out of each node
    m.a[0] = [1.0/(nnodes-1) for i in range(nnodes-1)]
    for i in range(nnodes):
        m.a[i+1][node1] += m.a[i+1][node2]
        if i not in [node1, node2]:
            m.a[node1+1][i] = (m.a[node2+1][i] + m.a[node1+1][i])/2
    m.a[node1+1][node1] = (m.a[node1+1][node1] + m.a[node2+1][node1])/2
    for i in range(nnodes):
        m.a[i+1].pop(node2)
    m.a.pop(node2+1)

    # average emission probabilities
    for i in range(20):
        m.e[i][node1] = (m.e[i][node1] + m.e[i][node2])/2
    for i in range(20):
        m.e[i].pop(node2)

    # remove the extra label
    m.labels.pop(node2)

    # make sure everything worked ok
    #check_model(m)

    return m

def add_edge(m):
    """Add an adjacency in m"""

    # find non-existant edges
    nnodes = len(m.labels)
    choices = []
    for i in range(nnodes):
        for j in range(nnodes):
            if (m.a[i][j] == 0):
                choices.append((i, j))

    # if there are no non-existant edges, propose the same model
    if (len(choices) == 0):
        return m

    # make a new model
    m = copy.deepcopy(m)

    # choose the edge to add, and its transition probability
    new_i, new_j = random.choice(choices)
    p = random.random()

    # update the transition matrix
    m.a[new_i][new_j] = p
    m.a[new_i] = divsum(m.a[new_i])

    # make sure it worked
    #check_model(m)

    return m

def delete_edge(m):
    """Remove an adjacency in m"""

    # make a copy of the graph
    m2 = copy.deepcopy(m)
    
    # find all the edges
    nnodes = len(m.labels)
    choices = []
    for i in range(nnodes):
        for j in range(nnodes):
            if (m.a[i+1][j] != 0):
                choices.append((i, j))

    # choose an edge to remove
    del_i, del_j = random.choice(choices)

    # update the transition matrix
    m2.a[del_i+1][del_j] = 0
    try:
        m2.a[del_i+1] = divsum(m2.a[del_i+1])
    except ZeroDivisionError: # we deleted the only edge out of a state
        return (m)

    # if the model isn't strongly connected anymore, re-propose the old model
    if not is_strongly_connected(m2.a[1:]):
        return (m)

    # make sure it went ok
    #check_model(m2)

    return m2

def edit_transition(m):
    """Alter a transition probability in m"""

    # make a copy of the graph
    m2 = copy.deepcopy(m)
    
    # find all the edges
    nnodes = len(m.labels)
    choices = []
    for i in range(nnodes):
        for j in range(nnodes):
            if (m.a[i][j] != 0):
                choices.append((i, j))

    # choose an edge to edit and its new probability
    i, j = random.choice(choices)
    p = random.random()

    # update the transition matrix
    m2.a[i][j] = p
    m2.a[i] = divsum(m2.a[i])

    # make sure it went ok
    #check_model(m2)

    return m2

def mcmc_move(m):
    """Perform a random MCMC move on m, and return a new model"""
    move_names = ["number.of.state.increase", "number.of.state.decrease",
                  "transition.prob.change", "add.edge", "remove.edge"]

    r = random.random()
    moves = [join, split, add_edge, delete_edge, edit_transition]
    for i, move in enumerate(moves, start=1):
        if r < i/len(moves):
            return move_names[i-1], move(m)
