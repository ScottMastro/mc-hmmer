#!/usr/bin/env python3
import math

# because python3 doesn't like xrange
try:
    xrange
except NameError:
    xrange = range
    
def run_viterbi(e, t, s, return_path) :
    viterbiDecoding(e, t, s)
    out = forwardAlgorithm(e, t, s)
    B = backwardAlgorithm(e, t, s, out[1])
    
    return posteriorPath(e, t, s, out[0], B, out[1], return_path)
    
class cell:
    #class to store matrix cell information for traceback
    value = -float("inf")
    predecessor = -1
    
def viterbiDecoding(e, t, s) :
    #Step 0: Initialize
    n = len(s)
    states = len(e[1])

    V = [[cell() for x in xrange(states+1)] for x in xrange(n+1)] 
    V[0][0].value = 0
    
    #Step 1: Fill matrix
    for i in xrange(1, n+1) :
        for j in xrange(1, states+1) :
            index = chars[s[i-1:i]]
            ej = math.log(e[index][j-1])
            maxk = -float("inf")
            pred = 0
            for k in xrange(states+1) :
                Vk = V[i-1][k].value
                akj = -float("inf")
                if t[k][j-1] != 0:
                    akj = math.log(t[k][j-1])
                else:
                    continue
                    
                if maxk < Vk + akj :
                    maxk = Vk + akj
                    pred = k
                    
            V[i][j].value = ej + maxk
            V[i][j].predecessor = pred
    
    score = -float("inf")
    lastState = -1
    for j in xrange(1,states+1) :
        if V[n][j].value > score :
            score = V[n][j].value
            lastState = j
            
    concat = []
    while lastState != -1 and lastState != 0 :
        concat.append(lastState)
        lastState = V[i][lastState].predecessor
        i = i - 1 
    return(concat)    
    
def forwardAlgorithm(e, t, s) :
    #Step 0: Initialize
    n = len(s)
    states = len(e[1])
 
    F = [[0 for x in xrange(states+1)] for x in xrange(n+1)] 
    F[0][0] = 1
    scale = [1 for x in xrange(n+1)]
    sum_logscale = 0

    #Step 1: Fill matrix
    for i in xrange(1, n+1) :
        sum = 0
        for j in xrange(1, states+1) :
            index = chars[s[i-1:i]]
            ej = e[index][j-1]
            
            sum_fkakl = 0
            for k in xrange(states+1) :
                sum_fkakl += F[i-1][k] * t[k][j-1]
                
            F[i][j] = ej * sum_fkakl
            sum += F[i][j]
        
        scale[i] = sum
        sum_logscale += math.log(sum)

        for l in xrange(states+1) :
            F[i][l] = F[i][l]/sum
   
    logprob = 0;
    for j in xrange(1, states+1) :
        logprob += F[n][j]
        
    logprob = -(math.log(logprob) - sum_logscale)

    return [F, scale, logprob]
    
def backwardAlgorithm(e, t, s, scale) :
    #Step 0: Initialize
    n = len(s)
    states = len(e[1])

    B = [[0 for x in xrange(states+1)] for x in xrange(n+1)] 
    
    for j in xrange(1, states+1) :
        B[n][j] = 1;
        
    #Step 1: Fill matrix
    for i in xrange(n-1, 0, -1) :
        for j in xrange(1, states+1) :
            index = chars[s[i:i+1]]
            sumk = 0
            for k in xrange(1, states+1) :
                sumk += B[i+1][k] * t[j][k-1] * e[index][k-1]
            
            B[i][j] = sumk / scale[i]
            
    return B
    
def posteriorPath(e, t, s, F, B, scale, return_path) :
    n = len(s)
    states = len(e[1])
    
    if return_path :
        concat = []
        for i in xrange(1, n+1) :    
            bestState = -1
            mostProb = -1

            for j in xrange(1, states +1) :    
                p = F[i][j] * B[i][j]

                if p > mostProb :
                    mostProb = p
                    bestState = j
            
            concat.append(bestState)
        
        #print(s)
        #print(concat)
        return concat
    
    #todo: return something?
    return 0

 #ACDEFGHIKLMNPQRSTVWYX
chars = {"A":0, "C":1, "D":2, "E":3, "F":4, "G":5, "H":6, "I":7, "K":8, 
"L":9, "M":10, "N":11, "P":12, "Q":13, "R":14, "S":15, "T":16, 
"V":17, "W":18,"Y":19}

def setupTest() :

    e = [
    [0.106134,0.059686,0.062796],
    [0.013200,0.019820,0.015719],
    [0.049837,0.043423,0.072972],
    [0.091491,0.050432,0.058141],
    [0.041080,0.052439,0.032509],
    [0.034454,0.067148,0.108530],
    [0.021396,0.023232,0.024904],
    [0.060069,0.085248,0.037289],
    [0.068427,0.050830,0.060637],
    [0.121081,0.093883,0.068519],
    [0.026845,0.020555,0.022728],
    [0.032872,0.032855,0.057896],
    [0.023658,0.027213,0.076311],
    [0.048216,0.029758,0.034650],
    [0.059596,0.045907,0.046746],
    [0.049883,0.055435,0.075226],
    [0.041552,0.064133,0.059066],
    [0.062080,0.117532,0.048200],
    [0.014381,0.015942,0.010063],
    [0.033735,0.044513,0.027079]]
    
    t = [[1.0/3.0, 1.0/3.0, 1.0/3.0],
    [0.9, 0.0, 0.1],
    [0.0, 0.89, 0.11],
    [0.1, 0.11, 0.79]]
    
    s = "MKSIGVVRKVDELGRIVMPIELRRALDIAIKDSIEFFVDGDKIILKKYKPHGVCLMTGEITSENKEYGNGKITLSPEGAQLLLEEIQAALKE"
    run_viterbi(e,t,s, True)
    
if __name__ == "__main__":
    setupTest()
