# mc-hmmer

## *Protein structure prediction with MCMC/HMM*


---



**Group: Rosemary McCloskey, James Topham, Scott Mastromatteo**

**Title: Protein Structure Prediction Using a Markov Chain Monte Carlo Approach**


Computational protein structure prediction is the subject of considerable research, with one of the most widely used approaches being hidden Markov models (HMMs). In its most basic form, an HMM for protein structure prediction would have one hidden state for each structural element (say, a-helices, ß-sheets, coils, and no structural element). More sophisticated models employ more than one hidden state per element, each with distinct emission probabilities. The structure and parameters of these models may be either biologically informed, automatically learned from a training data set, or a combination of both.


A drawback of many of these HMM-based methods is that they do not account for model uncertainty. The structure and parameters of a single HMM are learned, through maximum likelihood estimation or heuristics, and then that HMM is taken to be the “true” model when making predictions about new proteins. Even if the learned model has the globally optimal posterior probability, here may be other models which have nearly equal predictive power, but these are ignored in the output of the learning process.


Markov chain Monte Carlo (MCMC) sampling offers a way to overcome this limitation. MCMC (in particular, the Metropolis-Hastings algorithm) is an iterative process that begins at a random model. At each step, a new model is proposed by making some small change to the current model. If the proposed new model has a higher posterior probability, it is accepted with certainty; otherwise, it is accepted with probability proportional to the posterior probability ratio of the two models. It can be shown that this process will eventually revisit models at a frequency approximating the posterior probability distribution. By using the consensus of a sample of models visited by this procedure, we can obtain a protein structure prediction informed by the entire posterior probability distribution over all models, which may be more accurate than a single estimate. Moreover, because MCMC is a Bayesian method, it affords the opportunity to specify biologically-informed priors on model parameters.


We propose to investigate extensions of a three-state HMM for protein structure prediction via MCMC sampling. In particular, we will implement the Metropolis-Hastings sampler, using a proposal algorithm set forth by Ranganathan and Dellaert [1]. We will sample topologies only; the model parameters will be estimated by either Viterbi training or the Baum-Welch method (to be decided). 


We will train and test the model on a curated dataset from the RCSB Protein Data Bank. These are protein sequences which have been annotated with structural elements. Here is an example of input obtained from RCBS:


```
DIVLTQSPAIMSASLGSSVTLTCSASSSVSYMHWYQQKSGTSPVLLIYTT
NNEENSSSENEEENTTNNEEEEEEESSNNSNNEEEEENTTSNNEEEEETT
```

where
 
N = no assoc. structure motif
E = beta strand
T = turn
G = 3/10 helix	B = beta bridge
S = bend
H = alpha helix



#### **Group Contributions**


Scott: implement forward/backward and Viterbi training
James: data curation/import, additional scripts
Rosemary: implement MCMC sampling
Everybody: analyze results

#### **References**


[1] Ranganathan, Ananth, and Frank Dellaert. "Inference in the space of topological maps: An MCMC-based approach." Intelligent Robots and Systems, 2004.(IROS 2004). Proceedings. 2004 IEEE/RSJ International Conference on. Vol. 2. IEEE, 2004.

