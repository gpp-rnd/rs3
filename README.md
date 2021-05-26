# Rule Set 3
> Python package to predict the activity of CRISPR sgRNA sequences using Rule Set 3


## Install

`pip install git+ssh://git@github.com/gpp-rnd/rs3.git`

## How to use

Import packages

```
from rs3.seq import predict_seq
```

Create a list of context sequences you want to predict

```
context_seqs = ['GACGAAAGCGACAACGCGTTCATCCGGGCA', 'AGAAAACACTAGCATCCCCACCCGCGGACT']
```

Predict on sequences using the `predict_seq` function

```
predict_seq(context_seqs)
```




    array([0.28923641, 0.78714099])


