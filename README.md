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

You can predict on-target scores for sequences using the `predict_seq` function, specifying either
[Hsu2013](https://www.nature.com/articles/nbt.2647) or
[Chen2013](https://www.sciencedirect.com/science/article/pii/S0092867413015316?via%3Dihub)
as the tracrRNA to score with

```
predict_seq(context_seqs, sequence_tracr='Hsu2013')
```




    array([0.28923641, 0.78714099])


