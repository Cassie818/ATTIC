# ATTIC 

## Introduction ##

ATTIC is a novel stacked-ensemble learning model to identify A-to-I editing sites in three species accurately. We first comprehensively evaluated 37 RNA sequence-derived features combined with 14 popular machine learning algorithms. Then we selected the optimal base models to build a series of stacked ensemble models. The final ATTIC framework was developed upon the optimal models improved by feature selection strategy for specific species. Extensive cross-validation and independent tests demonstrate that ATTIC outperforms state-of-the-art tools for predicting A-to-I editing sites.

![](/Users/cassie/Desktop/Submit_new/Figures/Figure 1.png)

## Setting the environment

The environment on our computer is as follows:

* Python 3.8
* pycaret 2.3.5
* Pandas 1.4.1
* NumPy 1.21.2
* scikit-learn 0.23.2

## Usage

```python
python ATTIC.py -i "input fasta file" -s "species" ['H.sapiens','M.musculus','D.melanogaster'] -o "output file" -t "threshold"
```

For example:

```python
python ATTIC.py -i "Example/input.fasta" -s 'H.sapiens' -o Example -t 0.5
```

## Notes

The predicted sequence length of *H. sapiens* and *D. melanogaster* must be greater than 51, while the expected sequence length of *M. musculus* must be greater than 41.

## Cite

Please cite our paper if you use this code in your own work.

## Contact

Feel free to contact us if you have any questions. You can send your queries and suggestions to the developers:
[fuyi.li@nwafu.edu.cn ](mailto:fuyi.li@nwafu.edu.cn), [Jiangning.song@monash.edu ](mailto:Jiangning.song@monash.edu)and [lachlan.coin@unimelb.edu.au](mailto:lachlan.coin@unimelb.edu.au) .








