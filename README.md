# ATTIC

## Setting the environment

The environment on our computer is as follows:
* Python 3.8
* pycaret 2.3.5
* Pandas 1.4.1
* NumPy 1.21.2
* scikit-learn 0.23.2

## Usage

```
python ATTIC.py -i "input fasta file" -s "species" ['H.sapiens','M.musculus','D.melanogaster'] -o "output file" -t "threshold"
```

For example:
```python
python ATTIC.py -i "Example/input.fasta" -s 'H.sapiens' -o Example -t 0.5
```
## Cite

Please cite our paper if you use this code in your own work.

## Note
The length of the input sequence should be greater than 51.
