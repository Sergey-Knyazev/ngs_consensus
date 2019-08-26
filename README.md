# ngs_consensus
This tool extracts consensus from NGS reads of viral samples. It resolves ambiguities if asked or keep the ambiguities that occur in reads with frequency above threshold.
The threshold can be 5%, 10%, 20%, or any other number.

## running
**miseq_consensus.ipynb** takes raw reads as an input and generates four consensus sequences: one with ambuiguities resolved, and three with ambiguities with frequency above thresholds of 5%, 10%, and 20%.

To run miseq_consensus.ipynb, modify config file *miseq_consensus.config.yaml* and execute the following command from terminal:
```
jupyter nbconvert --to notebook --execute miseq_consensus.ipynb
```
