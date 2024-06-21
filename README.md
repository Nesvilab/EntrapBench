# EntrapBench
Generate entrapment database and calculate false discovery proportion (FDP)

### Generate a target+entrapment database given a target database
Given each protein in a target database, digest it, shuffle the peptides, and then put the peptides back into proteins. Each peptide is shuffled at most 10 times to get a unique sequence. Depending on the parameter, one target protein can generate multiple entrapment proteins.

Usage:
```shell
java -cp EntrapBench.jar entrapment.GenerateDatabase <UniProt fasta file path> <cut sites> <protect sites> <cleavage from C-term: 0=false, 1 = true> <number of entrapment proteins for each target protein> <entrapment prefix> <add prefix>
Example: java -cp EntrapBench.jar entrapment.GenerateDatabase uniprot_human.fasta KR P 1 10 entrapment 0 # Each target protein generates 10 shuffled entrapment proteins.
```

### Calculate false discovery proportion (FDP)
Given a target+entrapment database and DIA-NN's `report.tsv`, calculate the false discovery proportion related estimations using the equations in [Wen et al. (2024)](https://doi.org/10.1101/2024.06.01.596967)

"combined" method: $$FDP = \frac{E \times (1 + 1/r)}{T + E}$$

Lower bound: $$FDP = \frac{E}{T + E}$$

"sample" method: $$FDP = \frac{E \times (1/r)}{T}$$

where $T$ is the number of target proteins/precursors in the result, $E$ is the number of entrapment proteins/precursors in the result, $r$ is the ratio of entrapment and target proteins in the database.

`_Disclaimer: The equation may be slightly different depending on different target-decoy approaches and the interpretations of the false matches._`

Usage:
```shell
java -cp EntrapBench.jar entrapment.CalculateFDP <fasta file path> <entrapment prefix> <result file path> <run precursor FDR> <global precursor FDR> <run protein group FDR> <global protein group FDR>
Example: java -cp EntrapBench.jar entrapment.CalculateFDP uniprot_human.fasta entrapment report.tsv 0.01 0.01 0.01 0.01
```

```shell
java -cp EntrapBench.jar entrapment.DiannEntrapmentQValue <entrapment prefix> <entrapment to target ratio> <run-wise precursor q-value threshold> <global precursor q-value threshold> <run-wise protein q-value threshold> <global protein q-value threshold> <result file path> <output file path>
Example: java -cp EntrapBench.jar entrapment.DiannEntrapmentQValue entrapment 1 0.01 0.01 0.01 0.01 report.tsv entrapment_q_values.csv
```

__Note:__ the "target" here is different from the term "target" in the target-decoy database searching approach. To use this target+entrapment database in the target-decoy approach, need to generate decoy proteins (beforehand or on-the-fly by the tool itself) for both target and entrapment proteins.
