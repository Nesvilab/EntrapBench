# EntrapBench
Generate entrapment database and calculate false discovery proportion (FDP)

### Generate a target+entrapment database given a target database
Given each protein in a target database, digest it, shuffle the peptides, and then put the peptides back into proteins. Each peptide is shuffled at most 10 times to get a unique sequence. Depending on the parameter, one target protein can generate multiple entrapment proteins.

Usage:
```shell
java -cp EntrapBench.jar entrapment.GenerateDatabase <UniProt fasta file path> <cut sites> <protect sites> <cleavage from C-term: 0=false, 1 = true> <number of entrapment proteins for each target protein> <entrapment prefix>
Example: java -cp EntrapBench.jar entrapment.GenerateDatabase uniprot_human.fasta KR P 1 10 entrapment # Each target protein generates 10 shuffled entrapment proteins.
```

### Calculate false discovery proportion (FDP)
Given a target+entrapment database and DIA-NN's `report.tsv`, calculate the false discovery proportion using the equation

$$FDP = \frac{T \times e}{E \times t}$$

where $T$ is the number of target proteins in the database, $E$ is the number of entrapment proteins in the database, $t$ is the number of target precursors in the result, and $e$ is the number of entrapment precursors in the result.

_Disclaimer: The equation may be slightly different depending on different target-decoy approaches and the interpretations of the false matches. The equation listed above is the one I could find from the literature._

Usage:
```shell
java -cp EntrapBench.jar entrapment.CalculateFDP <fasta file path> <decoy prefix> <entrapment prefix> <result file path> <run precursor FDR> <global precursor FDR> <run protein group FDR> <global protein group FDR>
Example: java -cp EntrapBench.jar entrapment.CalculateFDP uniprot_human.fasta null entrapment diann-output/report.tsv 0.01 0.01 0.01 0.01
```

__Note:__ the "target" here is different from the term "target" in the target-decoy database searching approach. To use this target+entrapment database in the target-decoy approach, need to generate decoy proteins for both target and entrapment proteins.
