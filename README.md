# EntrapBench
Generate entrapment database

### Generate a target+entrapment database given a target database
Given each protein in a target database, digest it, shuffle the peptides, and then put the peptides back into proteins. Each peptide is shuffled at most 10 times to get a unique sequence. Depending on the parameter, one target protein can generate multiple entrapment proteins.

Usage:
```shell
java -cp EntrapBench.jar entrapment.GenerateDatabase <UniProt fasta file path> <cut sites> <protect sites> <cleavage from C-term: 0=false, 1 = true> <number of entrapment proteins for each target protein> <entrapment prefix> <add prefix>
Example: java -cp EntrapBench.jar entrapment.GenerateDatabase uniprot_human.fasta KR P 1 10 entrapment 0 # Each target protein generates 10 shuffled entrapment proteins.
```

__Note:__ the "target" here is different from the term "target" in the target-decoy database searching approach. To use this target+entrapment database in the target-decoy approach, need to generate decoy proteins for both target and entrapment proteins.
