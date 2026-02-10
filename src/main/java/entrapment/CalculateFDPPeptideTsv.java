/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package entrapment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class CalculateFDPPeptideTsv {

  public static void main(String[] args) {
    if (args.length != 3) {
      System.out.println("Usage: java -cp EntrapBench.jar entrapment.CalculateFDPPeptideTsv <fasta file path> <entrapment style> <peptide.tsv file path>");
      System.exit(1);
    }

    Path fastaPath = Paths.get(args[0]);
    int entrapmentStyle = Integer.parseInt(args[1]);
    Path resultPath = Paths.get(args[2]);

    if (entrapmentStyle != 0 && entrapmentStyle != 1) {
      System.out.println("Unknown entrapment style.");
      System.exit(1);
    }

    String entrapmentMarker = entrapmentStyle == 0 ? "entrapment_" : "_p_target";

    if (!Files.exists(fastaPath) || !Files.isReadable(fastaPath) || !Files.isRegularFile(fastaPath)) {
      System.out.println("The fasta file " + args[0] + " is not valid.");
      System.exit(1);
    }

    if (!Files.exists(resultPath) || !Files.isReadable(resultPath) || !Files.isRegularFile(resultPath)) {
      System.out.println("The result file " + args[2] + " is not valid.");
      System.exit(1);
    }

    try {
      CalculateFDP.Entry1 entry1 = summaryEntrapments(fastaPath, entrapmentMarker);
      PeptideTsvResult entry2 = peptideTsvParser(resultPath, entrapmentMarker);
      double r = (double) entry1.entrapmentProteinCount / (double) entry1.nonEntrapmentProteinCount;

      System.out.println("Non-entrapment proteins in the database: " + entry1.nonEntrapmentProteinCount);
      System.out.println("Entrapment proteins in the database: " + entry1.entrapmentProteinCount);
      System.out.println("r: " + r);
      System.out.println();
      System.out.println("Peptide level:");
      System.out.println("Target: " + entry2.targetPeptideCount);
      System.out.println("Entrapment: " + entry2.entrapmentPeptideCount);
      System.out.println("ET * (1 + 1/r) / (NT + ET): " + (entry2.entrapmentPeptideCount * (1 + 1 / r) * 100.0 / (entry2.targetPeptideCount + entry2.entrapmentPeptideCount)) + "%");
      System.out.println("ET / (NT + ET): " + (entry2.entrapmentPeptideCount * 100.0 / (entry2.targetPeptideCount + entry2.entrapmentPeptideCount)) + "%");
      System.out.println("ET * (1/r) / NT: " + (entry2.entrapmentPeptideCount * (1 / r) * 100.0 / entry2.targetPeptideCount) + "%");
    } catch (Exception ex) {
      ex.printStackTrace();
      System.exit(1);
    }
  }

  private static CalculateFDP.Entry1 summaryEntrapments(Path fastaPath, String entrapmentMarker) throws Exception {
    String line;
    BufferedReader reader = new BufferedReader(new FileReader(fastaPath.toFile()));
    long entrapmentProteinCount = 0;
    long nonEntrapmentProteinCount = 0;
    while ((line = reader.readLine()) != null) {
      line = line.trim();
      if (line.startsWith(">")) {
        if (line.contains(entrapmentMarker)) {
          ++entrapmentProteinCount;
        } else {
          ++nonEntrapmentProteinCount;
        }
      }
    }
    reader.close();

    return new CalculateFDP.Entry1(nonEntrapmentProteinCount, entrapmentProteinCount);
  }

  private static PeptideTsvResult peptideTsvParser(Path resultPath, String entrapmentMarker) throws Exception {
    long targetPeptideCount = 0, entrapmentPeptideCount = 0;
    String line;
    BufferedReader reader = new BufferedReader(new FileReader(resultPath.toFile()));
    int proteinColumnIdx = -1;
    int mappedProteinsColumnIdx = -1;

    while ((line = reader.readLine()) != null) {
      line = line.trim();
      if (line.isEmpty()) {
        continue;
      }

      String[] parts = line.split("\t", -1);
      if (line.startsWith("Peptide\t")) {
        for (int i = 0; i < parts.length; ++i) {
          String col = parts[i].trim();
          if (col.equalsIgnoreCase("Protein")) {
            proteinColumnIdx = i;
          } else if (col.equalsIgnoreCase("Mapped Proteins")) {
            mappedProteinsColumnIdx = i;
          }
        }
        if (proteinColumnIdx < 0) {
          System.out.println("Protein column is missing in the result file: " + resultPath.toAbsolutePath());
          System.exit(1);
        }
      } else {
        String protein = parts[proteinColumnIdx].trim();
        String mappedProteins = (mappedProteinsColumnIdx >= 0 && mappedProteinsColumnIdx < parts.length) ? parts[mappedProteinsColumnIdx].trim() : "";

        // Combine primary protein and mapped proteins for entrapment check
        boolean isEntrapment = protein.contains(entrapmentMarker);
        if (isEntrapment && !mappedProteins.isEmpty()) {
          String[] mapped = mappedProteins.split(",");
          for (String m : mapped) {
            if (!m.trim().contains(entrapmentMarker)) {
              isEntrapment = false;
              break;
            }
          }
        }

        if (isEntrapment) {
          ++entrapmentPeptideCount;
        } else {
          ++targetPeptideCount;
        }
      }
    }
    reader.close();

    return new PeptideTsvResult(targetPeptideCount, entrapmentPeptideCount);
  }


  static class PeptideTsvResult {

    final long targetPeptideCount;
    final long entrapmentPeptideCount;

    public PeptideTsvResult(long targetPeptideCount, long entrapmentPeptideCount) {
      this.targetPeptideCount = targetPeptideCount;
      this.entrapmentPeptideCount = entrapmentPeptideCount;
    }
  }
}
