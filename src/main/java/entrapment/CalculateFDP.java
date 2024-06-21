/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

package entrapment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Set;

public class CalculateFDP {

  public static void main(String[] args) {
    if (args.length != 7) {
      System.out.println("Usage: java -cp EntrapBench.jar entrapment.CalculateFDP <fasta file path> <entrapment prefix> <result file path> <run precursor FDR> <global precursor FDR> <run protein group FDR> <global protein group FDR>");
      System.exit(1);
    }

    Path fastaPath = Paths.get(args[0]);
    String entrapmentPrefix = args[1];
    Path resultPath = Paths.get(args[2]);
    double runPrecursorFdrT = Double.parseDouble(args[3]);
    double globalPrecursorFdrT = Double.parseDouble(args[4]);
    double runPGFdrT = Double.parseDouble(args[5]);
    double globalPGFdrT = Double.parseDouble(args[6]);

    if (!Files.exists(fastaPath) || !Files.isReadable(fastaPath) || !Files.isRegularFile(fastaPath)) {
      System.out.println("The fasta file " + args[0] + " is not valid.");
      System.exit(1);
    }

    if (!Files.exists(resultPath) || !Files.isWritable(resultPath) || !Files.isRegularFile(resultPath)) {
      System.out.println("The result file " + args[2] + " is not valid.");
      System.exit(1);
    }

    try {
      Entry1 entry1 = summaryEntrapments(fastaPath, entrapmentPrefix);
      Entry2 entry2 = diannParser(resultPath, entrapmentPrefix, runPrecursorFdrT, globalPrecursorFdrT, runPGFdrT, globalPGFdrT);
      double r = (double) entry1.entrapmentProteinCount / (double) entry1.nonEntrapmentProteinCount;

      System.out.println("Non-entrapment proteins in the database: " + entry1.nonEntrapmentProteinCount);
      System.out.println("Entrapment proteins in the database: " + entry1.entrapmentProteinCount);
      System.out.println("r: " + r);
      System.out.println();
      if (entry2.thereAreDecoyScoreLargerThanTargetScore) {
        System.out.println("WARNING: There are decoy scores larger than target scores.");
      }
      System.out.println("Precursor level filtered with " + runPrecursorFdrT + " run q-value and " + globalPrecursorFdrT + " global q-value:");
      System.out.println("Target: " + entry2.targetPrecursorCount);
      System.out.println("Decoy (not accurate because DIA-NN does not report all decoys and the decoys are not FDR filtered): " + entry2.decoyPrecursorCount);
      System.out.println("Entrapment: " + entry2.entrapmentPrecursorCount);
      System.out.println("Decoy entrapment (not accurate because DIA-NN does not report all decoys and the decoys are not FDR filtered): " + entry2.decoyEntrapmentPrecursorCount);
      System.out.println("ET * (1 + 1/r) / (NT + ET): " + (entry2.entrapmentPrecursorCount * (1 + 1 / r) * 100.0 / (entry2.targetPrecursorCount + entry2.entrapmentPrecursorCount)) + "%");
      System.out.println("ET / (NT + ET): " + (entry2.entrapmentPrecursorCount * 100.0 / (entry2.targetPrecursorCount + entry2.entrapmentPrecursorCount)) + "%");
      System.out.println("ET * (1/r) / NT: " + (entry2.entrapmentPrecursorCount * (1 / r) * 100.0 / entry2.targetPrecursorCount) + "%");
      System.out.println();
      System.out.println("Protein level filtered with " + runPGFdrT + " run q-value and " + globalPGFdrT + " global q-value:");
      System.out.println("Target: " + entry2.targetProteinCount);
      System.out.println("Entrapment: " + entry2.entrapmentProteinCount);
      System.out.println("ET * (1 + 1/r) / (NT + ET): " + (entry2.entrapmentProteinCount * (1 + 1 / r) * 100.0 / (entry2.targetProteinCount + entry2.entrapmentProteinCount)) + "%");
      System.out.println("ET / (NT + ET): " + (entry2.entrapmentProteinCount * 100.0 / (entry2.targetProteinCount + entry2.entrapmentProteinCount)) + "%");
      System.out.println("ET * (1/r) / NT: " + (entry2.entrapmentProteinCount * (1 / r) * 100.0 / entry2.targetProteinCount) + "%");
    } catch (Exception ex) {
      ex.printStackTrace();
      System.exit(1);
    }

  }

  private static Entry1 summaryEntrapments(Path fastaPath, String entrapmentPrefix) throws Exception {
    String line;
    BufferedReader reader = new BufferedReader(new FileReader(fastaPath.toFile()));
    long entrapmentProteinCount = 0;
    long nonEntrapmentProteinCount = 0;
    while ((line = reader.readLine()) != null) {
      line = line.trim();
      if (line.startsWith(">")) {
        if (line.contains(entrapmentPrefix)) {
          ++entrapmentProteinCount;
        } else {
          ++nonEntrapmentProteinCount;
        }
      }
    }
    reader.close();

    return new Entry1(nonEntrapmentProteinCount, entrapmentProteinCount);
  }

  private static Entry2 diannParser(Path resultPath, String entrapmentPrefix, double runPrecursorFdrT, double globalPrecursorFdrT, double runPGFdrT, double globalPGFdrT) throws Exception {
    long targetPrecursorCount = 0, entrapmentPrecursorCount = 0, decoyPrecursorCount = 0, decoyEntrapmentPrecursorCount = 0;
    Set<String> targetProteins = new HashSet<>(), entrapmentProteins = new HashSet<>();
    String line;
    BufferedReader reader = new BufferedReader(new FileReader(resultPath.toFile()));
    int runColumnIdx = -1;
    int pgColumnIdx = -1;
    int cscoreColumnIdx = -1;
    int decoyCscoreColumnIdx = -1;
    int runPrecursorFdrColumnIdx = -1;
    int globalPrecursorFdrColumnIdx = -1;
    int runPGFdrColumnIdx = -1;
    int globalPGFdrColumnIdx = -1;
    boolean thereAreDecoyScoreLargerThanTargetScore = false;

    while ((line = reader.readLine()) != null) {
      line = line.trim();
      if (line.isEmpty()) {
        continue;
      }

      String[] parts = line.split("\t");
      if (line.startsWith("File.Name")) {
        for (int i = 0; i < parts.length; ++i) {
          if (parts[i].trim().equalsIgnoreCase("Run")) {
            runColumnIdx = i;
          } else if (parts[i].trim().equalsIgnoreCase("Protein.Group")) {
            pgColumnIdx = i;
          } else if (parts[i].trim().equalsIgnoreCase("CScore")) {
            cscoreColumnIdx = i;
          } else if (parts[i].trim().equalsIgnoreCase("Decoy.CScore")) {
            decoyCscoreColumnIdx = i;
          } else if (parts[i].trim().equalsIgnoreCase("Q.Value")) {
            runPrecursorFdrColumnIdx = i;
          } else if (parts[i].trim().equalsIgnoreCase("Global.Q.Value")) {
            globalPrecursorFdrColumnIdx = i;
          } else if (parts[i].trim().equalsIgnoreCase("PG.Q.Value")) {
            runPGFdrColumnIdx = i;
          } else if (parts[i].trim().equalsIgnoreCase("Global.PG.Q.Value")) {
            globalPGFdrColumnIdx = i;
          }
        }
        if (runColumnIdx < 0 || pgColumnIdx < 0 || cscoreColumnIdx < 0 || decoyCscoreColumnIdx < 0 || runPrecursorFdrColumnIdx < 0 || globalPrecursorFdrColumnIdx < 0 || runPGFdrColumnIdx < 0 || globalPGFdrColumnIdx < 0) {
          System.out.println("Some columns are missing in the result file: " + resultPath.toAbsolutePath());
          System.exit(1);
        }
      } else {
        String run = parts[runColumnIdx].trim();
        String pg = parts[pgColumnIdx].trim();
        double cscore = Double.parseDouble(parts[cscoreColumnIdx]);
        double decoyCscore = Double.parseDouble(parts[decoyCscoreColumnIdx]);
        double runPrecursorFdr = Double.parseDouble(parts[runPrecursorFdrColumnIdx]);
        double globalPrecursorFdr = Double.parseDouble(parts[globalPrecursorFdrColumnIdx]);
        double runPGFdr = Double.parseDouble(parts[runPGFdrColumnIdx]);
        double globalPGFdr = Double.parseDouble(parts[globalPGFdrColumnIdx]);

        if (cscore <= decoyCscore) {
          thereAreDecoyScoreLargerThanTargetScore = true;
        }

        String[] parts2 = pg.split(";");
        boolean isEntrapment = true;
        for (String p : parts2) {
          if (!p.contains(entrapmentPrefix)) { // As long as there is a non-entrapment protein, it is not an entrapment.
            isEntrapment = false;
            break;
          }
        }

        if (runPrecursorFdr < runPrecursorFdrT && globalPrecursorFdr < globalPrecursorFdrT) {
          if (isEntrapment) {
            ++entrapmentPrecursorCount;
            if (decoyCscore > 0) {
              ++decoyEntrapmentPrecursorCount;
            }
          } else {
            ++targetPrecursorCount;
            if (decoyCscore > 0) {
              ++decoyPrecursorCount;
            }
          }
        }

        if (runPGFdr < runPGFdrT && globalPGFdr < globalPGFdrT) {
          if (isEntrapment) {
            entrapmentProteins.add(run + "_" + pg);
          } else {
            targetProteins.add(run + "_" + pg);
          }
        }
      }
    }
    reader.close();

    return new Entry2(targetPrecursorCount, decoyPrecursorCount, entrapmentPrecursorCount, decoyEntrapmentPrecursorCount, targetProteins.size(), entrapmentProteins.size(), thereAreDecoyScoreLargerThanTargetScore);
  }


  static class Entry1 {

    final long nonEntrapmentProteinCount;
    final long entrapmentProteinCount;

    public Entry1(long nonEntrapmentProteinCount, long entrapmentProteinCount) {
      this.nonEntrapmentProteinCount = nonEntrapmentProteinCount;
      this.entrapmentProteinCount = entrapmentProteinCount;
    }
  }


  static class Entry2 {

    final double targetPrecursorCount;
    final double decoyPrecursorCount;
    final double entrapmentPrecursorCount;
    final double decoyEntrapmentPrecursorCount;
    final double targetProteinCount;
    final double entrapmentProteinCount;
    final boolean thereAreDecoyScoreLargerThanTargetScore;

    public Entry2(double targetPrecursorCount, double decoyPrecursorCount, double entrapmentPrecursorCount, double decoyEntrapmentPrecursorCount, double targetProteinCount, double entrapmentProteinCount, boolean thereAreDecoyScoreLargerThanTargetScore) {
      this.targetPrecursorCount = targetPrecursorCount;
      this.decoyPrecursorCount = decoyPrecursorCount;
      this.entrapmentPrecursorCount = entrapmentPrecursorCount;
      this.decoyEntrapmentPrecursorCount = decoyEntrapmentPrecursorCount;
      this.targetProteinCount = targetProteinCount;
      this.entrapmentProteinCount = entrapmentProteinCount;
      this.thereAreDecoyScoreLargerThanTargetScore = thereAreDecoyScoreLargerThanTargetScore;
    }
  }
}
