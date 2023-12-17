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
import java.io.BufferedWriter;
import java.io.FileReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class DiannEntrapmentQValue {

  private static final double binSize = 0.0000001;

  public static void main(String[] args) {
    if (args.length != 7) {
      System.out.println("Usage: java -cp EntrapBench.jar entrapment.DiannEntrapmentQValue <entrapment prefix> <run-wise precursor q-value threshold> <global precursor q-value threshold> <run-wise protein q-value threshold> <global protein q-value threshold> <result file path> <output file path>\n");
      System.exit(1);
    }

    String entrapmentPrefix = args[0];
    double runPrecursorQValueT = Double.parseDouble(args[1]);
    double globalPrecursorQValueT = Double.parseDouble(args[2]);
    double runPGQValueT = Double.parseDouble(args[3]);
    double globalPGQValueT = Double.parseDouble(args[4]);
    Path resultPath = Paths.get(args[5]);
    Path outputPath = Paths.get(args[6]);

    if (!Files.exists(resultPath) || !Files.isWritable(resultPath) || !Files.isRegularFile(resultPath)) {
      System.out.println("The result file " + args[1] + " is not valid.");
      System.exit(1);
    }

    try {
      Entry entry = calculate(resultPath, entrapmentPrefix, runPrecursorQValueT, globalPrecursorQValueT, runPGQValueT, globalPGQValueT);

      BufferedWriter writer = Files.newBufferedWriter(outputPath);
      writer.write("cscore_threshold,target_count,decoy_count,entrapment_target_count,entrapment_decoy_count,reported_run_precursor_Q_value,reported_global_precursor_Q_value,entrapment_Q_value\n");
      for (int i = entry.entrapmentTargetCounts.length - 1; i >= 0; --i) {
        if (Double.isNaN(entry.reportedRunQValues[i]) || Double.isNaN(entry.reportedGlobalQValues[i])) {
          continue;
        }
        writer.write((i * binSize) + "," + entry.targetCounts[i] + "," + entry.decoyCounts[i] + "," + entry.entrapmentTargetCounts[i] + "," + entry.entrapmentDecoyCounts[i] + "," + entry.reportedRunQValues[i] + "," + entry.reportedGlobalQValues[i] + "," + entry.entrapmentQValues[i] + "\n");
      }
      writer.close();

      System.out.println("Precursor level filtered with " + runPrecursorQValueT + " run q-value and " + globalPrecursorQValueT + " global q-value:");
      System.out.println("Target: " + entry.reportedTargetPrecursorCount);
      System.out.println("Decoy (not accurate because DIA-NN does not report all decoys and the decoys are not FDR filtered): " + entry.reportedDecoyPrecursorCount);
      System.out.println("Entrapment: " + entry.reportedEntrapmentPrecursorCount);
      System.out.println("Entrapment decoy (not accurate because DIA-NN does not report all decoys and the decoys are not FDR filtered): " + entry.reportedEntrapmentDecoyPrecursorCount);
      System.out.println("(ND + ET) / (NT + ET) (not accurate because DIA-NN does not report all ND and the ND are not FDR filtered): " + ((entry.reportedDecoyPrecursorCount + entry.reportedEntrapmentPrecursorCount) * 100.0 / (entry.reportedTargetPrecursorCount + entry.reportedEntrapmentPrecursorCount)) + "%");
      System.out.println("ET / (NT + ET): " + (entry.reportedEntrapmentPrecursorCount * 100.0 / (entry.reportedTargetPrecursorCount + entry.reportedEntrapmentPrecursorCount)) + "%");
      System.out.println("DIA-NN reported " + (entry.reportedTargetPrecursorCount + entry.reportedEntrapmentPrecursorCount) + " precursors.");
      System.out.println("With global entrapment q-value (ET / (NT + ET)) threshold = " + (Math.min(runPrecursorQValueT, globalPrecursorQValueT) * 100) + "%, there are " + entry.entrapmentQValueFilteredPrecursorCount + " precursors.");

      System.out.println();
      System.out.println("Protein level filtered with " + runPGQValueT + " run q-value and " + globalPGQValueT + " global q-value:");
      System.out.println("Target: " + entry.targetProteinCount);
      System.out.println("Entrapment: " + entry.entrapmentProteinCount);
      System.out.println("ET / (NT + ET): " + (entry.entrapmentProteinCount * 100.0 / (entry.targetProteinCount + entry.entrapmentProteinCount)) + "%");
    } catch (Exception ex) {
      ex.printStackTrace();
      System.exit(1);
    }
  }

  private static Entry calculate(Path resultPath, String entrapmentPrefix, double runPrecursorQValueT, double globalPrecursorQValueT, double runPGQValueT, double globalPGQValueT) throws Exception {
    long[] targetCounts = new long[(int) (1 / binSize) + 1];
    long[] decoyCounts = new long[(int) (1 / binSize) + 1];
    long[] entrapmentTargetCounts = new long[(int) (1 / binSize) + 1];
    long[] entrapmentDecoyCounts = new long[(int) (1 / binSize) + 1];
    double[] reportedRunQValues = new double[(int) (1 / binSize) + 1];
    Arrays.fill(reportedRunQValues, Double.NaN);
    double[] reportedGlobalQValues = new double[(int) (1 / binSize) + 1];
    Arrays.fill(reportedGlobalQValues, Double.NaN);

    long reportedTargetPrecursorCount = 0;
    long reportedDecoyPrecursorCount = 0;
    long reportedEntrapmentPrecursorCount = 0;
    long reportedEntrapmentDecoyPrecursorCount = 0;
    Set<String> entrapmentProteins = new HashSet<>();
    Set<String> targetProteins = new HashSet<>();

    String line;
    BufferedReader reader = new BufferedReader(new FileReader(resultPath.toFile()));
    int runColumnIdx = -1;
    int pgColumnIdx = -1;
    int cscoreColumnIdx = -1;
    int decoyCscoreColumnIdx = -1;
    int runPrecursorQValueColumnIdx = -1;
    int globalPrecursorQValueColumnIdx = -1;
    int runPGQValueColumnIdx = -1;
    int globalPGQValueColumnIdx = -1;
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
            runPrecursorQValueColumnIdx = i;
          } else if (parts[i].trim().equalsIgnoreCase("Global.Q.Value")) {
            globalPrecursorQValueColumnIdx = i;
          } else if (parts[i].trim().equalsIgnoreCase("PG.Q.Value")) {
            runPGQValueColumnIdx = i;
          } else if (parts[i].trim().equalsIgnoreCase("Global.PG.Q.Value")) {
            globalPGQValueColumnIdx = i;
          }
        }
        if (runColumnIdx < 0 || pgColumnIdx < 0 || cscoreColumnIdx < 0 || decoyCscoreColumnIdx < 0 || runPrecursorQValueColumnIdx < 0 || globalPrecursorQValueColumnIdx < 0 || runPGQValueColumnIdx < 0 || globalPGQValueColumnIdx < 0) {
          System.out.println("Some columns are missing in the result file: " + resultPath.toAbsolutePath());
          System.exit(1);
        }
      } else {
        String run = parts[runColumnIdx].trim();
        String pg = parts[pgColumnIdx].trim();
        double cscore = Double.parseDouble(parts[cscoreColumnIdx]);
        double decoyCscore = Double.parseDouble(parts[decoyCscoreColumnIdx]);
        double runPrecursorQValue = Double.parseDouble(parts[runPrecursorQValueColumnIdx]);
        double globalPrecursorQValue = Double.parseDouble(parts[globalPrecursorQValueColumnIdx]);
        double runPGQValue = Double.parseDouble(parts[runPGQValueColumnIdx]);
        double globalPGQValue = Double.parseDouble(parts[globalPGQValueColumnIdx]);

        String[] parts2 = pg.split(";");
        boolean isEntrapment = true;
        for (String p : parts2) {
          if (!p.contains(entrapmentPrefix)) { // As long as there is a non-entrapment protein, it is not an entrapment.
            isEntrapment = false;
            break;
          }
        }

        if (runPrecursorQValue < runPrecursorQValueT && globalPrecursorQValue < globalPrecursorQValueT) {
          if (isEntrapment) {
            ++reportedEntrapmentPrecursorCount;
            if (decoyCscore > 0) {
              ++reportedEntrapmentDecoyPrecursorCount;
            }
          } else {
            ++reportedTargetPrecursorCount;
            if (decoyCscore > 0) {
              ++reportedDecoyPrecursorCount;
            }
          }
        }

        if (runPGQValue < runPGQValueT && globalPGQValue < globalPGQValueT) {
          if (isEntrapment) {
            entrapmentProteins.add(run + "_" + pg);
          } else {
            targetProteins.add(run + "_" + pg);
          }
        }

        if (Double.isNaN(reportedRunQValues[(int) (cscore / binSize)])) {
          reportedRunQValues[(int) (cscore / binSize)] = runPrecursorQValue;
        } else if (Math.abs(reportedRunQValues[(int) (cscore / binSize)] - runPrecursorQValue) > 1e-6) {
          // System.out.println("There are different reported run q-values for the same cscore: " + cscore + " : " + reportedRunQValues[(int) (cscore / binSize)] + " vs " + runPrecursorQValue);
          reportedRunQValues[(int) (cscore / binSize)] = Math.min(reportedRunQValues[(int) (cscore / binSize)], runPrecursorQValue);
        }

        if (Double.isNaN(reportedGlobalQValues[(int) (cscore / binSize)])) {
          reportedGlobalQValues[(int) (cscore / binSize)] = globalPrecursorQValue;
        } else if (Math.abs(reportedGlobalQValues[(int) (cscore / binSize)] - globalPrecursorQValue) > 1e-6) {
          // System.out.println("There are different reported global q-values for the same cscore: " + cscore + " : " + reportedGlobalQValues[(int) (cscore / binSize)] + " vs " + globalPrecursorQValue);
          reportedGlobalQValues[(int) (cscore / binSize)] = Math.min(reportedGlobalQValues[(int) (cscore / binSize)], globalPrecursorQValue);
        }

        if (isEntrapment) {
          ++entrapmentTargetCounts[(int) (cscore / binSize)];
          if (decoyCscore > 0) {
            ++entrapmentDecoyCounts[(int) (decoyCscore / binSize)];
          }
        } else {
          ++targetCounts[(int) (cscore / binSize)];
          if (decoyCscore > 0) {
            ++decoyCounts[(int) (decoyCscore / binSize)];
          }
        }
      }
    }
    reader.close();

//    fillMissingValues(reportedRunQValues);
//    fillMissingValues(reportedGlobalQValues);

    double[] entrapmentQValues = calculateQValue(targetCounts, decoyCounts, entrapmentTargetCounts, entrapmentDecoyCounts, 3);
    long entrapmentQValueFilteredPrecursors = filterPrecursors(entrapmentQValues, Math.min(runPrecursorQValueT, globalPrecursorQValueT), targetCounts, entrapmentTargetCounts);

    return new Entry(targetCounts, decoyCounts, entrapmentTargetCounts, entrapmentDecoyCounts, reportedRunQValues, reportedGlobalQValues, entrapmentQValues, reportedTargetPrecursorCount, reportedDecoyPrecursorCount, reportedEntrapmentPrecursorCount, reportedEntrapmentDecoyPrecursorCount, entrapmentQValueFilteredPrecursors, entrapmentProteins.size(), targetProteins.size());
  }

  private static void fillMissingValues(double[] qValues) {
    if (Double.isNaN(qValues[qValues.length - 1])) {
      qValues[qValues.length - 1] = 0;
    }
    for (int i = qValues.length - 2; i >= 0; --i) {
      if (Double.isNaN(qValues[i])) {
        qValues[i] = qValues[i + 1];
      }
    }
  }

  private static double[] calculateQValue(long[] targetCounts, long[] decoyCounts, long[] entrapmentTargetCounts, long[] entrapmentDecoyCounts, int equation) {
    if (targetCounts.length != decoyCounts.length || targetCounts.length != entrapmentTargetCounts.length || targetCounts.length != entrapmentDecoyCounts.length) {
      System.out.println("The length of the target and decoy arrays are not equal.");
      System.exit(1);
    }

    long decoyCount = 0;
    long targetCount = 0;
    double[] fdrs = new double[targetCounts.length];
    double fdr;
    for (int i = targetCounts.length - 1; i >= 0; --i) {
      if (equation == 1) { // not accurate because DIA-NN does not report all decoys and the decoys are not FDR filtered
        decoyCount += decoyCounts[i] + entrapmentDecoyCounts[i];
        targetCount += targetCounts[i] + entrapmentTargetCounts[i];
      } else if (equation == 2) { // not accurate because DIA-NN does not report all decoys and the decoys are not FDR filtered
        decoyCount += decoyCounts[i] + entrapmentTargetCounts[i];
        targetCount += targetCounts[i] + entrapmentTargetCounts[i];
      } else if (equation == 3) {
        decoyCount += entrapmentTargetCounts[i];
        targetCount += targetCounts[i] + entrapmentTargetCounts[i];
      } else {
        System.out.println("The equation " + equation + " is not supported.");
        System.exit(1);
      }

      if (targetCount == 0) {
        fdr = 0;
      } else {
        fdr = (double) decoyCount / (double) targetCount;
      }
      fdrs[i] = Math.min(fdr, 1);
    }

    double[] qValues = new double[fdrs.length];
    double lastQValue = fdrs[0];
    qValues[0] = fdrs[0];
    for (int i = 1; i < fdrs.length; ++i) {
      if (fdrs[i] > lastQValue) {
        qValues[i] = lastQValue;
      } else {
        qValues[i] = fdrs[i];
        lastQValue = fdrs[i];
      }
    }

    return qValues;
  }

  private static long filterPrecursors(double[] qValues, double qValueT, long[] targetCounts, long[] entrapmentTargetCounts) {
    long count = 0;
    for (int i = qValues.length - 1; i >= 0; --i) {
      if (qValues[i] < qValueT) {
        count += targetCounts[i] + entrapmentTargetCounts[i];
      } else {
        break;
      }
    }

    return count;
  }


  static class Entry {

    final long[] targetCounts;
    final long[] decoyCounts;
    final long[] entrapmentTargetCounts;
    final long[] entrapmentDecoyCounts;
    final double[] reportedRunQValues;
    final double[] reportedGlobalQValues;
    final double[] entrapmentQValues;
    final long reportedTargetPrecursorCount;
    final long reportedDecoyPrecursorCount;
    final long reportedEntrapmentPrecursorCount;
    final long reportedEntrapmentDecoyPrecursorCount;
    final long entrapmentQValueFilteredPrecursorCount;
    final long entrapmentProteinCount;
    final long targetProteinCount;

    public Entry(long[] targetCounts, long[] decoyCounts, long[] entrapmentTargetCounts, long[] entrapmentDecoyCounts, double[] reportedRunQValues, double[] reportedGlobalQValues, double[] entrapmentQValues, long reportedTargetPrecursorCount, long reportedDecoyPrecursorCount, long reportedEntrapmentPrecursorCount, long reportedEntrapmentDecoyPrecursorCount, long entrapmentQValueFilteredPrecursorCount, long entrapmentProteinCount, long targetProteinCount) {
      this.targetCounts = targetCounts;
      this.decoyCounts = decoyCounts;
      this.entrapmentTargetCounts = entrapmentTargetCounts;
      this.entrapmentDecoyCounts = entrapmentDecoyCounts;
      this.reportedRunQValues = reportedRunQValues;
      this.reportedGlobalQValues = reportedGlobalQValues;
      this.entrapmentQValues = entrapmentQValues;
      this.reportedTargetPrecursorCount = reportedTargetPrecursorCount;
      this.reportedDecoyPrecursorCount = reportedDecoyPrecursorCount;
      this.reportedEntrapmentPrecursorCount = reportedEntrapmentPrecursorCount;
      this.reportedEntrapmentDecoyPrecursorCount = reportedEntrapmentDecoyPrecursorCount;
      this.entrapmentQValueFilteredPrecursorCount = entrapmentQValueFilteredPrecursorCount;
      this.entrapmentProteinCount = entrapmentProteinCount;
      this.targetProteinCount = targetProteinCount;
    }
  }
}
