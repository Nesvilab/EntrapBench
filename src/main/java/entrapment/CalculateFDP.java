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

public class CalculateFDP {

  public static void main(String[] args) {
    if (args.length != 8) {
      System.out.println("Usage: java -cp EntrapBench.jar entrapment.CalculateFDP <fasta file path> <decoy prefix> <entrapment prefix> <result file path> <run precursor FDR> <global precursor FDR> <run protein group FDR> <global protein group FDR>\n"
          + "If there are no decoys, put null to the decoy prefix.");
      System.exit(1);
    }

    Path fastaPath = Paths.get(args[0]);
    String decoyPrefix = args[1].equals("null") ? null : args[1];
    String prefix = args[2];
    Path resultPath = Paths.get(args[3]);
    double runPrecursorFdrT = Double.parseDouble(args[4]);
    double globalPrecursorFdrT = Double.parseDouble(args[5]);
    double runPGFdrT = Double.parseDouble(args[6]);
    double globalPGFdrT = Double.parseDouble(args[7]);

    if (!Files.exists(fastaPath) || !Files.isReadable(fastaPath) || !Files.isRegularFile(fastaPath)) {
      System.out.println("The fasta file " + args[0] + " is not valid.");
      System.exit(1);
    }

    if (!Files.exists(resultPath) || !Files.isWritable(resultPath) || !Files.isRegularFile(resultPath)) {
      System.out.println("The result file " + args[2] + " is not valid.");
      System.exit(1);
    }

    try {
      double[] dbTE = getTargetEntrapmentProteins(fastaPath, decoyPrefix, prefix);
      double[] targetEntrapmentPrecursors = diannParser(resultPath, prefix, runPrecursorFdrT, globalPrecursorFdrT, runPGFdrT, globalPGFdrT);

      double fdp = (dbTE[0] * targetEntrapmentPrecursors[1]) / (dbTE[1] * targetEntrapmentPrecursors[0]);

      System.out.println("Target proteins in the database: " + dbTE[0]);
      System.out.println("Entrapment proteins in the database: " + dbTE[1]);
      System.out.println("Target precursors: " + targetEntrapmentPrecursors[0]);
      System.out.println("Entrapment precursors: " + targetEntrapmentPrecursors[1]);
      System.out.println("FDP: " + (fdp * 100) + "%");
    } catch (Exception ex) {
      ex.printStackTrace();
      System.exit(1);
    }

  }

  private static double[] getTargetEntrapmentProteins(Path fastaPath, String decoyPrefix, String prefix) throws Exception {
    String line;
    BufferedReader reader = new BufferedReader(new FileReader(fastaPath.toFile()));
    long entrapmentProteinCount = 0;
    long targetProteinCount = 0;
    while ((line = reader.readLine()) != null) {
      line = line.trim();
      if (line.startsWith(">") && (decoyPrefix == null || !line.startsWith(">" + decoyPrefix))) {
        if (line.contains(prefix)) {
          ++entrapmentProteinCount;
        } else {
          ++targetProteinCount;
        }
      }
    }
    reader.close();

    return new double[]{targetProteinCount, entrapmentProteinCount};
  }

  private static double[] diannParser(Path resultPath, String prefix, double runPrecursorFdrT, double globalPrecursorFdrT, double runPGFdrT, double globalPGFdrT) throws Exception {
    long targetCount = 0, entrapmentCount = 0;
    String line;
    BufferedReader reader = new BufferedReader(new FileReader(resultPath.toFile()));
    int pgColumnIdx = -1;
    int runPrecursorFdrColumnIdx = -1;
    int globalPrecursorFdrColumnIdx = -1;
    int runPGFdrColumnIdx = -1;
    int globalPGFdrColumnIdx = -1;
    while ((line = reader.readLine()) != null) {
      line = line.trim();
      if (line.isEmpty()) {
        continue;
      }

      String[] parts = line.split("\t");
      if (line.startsWith("File.Name")) {
        for (int i = 0; i < parts.length; ++i) {
          if (parts[i].trim().equalsIgnoreCase("Protein.Group")) {
            pgColumnIdx = i;
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
        if (pgColumnIdx < 0 || runPrecursorFdrColumnIdx < 0 || globalPrecursorFdrColumnIdx < 0 || runPGFdrColumnIdx < 0 || globalPGFdrColumnIdx < 0) {
          System.out.println("Some columns are missing in the result file: " + resultPath.toAbsolutePath());
          System.exit(1);
        }
      } else {
        String pg = parts[pgColumnIdx].trim();
        double runPrecursorFdr = Double.parseDouble(parts[runPrecursorFdrColumnIdx]);
        double globalPrecursorFdr = Double.parseDouble(parts[globalPrecursorFdrColumnIdx]);
        double runPGFdr = Double.parseDouble(parts[runPGFdrColumnIdx]);
        double globalPGFdr = Double.parseDouble(parts[globalPGFdrColumnIdx]);

        if (runPrecursorFdr < runPrecursorFdrT && globalPrecursorFdr < globalPrecursorFdrT && runPGFdr < runPGFdrT && globalPGFdr < globalPGFdrT) {
          String[] parts2 = pg.split(";");
          boolean isEntrapment = true;
          for (String p : parts2) {
            if (!p.startsWith(prefix)) { // As long as there is a non-entrapment protein, it is not entrapment.
              isEntrapment = false;
              break;
            }
          }
          if (isEntrapment) {
            ++entrapmentCount;
          } else {
            ++targetCount;
          }
        }
      }
    }
    reader.close();

    return new double[]{targetCount, entrapmentCount};
  }
}
