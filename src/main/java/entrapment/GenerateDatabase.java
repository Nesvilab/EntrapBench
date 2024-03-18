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
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GenerateDatabase {

  private static final Pattern pattern = Pattern.compile("(\\w+)\\|(\\w+)\\|(\\w+)\\s*(.*)(GN=\\w+)?(.*)");
  private static final Pattern pattern2 = Pattern.compile("GN=([^ ]+)");

  public static void main(String[] args) {
    if (args.length != 7) {
      System.out.println("Usage: java -cp EntrapBench.jar entrapment.GenerateDatabase <UniProt fasta file path> <cut sites> <protect sites> <cleavage from C-term: 0=false, 1 = true> <number of entrapment proteins for each target protein> <entrapment prefix> <add prefix>");
      System.exit(1);
    }

    Path fastaPath = Paths.get(args[0]).toAbsolutePath();
    String cutSites = args[1];
    String protectSites = args[2];
    boolean cleavageFromCTerm = args[3].contentEquals("1");
    int N = 10;
    String prefix = args[5];
    boolean addPrefix = args[6].contentEquals("1");

    try {
      N = Integer.parseInt(args[4]);
    } catch (Exception e) {
      System.out.println("The last argument must be an integer.");
      System.exit(1);
    }

    if (!Files.exists(fastaPath) || !Files.isReadable(fastaPath) || !Files.isRegularFile(fastaPath)) {
      System.out.println("The fasta file " + args[0] + " is not valid.");
      System.exit(1);
    }

    String outputFile = fastaPath.getParent().resolve("target_shuffle_" + fastaPath.getFileName()).toAbsolutePath().toString();

    try {
      String line;
      BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
      BufferedReader reader = new BufferedReader(new FileReader(fastaPath.toFile()));
      String header = null;
      StringBuilder sequence = new StringBuilder();
      while ((line = reader.readLine()) != null) {
        line = line.trim();
        if (line.isEmpty()) {
          continue;
        }

        if (line.startsWith(">")) {
          if (sequence.length() > 0) {
            writeProtein(writer, header, sequence, cutSites, protectSites, cleavageFromCTerm, N, prefix, addPrefix);
          }
          sequence = new StringBuilder();
          header = line.substring(1);
        } else {
          sequence.append(line);
        }
      }

      if (sequence.length() > 0) {
        writeProtein(writer, header, sequence, cutSites, protectSites, cleavageFromCTerm, N, prefix, addPrefix);
      }

      reader.close();
      writer.close();
    } catch (Exception ex) {
      ex.printStackTrace();
      System.exit(1);
    }

  }

  private static void writeProtein(BufferedWriter writer, String header, StringBuilder sequence, String cutSites, String protectSites, boolean cleavageFromCTerm, int N, String prefix, boolean addPrefix) throws Exception {
    writer.write(">" + header + "\n");
    writer.write(sequence + "\n");

    String part1;
    String part2;
    String part3 = null;
    String part4 = null;
    Matcher matcher = pattern.matcher(header);
    if (matcher.matches()) {
      part1 = matcher.group(1).trim();
      part2 = matcher.group(2).trim();
      part3 = matcher.group(3).trim();
      part4 = matcher.group(4).trim();
    } else {
      part1 = "sp";
      String[] parts = header.split("\\s");
      part2 = parts[0];
      if (parts.length > 1) {
        part3 = parts[1];
      }
      if (parts.length > 2) {
        part4 = parts[2];
      }
    }

    String[] shuffledProteins = shuffleSeqFY(sequence.toString(), cutSites, protectSites, cleavageFromCTerm, N);
    for (int i = 0; i < shuffledProteins.length; ++i) {
      writer.write(">" + (addPrefix ? appendPrefix(prefix, i, part1) : part1) + "|" + appendPrefix(prefix, i, part2) + (part3 == null ? "" : "|" + appendPrefix(prefix, i, part3)) + (part4 == null ? "" : " " + replaceGN(prefix, i, part4)) + "\n");
      writer.write(shuffledProteins[i] + "\n");
    }
  }

  public static String[] shuffleSeqFY(String sequence, String cleavageSite, String protectionSite, boolean cleavageFromCTerm, int N) {
    // todo: A protection site may be shuffled, which may result in "non-existing" peptides after digestion.
    // todo: A "potential" protection site may also be shuffled to the side of a cleavage site so that it prevents a peptide from being digested.
    String sequenceToBeShuffled;
    if (sequence.startsWith("M")) {
      sequenceToBeShuffled = sequence.substring(1);
    } else {
      sequenceToBeShuffled = sequence;
    }

    Pattern digestSitePattern = getDigestSitePattern(cleavageSite, protectionSite, cleavageFromCTerm);
    Set<Integer> cutSiteSet = new HashSet<>();
    Matcher matcher = digestSitePattern.matcher(sequenceToBeShuffled);
    while (matcher.find()) {
      cutSiteSet.add(matcher.start());
    }

    int time = 0;
    Integer[] cutSiteArray = cutSiteSet.toArray(new Integer[0]);
    Arrays.sort(cutSiteArray);

    StringBuilder[] shuffledProteins = new StringBuilder[N];
    for (int i = 0; i < N; ++i) {
      shuffledProteins[i] = new StringBuilder();
      if (sequence.startsWith("M")) {
        shuffledProteins[i].append("M");
      }
    }

    int startIdx;
    int endIdx;
    Set<String> generatedShuffles = new HashSet<>();
    for (int k = 0; k <= cutSiteArray.length; ++k) {
      startIdx = k == 0 ? 0 : cutSiteArray[k - 1] + 1;
      endIdx = k == cutSiteArray.length ? sequenceToBeShuffled.length() : cutSiteArray[k];
      if (endIdx - startIdx > 2) {
        String targetSequence = sequenceToBeShuffled.substring(startIdx, endIdx);
        char[] targetArray = targetSequence.toCharArray();
        char[] shuffleArray = new char[targetArray.length];
        System.arraycopy(targetArray, 0, shuffleArray, 0, targetArray.length);
        generatedShuffles.clear();
        for (int l = 0; l < N; ++l) {
          Random random = new Random(l);
          String shuffledSequence;
          do {
            for (int i = 0; i < shuffleArray.length; ++i) {
              int j = random.nextInt(shuffleArray.length);
              while (j == i) {
                j = random.nextInt(shuffleArray.length);
              }
              char temp = shuffleArray[i];
              shuffleArray[i] = shuffleArray[j];
              shuffleArray[j] = temp;
            }
            shuffledSequence = String.valueOf(shuffleArray);
            ++time;
          } while (time < 10 && (targetSequence.contentEquals(shuffledSequence) || generatedShuffles.contains(shuffledSequence)));
          generatedShuffles.add(shuffledSequence);
          shuffledProteins[l].append(shuffledSequence).append(endIdx == sequenceToBeShuffled.length() ? "" : sequenceToBeShuffled.charAt(endIdx));
        }
      } else {
        for (int l = 0; l < N; ++l) {
          shuffledProteins[l].append(sequenceToBeShuffled, startIdx, endIdx).append(endIdx == sequenceToBeShuffled.length() ? "" : sequenceToBeShuffled.charAt(endIdx));
        }
      }
    }

    String[] output = new String[N];
    for (int i = 0; i < N; ++i) {
      output[i] = shuffledProteins[i].toString();
    }

    return output;
  }

  private static String appendPrefix(String prefix, int i, String s) {
    return prefix + "_" + i + "_" + s;
  }

  private static String replaceGN(String prefix, int i, String s) {
    Matcher matcher = pattern2.matcher(s);
    if (matcher.find()) {
      return matcher.replaceFirst("GN=" + prefix + "_" + i + "_" + matcher.group(1));
    } else {
      return s;
    }
  }

  private static Pattern getDigestSitePattern(String cleavageSite, String protectionSite, boolean cleavageFromCTerm) {
    Pattern digestSitePattern;
    if (cleavageFromCTerm) {
      if (protectionSite.contentEquals("-")) {
        digestSitePattern = Pattern.compile("[" + cleavageSite + "]");
      } else {
        digestSitePattern = Pattern.compile("[" + cleavageSite + "](?![" + protectionSite + "])");
      }
    } else {
      if (protectionSite.contentEquals("-")) {
        digestSitePattern = Pattern.compile("[" + cleavageSite + "]");
      } else {
        digestSitePattern = Pattern.compile("(?<![" + protectionSite + "])" + "[" + cleavageSite + "]");
      }
    }
    return digestSitePattern;
  }

}
