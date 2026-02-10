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
    if (args.length != 6) {
      System.out.println("Usage: java -cp EntrapBench.jar entrapment.GenerateDatabase <UniProt fasta file path> <cut sites> <protect sites> <cleavage from C-term: 0=false, 1 = true> <number of entrapment proteins for each target protein> <entrapment style>");
      System.out.println("entrapment style: 0 = add \"entrapment_\" prefix to the protein ID, 1 = add \"_p_target\" suffix to the protein ID which is used by https://doi.org/10.1038/s41592-025-02719-x");
      System.exit(1);
    }

    Path fastaPath = Paths.get(args[0]).toAbsolutePath();
    String cutSites = args[1];
    String protectSites = args[2];
    boolean cleavageFromCTerm = args[3].contentEquals("1");
    int N = Integer.parseInt(args[4]);
    int entrapmentStyle = Integer.parseInt(args[5]);

    if (entrapmentStyle == 1 && N != 1) {
      N = 1;
      System.out.println("The number of entrapment proteins for each target protein is set to 1.");
    }

    if (!Files.exists(fastaPath) || !Files.isReadable(fastaPath) || !Files.isRegularFile(fastaPath)) {
      System.out.println("The fasta file " + args[0] + " is not valid.");
      System.exit(1);
    }
    String outputFile1 = fastaPath.getParent().resolve("target_shuffle_" + fastaPath.getFileName()).toAbsolutePath().toString();
    String outputFile2 = fastaPath.getParent().resolve("target_shuffle_pep_" + fastaPath.getFileName()).toAbsolutePath().toString();

    if (Files.exists(Paths.get(outputFile1))) {
      System.out.println("The output file " + outputFile1 + " already exists.");
      System.exit(1);
    }

    if (Files.exists(Paths.get(outputFile2))) {
      System.out.println("The output file " + outputFile2 + " already exists.");
      System.exit(1);
    }

    try {
      String line;
      BufferedWriter writer1 = new BufferedWriter(new FileWriter(outputFile1));
      BufferedWriter writer2 = new BufferedWriter(new FileWriter(outputFile2));

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
            writeProtein(writer1, writer2, header, sequence, cutSites, protectSites, cleavageFromCTerm, N, entrapmentStyle);
          }
          sequence = new StringBuilder();
          header = line.substring(1);
        } else {
          sequence.append(line);
        }
      }

      if (sequence.length() > 0) {
        writeProtein(writer1, writer2, header, sequence, cutSites, protectSites, cleavageFromCTerm, N, entrapmentStyle);
      }

      reader.close();
      writer1.close();
      writer2.close();
    } catch (Exception ex) {
      ex.printStackTrace();
      System.exit(1);
    }

  }

  public static void writePeptide(BufferedWriter writer2, String tt, String ee, String cleavageSite, String protectionSite, boolean cleavageFromCTerm) throws Exception {
    Pattern digestSitePattern = getDigestSitePattern(cleavageSite, protectionSite, cleavageFromCTerm);
    Set<Integer> cutSiteSet = new HashSet<>();
    cutSiteSet.add(0);
    cutSiteSet.add(tt.length());
    Matcher matcher = digestSitePattern.matcher(tt);
    while (matcher.find()) {
      cutSiteSet.add(matcher.start() + 1);
    }
    Integer[] cutSiteArray = cutSiteSet.toArray(new Integer[0]);
    Arrays.sort(cutSiteArray);

    for (int i = 0; i < cutSiteArray.length - 1; ++i) {
      String t = tt.substring(cutSiteArray[i], Math.min(cutSiteArray[i + 1], tt.length()));
      String e = ee.substring(cutSiteArray[i], Math.min(cutSiteArray[i + 1], ee.length()));
      // default peptide length range used by https://github.com/Noble-Lab/FDRBench
      if (t.length() >= 7 && t.length() <= 35) {
        writer2.write(">sp|" + t + "_target|" + t + "_target\n");
        writer2.write(t + "\n");
        writer2.write(">sp|" + e + "_p_target|" + e + "_p_target\n");
        writer2.write(e + "\n");
      }
      if (i + 2 < cutSiteArray.length) {
        t = tt.substring(cutSiteArray[i], Math.min(cutSiteArray[i + 2], tt.length()));
        e = ee.substring(cutSiteArray[i], Math.min(cutSiteArray[i + 2], ee.length()));
        // default peptide length range used by https://github.com/Noble-Lab/FDRBench
        if (t.length() >= 7 && t.length() <= 35) {
          writer2.write(">sp|" + t + "_target|" + t + "_target\n");
          writer2.write(t + "\n");
          writer2.write(">sp|" + e + "_p_target|" + e + "_p_target\n");
          writer2.write(e + "\n");
        }
      }
    }
  }

  private static void writeProtein(BufferedWriter writer1, BufferedWriter writer2, String header, StringBuilder sequence, String cutSites, String protectSites, boolean cleavageFromCTerm, int N, int entrapmentStyle) throws Exception {
    String sequence2 = sequence.toString().replaceAll("I", "L");

    writer1.write(">" + header + "\n");
    writer1.write(sequence2 + "\n");

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

    String[] shuffledProteins = shuffleSeqFY(sequence2, cutSites, protectSites, cleavageFromCTerm, N);
    for (int i = 0; i < shuffledProteins.length; ++i) {
      writer1.write(">" + appendEntrapmentMarker(entrapmentStyle, i, part1) + "|" + appendEntrapmentMarker(entrapmentStyle, i, part2) + (part3 == null ? "" : "|" + appendEntrapmentMarker(entrapmentStyle, i, part3)) + (part4 == null ? "" : " " + replaceGN(entrapmentStyle, i, part4)) + "\n");
      writer1.write(shuffledProteins[i] + "\n");
    }

    if (entrapmentStyle == 1) {
      writePeptide(writer2, sequence2, shuffledProteins[0], cutSites, protectSites, cleavageFromCTerm);
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

  private static String appendEntrapmentMarker(int entrapmentStyle, int i, String s) {
    return entrapmentStyle == 0 ? "entrapment_" + i + "_" + s : s + "_p_target";
  }

  private static String replaceGN(int entrapmentStyle, int i, String s) {
    Matcher matcher = pattern2.matcher(s);
    if (matcher.find()) {
      return entrapmentStyle == 0 ?
          matcher.replaceFirst("GN=entrapment_" + i + "_" + matcher.group(1)) :
          matcher.replaceFirst("GN=" + matcher.group(1) + "_p_target");
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
