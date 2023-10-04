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

import static org.junit.Assert.*;

import org.junit.Test;

public class GenerateDatabaseTest {

  @Test
  public void shuffleSeqFY() {
    String originalSequence = "AICQFFLQGR";
    String[] output = GenerateDatabase.shuffleSeqFY(originalSequence, "KR", "P", true, 3);
    System.out.println("Original: " + originalSequence);
    for (String s : output) {
      System.out.println("   Decoy: " + s);
    }
    System.out.println();

    originalSequence = "MDPLFQQTHKAICQFFLQGR";
    output = GenerateDatabase.shuffleSeqFY(originalSequence, "KR", "P", true, 5);
    System.out.println("Original: " + originalSequence);
    for (String s : output) {
      System.out.println("   Decoy: " + s);
    }
    System.out.println();

    originalSequence = "MDPLFQQTHKPAICQFFLQGR";
    output = GenerateDatabase.shuffleSeqFY(originalSequence, "KR", "P", true, 6);
    System.out.println("Original: " + originalSequence);
    for (String s : output) {
      System.out.println("   Decoy: " + s);
    }
    System.out.println();

    originalSequence = "MEYMAESTDRAADFQLHTHVNDGTEFGGSIYQKAAFVAYALAFPRAALEEANGEIEKAAMEALVVEVT";
    output = GenerateDatabase.shuffleSeqFY(originalSequence, "KR", "P", true, 10);
    System.out.println("Original: " + originalSequence);
    for (String s : output) {
      System.out.println("   Decoy: " + s);
    }
    System.out.println();
  }

  @Test
  public void main() {
    GenerateDatabase.main(new String[]{"G:\\dev\\2021-03-16-reviewed-contam-UP000002311.fas", "KR", "P", "1", "3"});
  }
}