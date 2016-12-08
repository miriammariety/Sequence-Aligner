/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pkg201340385anslab5;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author ty
 */
public class Aligner {

    String output = "";
    static ArrayList<String[]> alignments = new ArrayList<>();
    ArrayList<Sequence> sequences;
    static HashMap<String, Integer> freqOne;
    static HashMap<String, Integer> freqTwo;
    String seq1, seq2;
    int row, column, gap, match, mismatch;
    String type;
    Cell[][] matrix;
    ScoringMatrix blosum;
    ScoringMatrix pam;

    public Aligner(ArrayList<Sequence> sequences, int match, int mismatch, int gap) {
        this.sequences = sequences;
        this.seq1 = sequences.get(0).sequence;
        this.seq2 = sequences.get(1).sequence;
        this.row = seq1.length();
        this.column = seq2.length();
        this.match = match;
        this.mismatch = mismatch;
        this.gap = gap;

    }

    public Aligner(ArrayList<Sequence> sequences) {
        this.sequences = sequences;
        this.seq1 = sequences.get(0).sequence;
        this.seq2 = sequences.get(1).sequence;
        this.row = seq1.length();
        this.column = seq2.length();
    }
    
    public void runNucleotide() {
        int score;
        freqOne = initializeFrequencies(true);
        freqTwo = initializeFrequencies(true);
        matrix = initializeMatrix(seq1.length() + 1, seq2.length() + 1, true, gap);
        fillMatrixNucleotide(matrix, new int[]{match, mismatch, gap});
        backtrack(matrix[seq1.length()][seq2.length()], 1, 1, "", "");
        score = matrix[seq1.length()][seq2.length()].value;

        output += "NUCLEOTIDE PAIRWISE SEQUENCE ALIGNMENT\n";
        countFrequencies(sequences.get(0), sequences.get(1));
        printOutput(score);
        freqOne.clear();
        freqTwo.clear();
        alignments.clear();
        seq1 = "";
        seq2 = "";
        matrix = initializeMatrix(0, 0, true, 0);
    }

    public void runGlobalProtein() throws IOException {
        int score;
        freqOne = initializeFrequencies(false);
        freqTwo = initializeFrequencies(false);
        matrix = initializeMatrix(seq1.length() + 1, seq2.length() + 1, true, -8);
        globalFillMatrixProtein(matrix);
        backtrack(matrix[seq1.length()][seq2.length()], 1, 1, "", "");
        score = matrix[seq1.length()][seq2.length()].value;
        
        output += "Protein Global Pairwise Sequence Alignment\n";
        countFrequencies(sequences.get(0), sequences.get(1));
        printOutput(score);
        freqOne.clear();
        freqTwo.clear();
        alignments.clear();
        seq1 = "";
        seq2 = "";
        matrix = initializeMatrix(0, 0, true, 0);
    }
    
    public void runLocalProtein() throws IOException {
        int score;
        ArrayList<Cell> cellArray = new ArrayList<>();
        freqOne = initializeFrequencies(false);
        freqTwo = initializeFrequencies(false);
        matrix = initializeMatrix(seq1.length()+1, seq2.length()+1, false, -4);
        localFillMatrixProtein(matrix);
        score = findMaxScore(matrix);
        cellArray = findMaxCell(score);
        for (int i = 0; i < cellArray.size(); i++) {
            backtrack(cellArray.get(i), seq1.length()+1 - cellArray.get(i).j, seq2.length()+1 - cellArray.get(i).i, "", "");
        }
        
        output += "Protein Local Pairwise Sequence Alignment\n";
        countFrequencies(sequences.get(0), sequences.get(1));
        printOutput(score);
        freqOne.clear();
        freqTwo.clear();
        alignments.clear();
        seq1 = "";
        seq2 = "";
        matrix = initializeMatrix(0, 0, true, 0);
    }

    public int findMaxScore(Cell[][] matrix) {
        Cell currentCell;
        int max = -9999;

        for (int i = 1; i < seq2.length() + 1; i++) {
            for (int j = 1; j < seq1.length() + 1; j++) {
                currentCell = matrix[j][i];

                if (currentCell.value > max) {
                    max = currentCell.value;
                }
            }
        }

        return max;
    }

    public ArrayList<Cell> findMaxCell(int max) {
        ArrayList<Cell> cellArray = new ArrayList<Cell>();
        Cell currentCell;

        for (int i = 1; i < seq2.length() + 1; i++) {
            for (int j = 1; j < seq1.length() + 1; j++) {
                currentCell = matrix[j][i];

                if (currentCell.value == max) {
                    cellArray.add(currentCell);
                }
            }
        }

        return cellArray;
    }

    public void fillMatrixNucleotide(Cell[][] matrix, int[] scoring) {
        Cell currentCell;
        int match = scoring[0];
        int mismatch = scoring[1];
        int gap = scoring[2];

        for (int i = 1; i < seq2.length() + 1; i++) {
            for (int j = 1; j < seq1.length() + 1; j++) {
                currentCell = matrix[j][i];

                if (seq1.substring(j - 1, j).matches(
                        seq2.substring(i - 1, i))) {
                    currentCell.value = currentCell.max(matrix[j - 1][i - 1],
                            matrix[j - 1][i], matrix[j][i - 1], match, gap, -9999);
                } else {
                    currentCell.value = currentCell.max(matrix[j - 1][i - 1],
                            matrix[j - 1][i], matrix[j][i - 1], mismatch, gap, -9999);
                }

                matrix[j][i] = currentCell;
            }
        }
    }

    public void localFillMatrixProtein(Cell[][] matrix) throws IOException {
        Cell currentCell;
        blosum = new ScoringMatrix("blosum62.txt");

        for (int i = 1; i < seq2.length() + 1; i++) {
            for (int j = 1; j < seq1.length() + 1; j++) {
                currentCell = matrix[j][i];

                currentCell.value = currentCell.max(matrix[j - 1][i - 1],
                        matrix[j - 1][i], matrix[j][i - 1],
                        blosum.getScore(seq1.charAt(j-1), seq2.charAt(i-1)),
                        -4, 0);

                matrix[j][i] = currentCell;
            }
        }
    }

    public void globalFillMatrixProtein(Cell[][] matrix) throws IOException {
        Cell currentCell;
        pam = new ScoringMatrix("pam120.txt");

        for (int i = 1; i < seq2.length() + 1; i++) {
            for (int j = 1; j < seq1.length() + 1; j++) {
                currentCell = matrix[j][i];

                currentCell.value = currentCell.max(matrix[j - 1][i - 1],
                        matrix[j - 1][i], matrix[j][i - 1],
                        pam.getScore(seq1.charAt(j-1), seq2.charAt(i-1)),
                        -8, -9999);

                matrix[j][i] = currentCell;
            }
        }
    }

    public void backtrack(Cell cell, int s1, int s2, String seqOne, String seqTwo) {
        boolean test = false;

        if (cell.diagonal != null) {
            if (test == true) {
                seqOne = seqOne.substring(1);
                seqTwo = seqTwo.substring(1);
            }

            seqOne = seq1.substring(
                    seq1.length() - s1, seq1.length() - s1 + 1) + seqOne;
            seqTwo = seq2.substring(
                    seq2.length() - s2, seq2.length() - s2 + 1) + seqTwo;
            backtrack(cell.diagonal, s1 + 1, s2 + 1, seqOne, seqTwo);
            test = true;
        }

        if (cell.left != null) {
            if (test == true) {
                seqOne = seqOne.substring(1);
                seqTwo = seqTwo.substring(1);
            }

            seqOne = seq1.substring(
                    seq1.length() - s1, seq1.length() - s1 + 1) + seqOne;
            seqTwo = "-" + seqTwo;
            backtrack(cell.left, s1 + 1, s2, seqOne, seqTwo);
            test = true;
        }

        if (cell.top != null) {
            if (test == true) {
                seqOne = seqOne.substring(1);
                seqTwo = seqTwo.substring(1);
            }

            seqOne = "-" + seqOne;
            seqTwo = seq2.substring(
                    seq2.length() - s2, seq2.length() - s2 + 1) + seqTwo;
            backtrack(cell.top, s1, s2 + 1, seqOne, seqTwo);
            test = true;
        }

        if (cell.diagonal == null && cell.left == null && cell.top == null) {
            alignments.add(new String[]{seqOne, seqTwo});
        } else {
            return;
        }
    }

    public void countFrequencies(Sequence sequence1, Sequence sequence2) {
        String key;

        for (int i = 0; i < sequence1.sequence.length(); i++) {
            key = sequence1.sequence.substring(i, i + 1);

            if (freqOne.containsKey(key)) {
                freqOne.put(key, freqOne.get(key) + 1);
            } else {
                freqOne.put(key, 1);
            }

        }

        for (int i = 0; i < sequence2.sequence.length(); i++) {
            key = sequence2.sequence.substring(i, i + 1);

            if (freqTwo.containsKey(key)) {
                freqTwo.put(key, freqTwo.get(key) + 1);
            } else {
                freqTwo.put(key, 1);
            }
        }
    }

    public void printFrequencies() {
        String format = "%-4s";
        String keys = "";
        String frequencies = "";
        Iterator f1 = freqOne.entrySet().iterator();
        Iterator f2 = freqTwo.entrySet().iterator();

        while (f1.hasNext()) {
            Map.Entry pair = (Map.Entry) f1.next();
            keys = keys + pair.getKey();
            frequencies = frequencies + pair.getValue();
            f1.remove();
        }

        keys = keys.replaceAll(".(?=.)", String.format(format, "$0"));
        frequencies = frequencies.replaceAll(".(?=.)", String.format(format, "$0"));
        output += "     " + keys + "\tTOTAL" + "\n";
        output += "     " + frequencies + "\t" + seq1.length() + "\n\n";
        keys = "";
        frequencies = "";

        while (f2.hasNext()) {
            Map.Entry pair = (Map.Entry) f2.next();
            keys = keys + pair.getKey();
            frequencies = frequencies + pair.getValue();
            f2.remove();
        }

        keys = keys.replaceAll(".(?=.)", String.format(format, "$0"));
        frequencies = frequencies.replaceAll(".(?=.)", String.format(format, "$0"));
        output += "     " + keys + "\tTOTAL" + "\n";
        output += "     " + frequencies + "\t" + seq2.length() + "\n";
    }

    public void printOutput(int score) {
        output += "Pairwise Sequence Alignment ver. 1.0 by Miriam Ty\n"; 
        output += "Date Created: " + new Date() + "\n\n";
        output += "Sequences: \n";
        output += "     " + sequences.get(0).comment + "\n" + "     " + formatString(seq1) + "\n";
        output += "     " + sequences.get(1).comment + "\n" + "     " + formatString(seq2) + "\n\n";
        output += "Lengths:\n" + 
                "     First Seqeunce: " + seq1.length() + "\n" + 
                "     Second Sequence: "+ seq2.length() + "\n\n";
        output += "Frequency Occurence:\n";
        printFrequencies();
        output += "\n";
        output += ("Optimal alignment result(s):\n\n");
        printAlignment(sequences.get(0), sequences.get(1));
        output += ("Score: " + score) + "\n";
    }

    public void printMatrix(Cell[][] matrix) {
        String format = "%-5s";

        for (int i = 0; i < seq2.length() + 1; i++) {
            for (int j = 0; j < seq2.length() + 1; j++) {
                System.out.printf(format, matrix[j][i].value);
            }
            System.out.println("");
        }
    }

    public void printAlignment(Sequence s1, Sequence s2) {
        String matches = "";
        int k = 0;

        for (int i = 0; i < alignments.size(); i++) {
            for (int j = 0; j < alignments.get(i)[0].length(); j++, k++) {
                if (k == 10) {
                    matches += " ";
                    k = 0;
                }

                if (alignments.get(i)[0].substring(j, j + 1).equals("-") || alignments.get(i)[1].substring(j, j + 1).equals("-")) {
                    matches += " ";
                } else if (alignments.get(i)[0].substring(j, j + 1).matches(alignments.get(i)[1].substring(j, j + 1))) {
                    matches += "*";
                } else {
                    matches += ".";
                }
            }

            output += formatString(alignments.get(i)[0]) + "\n";
            output += (matches) + "\n";
            output += (formatString(alignments.get(i)[1])) + "\n";
            output += "\n";
            matches = "";
            k = 0;
        }

    }

    public String formatString(String alignment) {
        String formattedAlignment = "";
        int k = 0;
        for (int i = 0; i < alignment.length(); i++, k++) {
            if (k == 10) {
                formattedAlignment += " ";
                k = 0;
            }
            formattedAlignment += alignment.substring(i, i + 1);
        }

        return formattedAlignment;
    }

    public Cell[][] initializeMatrix(int m, int n, boolean isGlobal, int gap) {
        Cell[][] matrix = new Cell[m][n];

        if (isGlobal) {
            for (int i = 0; i < n; i++) {
                matrix[0][i] = new Cell(i * gap, i, 0);
            }

            for (int i = 0; i < m; i++) {
                matrix[i][0] = new Cell(i * gap, 0, i);
            }

            for (int i = 1; i < n; i++) {
                for (int j = 1; j < m; j++) {
                    matrix[j][i] = new Cell(0, i, j);
                }
            }
        } else {
            for (int i = 0; i < n; i++) {
                matrix[0][i] = new Cell(0, i, 0);
            }

            for (int i = 0; i < m; i++) {
                matrix[i][0] = new Cell(0, 0, i);
            }

            for (int i = 1; i < n; i++) {
                for (int j = 1; j < m; j++) {
                    matrix[j][i] = new Cell(0, i, j);
                }
            }
        }

        return matrix;
    }

    public static HashMap<String, Integer> initializeFrequencies(boolean isNucleotide) {
        HashMap<String, Integer> map = new HashMap<>();

        if (isNucleotide) {
            map.put("A", 0);
            map.put("C", 0);
            map.put("G", 0);
            map.put("T", 0);
        } else {
            map.put("A", 0);
            map.put("R", 0);
            map.put("N", 0);
            map.put("D", 0);
            map.put("C", 0);
            map.put("Q", 0);
            map.put("E", 0);
            map.put("G", 0);
            map.put("H", 0);
            map.put("I", 0);
            map.put("L", 0);
            map.put("K", 0);
            map.put("M", 0);
            map.put("F", 0);
            map.put("P", 0);
            map.put("S", 0);
            map.put("T", 0);
            map.put("W", 0);
            map.put("Y", 0);
            map.put("V", 0);
            map.put("B", 0);
            map.put("Z", 0);
            map.put("X", 0);
        }

        return map;
    }

}
