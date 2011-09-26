package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.database.Database;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.MassUtil;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;

import static ru.spbau.ablab.tagfinder.database.Database.*;

public class KDCalculator implements Runnable {
    private static TreeMap<Result, ArrayList<Integer>> bestPairs = new TreeMap<>();

    private Spectrum spectrum;
    private String sequence;
    private Protein protein;
    private Envelope[] envelopes;
    private String[] longestPath;
    private double[] beginMasses;
    private DisjointSetUnion union;

    public static void main(String[] args) throws FileNotFoundException {
        Database database = Database.getInstance();
        for (int id : database.filter(true)) {
            assert !Database.USE_VIRTUAL_SPECTRA;
            new KDCalculator(database.getSpectraDb().getSpectrum(id), database.getProteinPredictedByAlign(id)).run();
        }
        PrintWriter writer = new PrintWriter("kd.txt");
        for (Result pair : bestPairs.keySet()) {
            writer.println("(" + pair.longestPath.length() + ", " + pair.longestCorrectPath.length() + ") " + bestPairs.get(pair).size() + " times : " + pair.longestPath + " " + pair.longestCorrectPath + " " + bestPairs.get(pair));
        }
        writer.close();
    }

    public KDCalculator(Spectrum spectrum, Protein protein) {
        this.spectrum = spectrum;
        this.sequence = protein.getString();
        this.protein = protein;
        envelopes = spectrum.getEnvelopes();
        union = new DisjointSetUnion(envelopes.length);
    }

    @Override
    public void run() {
        computeKs();
        Result best = getResult();
        if (best == null) {
            return;
        }
        ArrayList<Integer> list = bestPairs.get(best);
        if (list == null) {
            bestPairs.put(best, list = new ArrayList<>());
        }
        list.add(spectrum.id);
    }

    private void computeKs() {
        String[] maxPath = new String[envelopes.length];
        longestPath = new String[envelopes.length];
        beginMasses = new double[envelopes.length];
        Arrays.fill(maxPath, "");
        Arrays.fill(longestPath, "");
        Arrays.fill(beginMasses, Double.NaN);

        for (int i = envelopes.length - 1; i >= 0; --i) {
            String best = "";
            for (Edge[] edges : TagGenerator.EDGES) {
                for (Edge edge : edges) {
                    double currentMass = envelopes[i].getMass();
                    double needMass = envelopes[i].getMass() + edge.getMass();
                    for (int next = spectrum.getFirstMatchingEnvelopeIndex(currentMass, needMass, edge.getMass()); next < spectrum.envelopes.length && MassUtil.edgeMatches(currentMass, spectrum.envelopes[next].getMass(), edge.getMass()); ++next) {
                        union.unite(i, next);
                        String s = "" + edge + maxPath[next];
                        if (best.length() < s.length()) {
                            best = s;
                        }
                    }
                }
            }
            maxPath[i] = best;
            if (best.indexOf('I') >= 0) {
                throw new AssertionError();
            }
        }
        for (int i = 0; i < envelopes.length; ++i) {
            int setId = union.findSet(i);
            if (longestPath[setId].length() < maxPath[i].length()) {
                longestPath[setId] = maxPath[i];
                beginMasses[setId] = envelopes[i].getMass();
            }
        }
    }

    public Result getResult() {
        Result bestPair = null;
        bestPair = getBestPair(bestPair, sequence, false);
        bestPair = getBestPair(bestPair, new StringBuilder(sequence).reverse().toString(), true);
//        System.out.println(spectrum.id + " (" + bestPair.longestPath.length() + ", " + bestPair.longestCorrectPath.length() + ") " + String.format("%.2f", bestPair.beginMass) + " " + bestPair.longestPath + " " + String.format("%.2f", bestPair.beginCorrectMass) + " " + bestPair.longestCorrectPath + " " + protein.getId() + " " + protein.getName() + " " + sequence);//        System.out.println(spectrum.id + " (" + bestPair.longestPath.length() + ", " + bestPair.longestCorrectPath.length() + ")");
        System.out.println(spectrum.id + " (" + bestPair.longestPath.length() + ", " + bestPair.longestCorrectPath.length() + ") ");// + String.format("%.2f", bestPair.beginMass) + " " + bestPair.longestPath + " " + String.format("%.2f", bestPair.beginCorrectMass) + " " + bestPair.longestCorrectPath + " " + protein.getId() + " " + protein.getName() + " " + sequence);//        System.out.println(spectrum.id + " (" + bestPair.longestPath.length() + ", " + bestPair.longestCorrectPath.length() + ")");
        return bestPair;
    }

    private int[][] dp;

    private int getD(int aaIndex, int peakIndex, String sequence) {
        if (aaIndex == sequence.length() || peakIndex == envelopes.length) {
            return 0;
        }
        if (dp[aaIndex][peakIndex] >= 0) {
            return dp[aaIndex][peakIndex];
        }
        int ans = 0;
        double currentMass = envelopes[peakIndex].getMass();
        double edgeMass = AA_MASS_ARRAY[sequence.charAt(aaIndex)];
        double needMass = currentMass + edgeMass;
        for (int next = spectrum.getFirstMatchingEnvelopeIndex(currentMass, needMass, edgeMass); next < spectrum.envelopes.length && MassUtil.edgeMatches(currentMass, spectrum.envelopes[next].getMass(), edgeMass); ++next) {
            ans = Math.max(ans, getD(aaIndex + 1, next, sequence) + 1);
        }
        return dp[aaIndex][peakIndex] = ans;
    }

    private Result getBestPair(Result bestPair, String sequence, boolean reversed) {
        reversed = false;
        dp = new int[sequence.length()][envelopes.length];
        for (int [] a : dp) {
            Arrays.fill(a, -1);
        }
        for (int i = 0; i < sequence.length(); ++i) {
            for (int j = 0; j < envelopes.length; ++j) {
                int matched = getD(i, j, sequence);
                int set = union.findSet(j);
                String k = longestPath[set];
                String matchedTag = sequence.substring(i, i + matched);
                if (reversed) {
                    matchedTag = new StringBuilder(matchedTag).reverse().toString();
                }
                if (matched > k.length()) {
                    throw new AssertionError(spectrum.id + " " + matched + " " + k + " : " + matchedTag);
                }
                Result pair = new Result(k, matchedTag, beginMasses[set], envelopes[j].getMass() - Protein.getStringMass(matchedTag));
                if (bestPair == null || pair.compareTo(bestPair) < 0) {
                    bestPair = pair;
                }
            }
        }
        return bestPair;
    }

    private static class Result implements Comparable<Result> {
        private String longestPath;
        private String longestCorrectPath;
        private double beginMass;
        private double beginCorrectMass;

        private Result(String longestPath, String longestCorrectPath, double beginMass, double beginCorrectMass) {
            this.longestPath = longestPath;
            this.longestCorrectPath = longestCorrectPath;
            this.beginMass = beginMass;
            this.beginCorrectMass = beginCorrectMass;
        }

        @Override
        public int compareTo(Result o) {
            int d = longestPath.length() - o.longestPath.length();
            if (d != 0) {
                return -d;
            }
            d = longestCorrectPath.length() - o.longestCorrectPath.length();
            if (d != 0) {
                return -d;
            }
            return (longestPath + " " + longestCorrectPath).compareTo(o.longestPath + " " + longestCorrectPath);
        }
    }

    private static class DisjointSetUnion {
        private int[] parent;

        public int findSet(int v) {
            if (parent[v] == v) {
                return v;
            }
            return parent[v] = findSet(parent[v]);
        }

        public void unite(int u, int v) {
            parent[findSet(u)] = findSet(v);
        }

        private DisjointSetUnion(int n) {
            parent = new int[n];
            for (int i = 0; i < n; ++i) {
                parent[i] = i;
            }
        }
    }
}
