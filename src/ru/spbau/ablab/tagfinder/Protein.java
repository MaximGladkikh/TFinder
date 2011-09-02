package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.*;
import ru.spbau.ablab.tagfinder.util.trie.AKAutomaton;

import java.util.Arrays;

import static ru.spbau.ablab.tagfinder.TagGenerator.AA_LET;
import static ru.spbau.ablab.tagfinder.TagGenerator.AA_MONO_MASS;

public class Protein {
    public static final double[] AA_MASS_ARRAY = new double[256];

    static {
        for (int i = 0; i < AA_LET.length; ++i) {
            AA_MASS_ARRAY[AA_LET[i]] = AA_MONO_MASS[i];
        }
    }

    private final String shortName;
    private final String name;
    private final String protein;
    private final String revProtein;
    private Integer hashCode;
    public double[] masses;
    public final double parentMass;
    private double bestAbs;
    private static final double MAX_ACCEPTED_DISTANCE = ConfigReader.getDoubleProperty("MAX_ACCEPTED_DISTANCE");

    public Protein(String s, String name) {
        protein = s;
        masses = getMasses(s);
        parentMass = getStringMass(s);// + Database.WATER_MASS;
        this.name = name;
        revProtein = new StringBuilder(s).reverse().toString();
        String subName = name.substring(name.indexOf("NP_"));
        shortName = subName.substring(0, subName.indexOf('|'));
    }

    public String getString() {
        return protein;
    }

    public String getName() {
        return name;
    }

    public String getShortName() {
        return shortName;
    }

    @Override
    public int hashCode() {
        return hashCode == null ? hashCode = (int) StringUtil.getHash(name) : hashCode;
    }

    private double[] getMasses(String s) {
//		double[] ans = new double[s.length() + 2];
//		for (int i = 0; i < s.length(); ++i) {
//			ans[i + 1] = ans[i] + AA_MASS_ARRAY[s.charAt(i)];
//		}
        double[] ans = new double[s.length() * 2 + 2];
        for (int i = 0; i < s.length(); ++i) {
            ans[i + 1] = ans[i] + AA_MASS_ARRAY[s.charAt(i)];
        }
        double sum = Database.WATER_MASS;
        ans[s.length() + 1] = sum;
        for (int i = 0; i < s.length(); ++i) {
            ans[s.length() + 2 + i] = sum += AA_MASS_ARRAY[s.charAt(s.length() - 1 - i)];
        }
        Arrays.sort(ans);
        return ans;
    }

    public static double getStringMass(String s) {
        double sum = 0;
        for (char c : s.toCharArray()) {
            sum += TagGenerator.AA_MONO_MASS[AKAutomaton.AA_INDEX[c]];
        }
        return sum;
    }

    public boolean contains(Path path) {
        boolean ans = getMaxMatch(path) == path.length();
        assert TagGenerator.EDGE_OF_TWO_AA || ans == (protein.contains(path.toString()) || revProtein.contains(path.toString()));
        return ans;
    }

    public int getMaxMatch(Path path) {
        bestAbs = Double.POSITIVE_INFINITY;
        return Math.max(getMaxMatch(path.getEdges(), path.beginMass, path.length()), getMaxMatch(path.getReversedEdges(), parentMass - path.beginMass - path.getMass() + Database.WATER_MASS, path.length()));
    }

    private int getMaxMatch(Edge[] edges, double beginMass, int length) {
        int maxScore = 0;
        double mass = 0;
        for (int i = 0; i <= protein.length(); ++i) {
            double abs = Math.abs(beginMass - mass);
            if (abs > MAX_ACCEPTED_DISTANCE) {
                continue;
            }
            int pos = i;
            int matchedEdges = 0;
            j:
            for (int j = 0; j < edges.length && pos < protein.length(); ++j) {
                int len = 0;
                for (char[] decoding : edges[j].getDecodings()) {
                    boolean ok = true;
                    len = decoding.length;
                    for (int k = 0; k < decoding.length && pos + decoding.length <= protein.length(); ++k) {
                        if (protein.charAt(pos + k) != decoding[k]) {
                            ok = false;
                            break;
                        }
                    }
                    if (ok) {
                        pos += len;
                        ++matchedEdges;
                        continue j;
                    }
                }
                pos += len;
            }
            maxScore = Math.max(matchedEdges, maxScore);
            if (maxScore == length) {
                bestAbs = Math.min(bestAbs, abs);
            }
            if (i < protein.length()) {
                mass += AA_MONO_MASS[AKAutomaton.AA_INDEX[protein.charAt(i)]];
            }
        }
        return maxScore;
    }

    public double getBestLastAlignment() {
        return bestAbs;
    }

    public int getMatchScore(Spectrum spectrum, double mult, double shift) {
        int ans = 0;
        for (Envelope envelope : spectrum.envelopes) {
            double needMass = envelope.getMass() * mult + shift;
            int index = ArrayUtil.getClosestIndex(masses, needMass);
            if (MassComparator.sameForDeletion(needMass, masses[index])) {
                ++ans;
            }
        }
        return ans;
    }
}
