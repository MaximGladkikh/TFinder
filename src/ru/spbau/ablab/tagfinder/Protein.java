package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.*;

import java.util.ArrayList;
import java.util.Collections;

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
        protein = s = clearFromTruncated(s);
//        System.err.println(s);
        if (s.indexOf('I') >= 0) {
            throw new AssertionError();
        }
        parentMass = getStringMass(s) + Database.WATER_MASS;
        masses = getMasses(s);
        this.name = name;
        revProtein = new StringBuilder(s).reverse().toString();
//        System.err.println("name=" + name + " , index=" + name.indexOf("NP_"));
        String subName = name.substring(name.indexOf("NP_"));
        shortName = subName.substring(0, subName.indexOf('|'));
    }

    private String clearFromTruncated(String s) {
        int l = 0;
        int r = s.length();
        for (int i = 0; i + 3 < s.length(); ++i) {
           char c = s.charAt(i);
           if (c == '(') {
               break;
           }
           if (c == '.') {
               l = i + 1;
           }
        }
        for (int i = s.length() - 1; i > l && i > 3; --i) {
            char c = s.charAt(i);
            if (c == ']') {
                break;
            }
            if (c == '.') {
                r = i;
            }
        }
        return s.substring(l, r);
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
//        double[] ans = new double[s.length() * 2 + 2];
        ArrayList<Double> list = new ArrayList<Double>();
        list.add(0.);
        list.add(parentMass);
        double sum = 0;
        for (int i = 0; i < s.length(); ++i) {
            char c = s.charAt(i);
            if (Character.isLetter(c)) {
                sum += TagGenerator.AA_MONO_MASS[AKAutomaton.AA_INDEX[c]];
            } else if (c == '[') {
//                System.err.println(string.substring(i));
                StringBuilder builder = new StringBuilder();
                for (++i; s.charAt(i) != ']'; ++i) {
                    builder.append(s.charAt(i));
                }
                sum += Double.parseDouble(builder.toString());
            }
            list.add(sum);
            list.add(parentMass - sum);
        }
        Collections.sort(list);
        double [] ans = new double[list.size()];
        for (int i = 0; i < ans.length; ++i) {
            ans[i] = list.get(i);
        }
//        double sum = Database.WATER_MASS;
//        ans[s.length() + 1] = sum;
//        for (int i = 0; i < s.length(); ++i) {
//            ans[s.length() + 2 + i] = sum += AA_MASS_ARRAY[s.charAt(s.length() - 1 - i)];
//        }
//        Arrays.sort(ans);
        return ans;
    }

    public static double getStringMass(String string) {
        double sum = 0;
        char[] s = string.toCharArray();
        for (int i = 0; i < s.length; ++i) {
            char c = s[i];
            if (Character.isLetter(c)) {
                sum += TagGenerator.AA_MONO_MASS[AKAutomaton.AA_INDEX[c]];
            } else if (c == '[') {
//                System.err.println(string.substring(i));
                StringBuilder builder = new StringBuilder();
                for (++i; s[i] != ']'; ++i) {
                    builder.append(s[i]);
                }
                sum += Double.parseDouble(builder.toString());
            }
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
        return Math.max(getMaxMatch(path.getEdges(), path.beginMass, path.length()), getMaxMatch(path.getReversedEdges(), MassUtil.convertIonsType(path.beginMass + path.getMass(), parentMass), path.length()));
    }

    private int getMaxMatch(Edge[] edges, double beginMass, int length) {
        int maxScore = 0;
        double mass = 0;
        for (int i = 0; i <= protein.length(); ++i) {
            if (i < protein.length() && protein.charAt(i) == '[') {
                StringBuilder builder = new StringBuilder();
                for (++i; protein.charAt(i) != ']'; ++i) {
                    builder.append(protein.charAt(i));
                }
                mass += Double.parseDouble(builder.toString());
                continue;
            }
            double diff = beginMass - mass;
            if (diff > MAX_ACCEPTED_DISTANCE) {
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
                    int addDueToNonAInString = 0;
                    for (int k = 0; k < decoding.length && pos + decoding.length + addDueToNonAInString <= protein.length(); ++k) {
                        while (pos + k + addDueToNonAInString < protein.length() && !Character.isLetter(protein.charAt(pos + k + addDueToNonAInString))) {
                            ++addDueToNonAInString;
                        }
                        if (pos + k + addDueToNonAInString >= protein.length()) {
                            ok = false;
                            break;
                        }
                        if (protein.charAt(pos + k + addDueToNonAInString) != decoding[k]) {
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
            if (matchedEdges == length && Math.abs(bestAbs) > Math.abs(diff)) {
                bestAbs = diff;
//                bestAbs = Math.min(bestAbs, abs);
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
            if (MassUtil.sameForDeletion(needMass, masses[index])) {
                ++ans;
            }
        }
        return ans;
    }
}
