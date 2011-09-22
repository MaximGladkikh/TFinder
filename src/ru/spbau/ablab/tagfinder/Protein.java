package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.util.MassUtil;
import ru.spbau.ablab.tagfinder.util.StringUtil;

import java.util.ArrayList;
import java.util.Collections;

import static ru.spbau.ablab.tagfinder.database.Database.AA_MASS_ARRAY;
import static ru.spbau.ablab.tagfinder.database.Database.WATER_MASS;

public class Protein {
    private final String name;
    private final String fullname;
    private final String protein;
    private final String revProtein;
    private Integer hashCode;
    public double[] masses;
    public final double parentMass;
    private double bestShift;
    private int id;

    public Protein(String s, String name, int id) {
        protein = s = clearFromTruncated(s);
        if (s.indexOf('I') >= 0/* || s.indexOf('.') >= '.'*/) {
            throw new AssertionError(s);
        }
        parentMass = getStringMass(s) + WATER_MASS;
        masses = getMasses(s);
        this.fullname = name;
        revProtein = new StringBuilder(s).reverse().toString();
        this.name = getShortName(name);
        this.id = id;
    }

    private String clearFromTruncated(String s) {
        int l = s.length();
        int r = -1;
        for (int i = 0; i < s.length(); ++i) {
            char c = s.charAt(i);
            if (c == '.') {
                l = i + 1;
                break;
            }
        }
        for (int i = s.length() - 1; i >= 0; --i) {
            char c = s.charAt(i);
            if (c == '.') {
                r = i;
                break;
            }
        }
        return l < r ? s.substring(l, r) : s;
    }

    public String getString() {
        return protein;
    }

    public String getFullname() {
        return fullname;
    }

    public String getName() {
        return name;
    }

    @Override
    public int hashCode() {
        return hashCode == null ? hashCode = (int) StringUtil.getHash(fullname) : hashCode;
    }

    private double[] getMasses(String s) {
        ArrayList<Double> list = new ArrayList<Double>();
        list.add(0.);
        list.add(parentMass);
        double sum = 0;
        for (int i = 0; i < s.length(); ++i) {
            char c = s.charAt(i);
            if (Character.isLetter(c)) {
                sum += AA_MASS_ARRAY[c];
            } else if (c == '[') {
                StringBuilder builder = new StringBuilder();
                while (s.charAt(++i) != ']') {
                    builder.append(s.charAt(i));
                }
                sum += Double.parseDouble(builder.toString());
            }
            list.add(sum);
            list.add(parentMass - sum);
        }
        Collections.sort(list);
        double[] ans = new double[list.size()];
        for (int i = 0; i < ans.length; ++i) {
            ans[i] = list.get(i);
        }
        return ans;
    }

    public static double getStringMass(String string) {
        double sum = 0;
        char[] s = string.toCharArray();
        for (int i = 0; i < s.length; ++i) {
            char c = s[i];
            if (Character.isLetter(c)) {
                sum += AA_MASS_ARRAY[c];
            } else if (c == '[') {
                StringBuilder builder = new StringBuilder();
                while (s[++i] != ']') {
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
        bestShift = 1e200;
        return Math.max(getMaxMatch(path.getEdges(), path.beginMass, path.length()), getMaxMatch(path.getReversedEdges(), MassUtil.convertIonsType(path.beginMass + path.getMass(), parentMass), path.length()));
    }

    private int getMaxMatch(Edge[] edges, double beginMass, int length) {
        int maxScore = 0;
        double mass = 0;
        for (int i = 0; i < protein.length(); ++i) {
            if (protein.charAt(i) == '[') {
                StringBuilder builder = new StringBuilder();
                while (protein.charAt(++i) != ']') {
                    builder.append(protein.charAt(i));
                }
                mass += Double.parseDouble(builder.toString());
                continue;
            }
            int matchedEdges = 0;
            j:
            for (int j = 0, pos = i; j < edges.length && pos < protein.length(); ++j) {
                int len = 0;
                decoding:
                for (char[] decoding : edges[j].getDecodings()) {
                    len = decoding.length;
                    for (int k = 0; k < decoding.length && pos + decoding.length <= protein.length(); ++k) {
                        if (protein.charAt(pos + k) != decoding[k]) {
                            continue decoding;
                        }
                    }
                    pos += len;
                    ++matchedEdges;
                    continue j;
                }
                pos += len;
            }
            maxScore = Math.max(matchedEdges, maxScore);
            double diff = beginMass - mass;
            if (matchedEdges == length && Math.abs(bestShift) > Math.abs(diff)) {
                bestShift = diff;
            }
            if (Character.isLetter(protein.charAt(i))) {
                mass += AA_MASS_ARRAY[protein.charAt(i)];
            }
        }
        return maxScore;
    }

    public double getLastBestShift() {
        return bestShift;
    }

    public int getId() {
        return id;
    }

    public static String getShortName(String name) {
        String subName = name.substring(name.indexOf("NP_"));
        return subName.substring(0, subName.indexOf('|'));
    }
}
