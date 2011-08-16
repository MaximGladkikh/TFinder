package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.AAEdge;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.path.edges.GapEdge;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.Database;
import ru.spbau.ablab.tagfinder.util.MassComparator;
import ru.spbau.ablab.tagfinder.util.io.FastScanner;

import java.util.ArrayList;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class TagGenerator {
    public static final boolean DOUBLE_MASSES = ConfigReader.getBooleanProperty("DOUBLE_MASSES");
    public static final int MAX_TAG_LENGTH = ConfigReader.getIntProperty("MAX_TAG_LENGTH");
    public static final int MIN_TAG_LENGTH = ConfigReader.getIntProperty("MIN_TAG_LENGTH");
    public static final boolean EDGE_OF_TWO_AA = ConfigReader.getBooleanProperty("EDGE_OF_TWO_AA");

    public static final char[] AA_LET;
    public static final double[] AA_MONO_MASS;
    public static final double[] AA_AVG_MASS;

    private static final AAEdge[] AA_EDGES;
    private static final GapEdge[][] GAP_EDGES;

    private static final String MASS_LIST = ConfigReader.getProperty("MASS_LIST");

    static {
        AA_LET = new char[19];
        AA_MONO_MASS = new double[AA_LET.length];
        AA_AVG_MASS = new double[AA_LET.length];
        FastScanner scanner = new FastScanner(MASS_LIST);
        for (int i = 0; i < AA_LET.length; ++i) {
            AA_LET[i] = scanner.nextToken().charAt(0);
            AA_MONO_MASS[i] = scanner.nextDouble();
            AA_AVG_MASS[i] = scanner.nextDouble();
        }
        AA_EDGES = new AAEdge[AA_LET.length];
        GAP_EDGES = new GapEdge[AA_LET.length][AA_LET.length];
        for (int i = 0; i < AA_EDGES.length; ++i) {
            AA_EDGES[i] = new AAEdge(AA_LET[i]);
            for (int j = 0; j < AA_EDGES.length; ++j) {
                GAP_EDGES[i][j] = new GapEdge(AA_MONO_MASS[i] + AA_MONO_MASS[j]);
            }
        }
    }

    private static final double[] RATIO_THRESHOLDS;

    static {
        String[] s = ConfigReader.getProperty("RATIO_THRESHOLDS").split(" +");
        RATIO_THRESHOLDS = new double[s.length];
        for (int i = 0; i < RATIO_THRESHOLDS.length; ++i) {
            RATIO_THRESHOLDS[i] = Double.parseDouble(s[i]);
        }
    }

    private static double lastRatio;

    public static double getLastRatio() {
        return lastRatio;
    }

    public static ArrayList<Path> getTopTags(Database database, int id, int n) {
        ArrayList<Path> ans = new ArrayList<Path>();
        for (int i = 0; i < RATIO_THRESHOLDS.length && ans.size() < n; ++i) {
            Set<Path> paths = getAllPaths(database, id, RATIO_THRESHOLDS[i]);
            lastRatio = RATIO_THRESHOLDS[i];
            int pathN = 0;
            paths:
            for (Path path : paths) {
                if (ans.size() >= n || pathN >= StatisticsGenerator.MAX_PATHS) {
                    break;
                }
                ++pathN;
                for (Path stored : ans) {
                    if (stored.toString().equals(path.toString())) {
                        continue paths;
                    }
                }
                ans.add(path);
            }
        }
        return ans;
    }

    public static TreeSet<Path> getAllPaths(Database database, int id) {
        return getAllPaths(database, id, 1);
    }


    public static TreeSet<Path> getAllPaths(Database database, int id, double ratio) {
        TreeMap<Path, Double> bestScore = new TreeMap<Path, Double>(Path.LENGTH_FIRST_COMPARATOR);
        Spectrum spectrum = database.getSpectrum(id);
        double parentMass = spectrum.parentMass;
        Envelope[] envelopes = spectrum.envelopes;
        double minScore = getMinScore(envelopes, ratio);
        Envelope[] reversedEnvelopes = new Envelope[envelopes.length];
        for (int i = 0; i < envelopes.length; ++i) {
            reversedEnvelopes[i] = envelopes[envelopes.length - 1 - i].getReversed(parentMass);
        }
        Spectrum reversedSpectrum = new Spectrum(id, reversedEnvelopes, parentMass);
        boolean[] usedEnvelopes = new boolean[envelopes.length];
        for (int i = 0; i < envelopes.length; ++i) {
            if (envelopes[i].score < minScore) {
                continue;
            }
            ArrayList<Double> list = new ArrayList<Double>();
            addTags(spectrum, reversedSpectrum, i, new Path(new Edge[0], envelopes[i].score), bestScore, list, null, usedEnvelopes, minScore);
        }
        return new TreeSet<Path>(bestScore.keySet());
    }

    private static double getMinScore(Envelope[] envelopes, double ratio) {
        double l = 0;
        double r = 100000;
        for (int it = 0; it < 150; ++it) {
            double med = (l + r) * 0.5;
            int count = 0;
            for (Envelope envelope : envelopes) {
                if (envelope.score >= med) {
                    ++count;
                }
            }
            if (count >= ratio * envelopes.length - 1e-13) {
                l = med;
            } else {
                r = med;
            }
        }
        return r;
    }

    private static void addTags(Spectrum spectrum, Spectrum reversedSpectrum, int envelopeId, Path path, TreeMap<Path, Double> bestScore, ArrayList<Double> peaks, Double parentMassCorrection, boolean[] usedEnvelopes, double minScore) {
        if (path.length() >= MIN_TAG_LENGTH) {
            Double d = (d = bestScore.get(path)) == null ? Double.NEGATIVE_INFINITY : d;
            if (d < path.score) {
                bestScore.put(path, Math.max(d, path.score));
            }
        }
        if (path.length() >= MAX_TAG_LENGTH) {
            return;
        }
        Envelope v = (envelopeId >= 0 ? spectrum.envelopes[envelopeId] : reversedSpectrum.envelopes[-envelopeId - 1]);
        double currentMass = v.getMass(parentMassCorrection);
        peaks.add(currentMass);
        int arrayIndex = envelopeId >= 0 ? envelopeId : (spectrum.envelopes.length + envelopeId);
        usedEnvelopes[arrayIndex] = true;
        for (Edge edge : AA_EDGES) {
            findEdges(spectrum, reversedSpectrum, path, bestScore, peaks, parentMassCorrection, currentMass, edge, usedEnvelopes, minScore);
        }
        if (EDGE_OF_TWO_AA) {
            for (GapEdge[] edges : GAP_EDGES) {
                for (Edge edge : edges) {
                    findEdges(spectrum, reversedSpectrum, path, bestScore, peaks, parentMassCorrection, currentMass, edge, usedEnvelopes, minScore);
                }
            }
        }
        usedEnvelopes[arrayIndex] = false;
        peaks.remove(peaks.size() - 1);
    }

    private static void findEdges(Spectrum spectrum, Spectrum reversedSpectrum, Path path, TreeMap<Path, Double> bestScore, ArrayList<Double> peaks, Double parentMassCorrection, final double currentMass, Edge edge, boolean[] usedEnvelopes, double minScore) {
        double needMass = currentMass + edge.getMass();
        for (int next = spectrum.getFirstMatchingEnvelopeIndex(currentMass, needMass, edge.getMass()); next < spectrum.envelopes.length && MassComparator.edgeMatches(currentMass, spectrum.envelopes[next].getMass(), edge.getMass()); ++next) {
            Envelope nextEnvelope = spectrum.envelopes[next];
            if (usedEnvelopes[next] || nextEnvelope.score < minScore) {
                continue;
            }
            Path newPath = path.append(edge, nextEnvelope.score);
            addTags(spectrum, reversedSpectrum, next, newPath, bestScore, peaks, parentMassCorrection, usedEnvelopes, minScore);
        }
        if (DOUBLE_MASSES) {
            for (int next = reversedSpectrum.getFirstMatchingEnvelopeIndex(currentMass, needMass, edge.getMass(), parentMassCorrection); next < reversedSpectrum.envelopes.length && MassComparator.edgeMatches(currentMass, reversedSpectrum.envelopes[next].getMass(parentMassCorrection), edge.getMass(), spectrum.parentMass, parentMassCorrection); ++next) {
                Envelope nextEnvelope = reversedSpectrum.envelopes[next];
                if (usedEnvelopes[spectrum.envelopes.length - 1 - next] || nextEnvelope.score < minScore) {
                    continue;
                }
                Path newPath = path.append(edge, nextEnvelope.score * (parentMassCorrection == null ? 0.1 : 1));
                Double newMassCorrection = parentMassCorrection;
                if (newMassCorrection == null) {
                    newMassCorrection = -(nextEnvelope.getMass() - needMass);
                }
                addTags(spectrum, reversedSpectrum, -next - 1, newPath, bestScore, peaks, newMassCorrection, usedEnvelopes, minScore);
            }
        }
    }
}
