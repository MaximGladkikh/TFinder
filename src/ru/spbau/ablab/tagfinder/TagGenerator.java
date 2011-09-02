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
import ru.spbau.ablab.tagfinder.util.pairs.ComparablePair;
import ru.spbau.ablab.tagfinder.util.pairs.Pair;

import java.util.*;

public class TagGenerator {
    public static boolean DOUBLE_MASSES = ConfigReader.getBooleanProperty("DOUBLE_MASSES");
    private static boolean doubleMasses = DOUBLE_MASSES;
    public static int MAX_TAG_LENGTH = ConfigReader.getIntProperty("MAX_TAG_LENGTH");
    public static int MIN_TAG_LENGTH = ConfigReader.getIntProperty("MIN_TAG_LENGTH");
    public static final boolean EDGE_OF_TWO_AA = ConfigReader.getBooleanProperty("EDGE_OF_TWO_AA");
    public static final boolean EDGE_OF_THREE_AA = ConfigReader.getBooleanProperty("EDGE_OF_THREE_AA");

    public static final char[] AA_LET;
    public static final double[] AA_MONO_MASS;
    public static final int ALPHABET_SIZE = 19;

    private static final AAEdge[] AA_EDGES;
    private static final GapEdge[] GAP2_EDGES;
    private static final GapEdge[] GAP3_EDGES;
    private static final Edge[][] EDGES;

    private static final ArrayList<Pair<Integer, Pair<Integer, Integer>>> replaceRules = new ArrayList<Pair<Integer, Pair<Integer, Integer>>>();

    private static final String MASS_LIST = ConfigReader.getProperty("MASS_LIST");

    static {
        @SuppressWarnings("unchecked")
        ComparablePair<Double, Character>[] acids = new ComparablePair[ALPHABET_SIZE];
        FastScanner scanner = new FastScanner(MASS_LIST);
        for (int i = 0; i < acids.length; ++i) {
            char c = scanner.nextToken().charAt(0);
            double mass = scanner.nextDouble();
            acids[i] = new ComparablePair<Double, Character>(mass, c);
        }
        Arrays.sort(acids);
        AA_LET = new char[ALPHABET_SIZE];
        AA_MONO_MASS = new double[ALPHABET_SIZE];
        for (int i = 0; i < acids.length; ++i) {
            AA_LET[i] = acids[i].b;
            AA_MONO_MASS[i] = acids[i].a;
        }
        AA_EDGES = new AAEdge[ALPHABET_SIZE];
        ArrayList<GapEdge> gap2Edges = new ArrayList<GapEdge>();
        ArrayList<GapEdge> gap3Edges = new ArrayList<GapEdge>();
        for (int i = 0; i < ALPHABET_SIZE; ++i) {
            AA_EDGES[i] = new AAEdge(AA_LET[i]);
            for (int j = 0; j < ALPHABET_SIZE; ++j) {
                double mass = AA_MONO_MASS[i] + AA_MONO_MASS[j];
                String decoding = "" + AA_LET[i] + AA_LET[j];
                if (!findEdge(gap2Edges, mass, decoding)) {
                    gap2Edges.add(new GapEdge(mass, decoding));
                }
                for (int k = 0; k < ALPHABET_SIZE; ++k) {
                    mass = AA_MONO_MASS[i] + AA_MONO_MASS[j] + AA_MONO_MASS[k];
                    decoding = "" + AA_LET[i] + AA_LET[j] + AA_LET[k];
                    if (!findEdge(gap3Edges, mass, decoding)) {
                        gap3Edges.add(new GapEdge(mass, decoding));
                    }
                }
            }
        }
        GAP2_EDGES = gap2Edges.toArray(new GapEdge[gap2Edges.size()]);
        GAP3_EDGES = gap3Edges.toArray(new GapEdge[gap3Edges.size()]);
        EDGES = EDGE_OF_TWO_AA ? (EDGE_OF_THREE_AA ? new Edge[][]{AA_EDGES, GAP2_EDGES, GAP3_EDGES} : new Edge[][]{AA_EDGES, GAP2_EDGES}) : new Edge[][]{AA_EDGES};

        for (int i = 0; i < ALPHABET_SIZE; ++i) {
            for (int j = 0; j < ALPHABET_SIZE; ++j) {
                for (int k = 0; k < ALPHABET_SIZE; ++k) {
                    if (Math.abs(AA_MONO_MASS[i] - AA_MONO_MASS[j] - AA_MONO_MASS[k]) < 1e-1) {
                        replaceRules.add(new Pair<Integer, Pair<Integer, Integer>>(i, new Pair<Integer, Integer>(j, k)));
                    }
                }
            }
        }
    }

    private static boolean findEdge(ArrayList<GapEdge> edges, double mass, String decoding) {
        boolean found = false;
        for (GapEdge edge : edges) {
            if (MassComparator.compare(edge.getMass(), mass) == 0) {
                edge.addDecoding(decoding);
                found = true;
                break;
            }
        }
        return found;
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

    private static void setDoubleMasses(boolean b) {
        doubleMasses = b;
    }

    public static ArrayList<Path> getTopTags(Database database, int id, int n) {
        ArrayList<Path> ans = new ArrayList<Path>();
        for (int i = 0; i < RATIO_THRESHOLDS.length && ans.size() < n; ++i) {
            fillTopTags(database, id, n, ans, i, false);
            if (DOUBLE_MASSES && ans.size() < n) {
                fillTopTags(database, id, n, ans, i, true);
            }
        }
        applyReplaceRules(database, id, ans);
        return ans;
    }

    private static void applyReplaceRules(Database database, int id, ArrayList<Path> ans) {
        for (int i = 0; i < ans.size(); ++i) {
            Edge[] edges = ans.get(i).getEdges();
            ArrayList<Edge> resultEdges = new ArrayList<Edge>();
            double mass = ans.get(i).beginMass;
            j:
            for (Edge edge : edges) {
                for (Pair<Integer, Pair<Integer, Integer>> rule : replaceRules) {
                    double needMass = mass + AA_EDGES[rule.b.a].getMass();
                    Spectrum spectrum = database.getSpectrum(id);
                    Envelope closest = spectrum.getClosest(needMass);
                    Envelope reversedClosest = spectrum.getClosest(spectrum.parentMass - needMass + Database.WATER_MASS);
                    double mass1 = closest.getMass();
                    double mass2 = spectrum.parentMass - reversedClosest.getMass() + Database.WATER_MASS;
                    if (Math.abs(mass1 - needMass) > Math.abs(mass2 - needMass)) {
                        mass1 = mass2;
                    }
                    if (Character.valueOf(AA_LET[rule.a]).equals(edge.getLetter()) && Math.abs(mass1 - needMass) < MassComparator.ERROR_THRESHOLD * 2) {//MassComparator.compare(closest.getMass(), needMass) == 0) {
                        resultEdges.add(AA_EDGES[rule.b.a]);
                        resultEdges.add(AA_EDGES[rule.b.b]);
                        mass += AA_EDGES[rule.a].getMass();
                        continue j;
                    }
                }
                resultEdges.add(edge);
                mass += edge.getMass();
            }
            ans.set(i, new Path(resultEdges.toArray(new Edge[resultEdges.size()]), ans.get(i).beginMass, ans.get(i).spectrum, ans.get(i).isReversed()));
        }
    }

    private static void fillTopTags(Database database, int id, int n, ArrayList<Path> ans, int thresholdIndex, boolean doubleMasses) {
        setDoubleMasses(doubleMasses);
        Set<Path> paths = getAllPaths(database, id, RATIO_THRESHOLDS[thresholdIndex]);
        lastRatio = RATIO_THRESHOLDS[thresholdIndex];
        int pathN = 0;
        paths:
        for (Path path : paths) {
            if (ans.size() >= n || pathN >= StatisticsGenerator.MAX_PATHS) {
                break;
            }
            ++pathN;
            for (Path stored : ans) {
                if (stored.canBeReversedTo(path)) {
                    continue paths;
                }
            }
            ans.add(path);
        }
        setDoubleMasses(DOUBLE_MASSES);
    }

    public static TreeSet<Path> getAllPaths(Database database, int id) {
        return getAllPaths(database, id, 1);
    }

    private static double getMinScore(Envelope[] envelopes, double ratio) {
        double l = 0;
        double r = 1e5;//Max envelope score
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

    public static TreeSet<Path> getAllPaths(Database database, int id, double ratio) {
        TreeMap<Path, Double> bestScore = new TreeMap<Path, Double>(Path.LENGTH_FIRST_COMPARATOR);//best score for each string
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
            addTags(spectrum, reversedSpectrum, i, new Path(new Edge[0], envelopes[i].score, envelopes[i].getMass(), spectrum, false), bestScore, list, null, usedEnvelopes, minScore);
        }
        return new TreeSet<Path>(bestScore.keySet());
    }

    private static void addTags(Spectrum spectrum, Spectrum reversedSpectrum, int envelopeId, Path path, TreeMap<Path, Double> bestScore, ArrayList<Double> peaks, Double parentMassCorrection, boolean[] usedEnvelopes, double minScore) {
        if (path.length() >= MIN_TAG_LENGTH) {
            Double d = (d = bestScore.get(path)) == null ? Double.NEGATIVE_INFINITY : d;
            if (d < path.score) {
                bestScore.put(path, path.score);
            }
        }
        if (path.toString().equals("SSA")) {
            System.err.print("");
        }
        if (path.length() >= MAX_TAG_LENGTH) {
            return;
        }
        Envelope v = (envelopeId >= 0 ? spectrum.envelopes[envelopeId] : reversedSpectrum.envelopes[-envelopeId - 1]);
        double currentMass = v.getMass(parentMassCorrection);
        peaks.add(currentMass);
        int arrayIndex = envelopeId >= 0 ? envelopeId : (spectrum.envelopes.length + envelopeId);
        usedEnvelopes[arrayIndex] = true;
        for (Edge[] edges : EDGES) {
            for (Edge edge : edges) {
                findEdges(spectrum, reversedSpectrum, path, bestScore, peaks, parentMassCorrection, currentMass, edge, usedEnvelopes, minScore);
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
        if (doubleMasses) {
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
