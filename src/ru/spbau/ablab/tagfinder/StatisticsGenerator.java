package ru.spbau.ablab.tagfinder;

import java.io.FileNotFoundException;
import java.util.*;

import ru.spbau.ablab.tagfinder.path.edges.AAEdge;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.GapEdge;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.*;
import ru.spbau.ablab.tagfinder.util.io.FastScanner;
import ru.spbau.ablab.tagfinder.util.io.HtmlWriter;

public class StatisticsGenerator implements Runnable {
    public static final boolean USE_DEFAULT_FILENAME = ConfigReader.getBooleanProperty("USE_DEFAULT_FILENAME");
    public static final boolean SCORE_BY_LENGTH = ConfigReader.getBooleanProperty("SCORE_BY_LENGTH");
    public static final boolean EDGE_OF_TWO_AA = ConfigReader.getBooleanProperty("EDGE_OF_TWO_AA");
    public static final int MAX_TAG_LENGTH = ConfigReader.getIntProperty("MAX_TAG_LENGTH");
    public static final double MASS_EPS = ConfigReader.getDoubleProperty("MASS_EPS");
    public static final int MAX_PATHS = ConfigReader.getIntProperty("MAX_PATHS");
    public static final boolean DOUBLE_MASSES = ConfigReader.getBooleanProperty("DOUBLE_MASSES");
    public static final String OUTPUT_FILE = ConfigReader.getProperty("OUTPUT_FILE");


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

    private Database database;

    public static void main(String[] args) {
        new StatisticsGenerator().run();
    }

    @Override
    public void run() {
        runFor(null);
    }

    public void runFor(ArrayList<Integer> scanIds) {
        try {
            database = Database.getDatabase();
            scanIds = database.filter(scanIds);
            HtmlWriter writer = new HtmlWriter(USE_DEFAULT_FILENAME ? OUTPUT_FILE : getOutputFilename());
            printStatistics(writer, scanIds);
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(13);
        }
        System.out.println("done");
    }

    private void printStatistics(HtmlWriter writer, ArrayList<Integer> scanIds) {
        writer.printOpenTag("table", "cellpadding=0 cellspacing=20");
        writer.printHeader();
        int[] found = new int[MAX_PATHS];
        int[][] count = new int[MAX_PATHS][MAX_TAG_LENGTH + 1];
        for (int id : scanIds) {
            processScan(id, writer, count, found);
        }
        writer.printFrequences(count);
        writer.printRatio(found, scanIds.size());
        writer.printCloseTag("table");
    }


    private void processScan(int id, HtmlWriter writer, int[][] count, int[] found) {
        writer.printOpenTag("tr");
        writer.printThTaggedValue(id);
        writer.printThTaggedValue(database.getSpectrum(id).envelopes.length);
        Protein protein = database.getProtein(id);
        TreeSet<Path> paths = getAllPaths(id);
        writer.printThTaggedValue(paths.size());
        printProteinDBMatches(writer, paths);
        writer.printOpenTh();
        writer.printf("%.2E", database.getEValue(id));
        writer.printCloseTh();
        printTags(writer, count, found, protein, paths);
        writer.printCloseTag("tr");
        writer.flush();
        System.out.println(id + " ok");
    }

    private void printTags(HtmlWriter writer, int[][] count, int[] found, Protein protein, TreeSet<Path> paths) {
        int pathN = 0;
        boolean foundFirst = false;
        for (Path path : paths) {
            if (pathN == MAX_PATHS) {
                break;
            }
            writer.printOpenTag("td");
            writer.printOpenTag("div", "align=center");
            writer.printTagPrefix("span");
            boolean notPrint = false;
            if (protein.contains(path)) {
                if (!foundFirst) {
                    writer.printTagSuffix("style=\"color:red;font-weight:bold\"");
                    ++found[pathN];
                    foundFirst = true;
                } else {
                    writer.printTagSuffix("style=\"color:magenta\"");
                }
            } else {
                for (int mask = 1; mask < 4; ++mask) {
                    Path subPath = path.subPath(mask & 1, path.edges - mask / 2);
                    Edge[] edges = path.getEdges();
                    if (subPath.length() > 0 && protein.contains(subPath)) {
                        writer.printTagSuffix("style=\"color:blue\"");
                        if (mask % 2 == 1) {
                            writer.println(edges[0]);
                        }
                        writer.printCloseTag("span");
                        writer.printTaggedValue("span", subPath, "style=\"color:red\"");
                        if (mask / 2 == 1) {
                            writer.printTaggedValue("span", edges[path.edges - 1], "style=\"color:blue\"");
                        }
                        writer.printTaggedValue("span", protein.getMaxMatch(subPath), "style=\"color:blue\"");
                        notPrint = true;
                        break;
                    }
                }
                if (!notPrint) {
                    writer.printTagSuffix("style=\"color:blue\"");
                }
            }
            ++count[pathN][path.length()];
            if (!notPrint) {
                writer.println(path + " " + protein.getMaxMatch(path));
                writer.printCloseTag("span");
            }
            writer.printCloseTag("div");
            writer.printCloseTag("td");
            ++pathN;
        }
    }

    private void printProteinDBMatches(HtmlWriter writer, TreeSet<Path> paths) {
        writer.printThTaggedValue(database.getMatchedProteinsNumber(paths));
    }

    private TreeSet<Path> getAllPaths(int id) {
        TreeMap<Path, Double> bestScore = new TreeMap<Path, Double>(Path.LENGTH_FIRST_COMPARATOR);
        Spectrum spectrum = database.getSpectrum(id);
        double parentMass = spectrum.parentMass;
        Envelope[] envelopes = spectrum.envelopes;
        Envelope[] reversedEnvelopes = new Envelope[envelopes.length];
        for (int i = 0; i < envelopes.length; ++i) {
            reversedEnvelopes[i] = envelopes[envelopes.length - 1 - i].getReversed(parentMass);
        }
        Spectrum reversedSpectrum = new Spectrum(id, reversedEnvelopes, parentMass);
        boolean[] usedEnvelopes = new boolean[envelopes.length];
        for (int i = 0; i < envelopes.length; ++i) {
            ArrayList<Double> list = new ArrayList<Double>();
            addTags(spectrum, reversedSpectrum, i, new Path(new Edge[0], envelopes[i].score), bestScore, list, null, usedEnvelopes);
        }

        return new TreeSet<Path>(bestScore.keySet());
    }

    private void addTags(Spectrum spectrum, Spectrum reversedSpectrum, int envelopeId, Path path, TreeMap<Path, Double> bestScore, ArrayList<Double> peaks, Double parentMassCorrection, boolean[] usedEnvelopes) {
        if (path.length() > 2) {
            Double d = (d = bestScore.get(path)) == null ? Double.NEGATIVE_INFINITY : d;
            bestScore.put(path, Math.max(d, path.score));
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
            findEdges(spectrum, reversedSpectrum, path, bestScore, peaks, parentMassCorrection, currentMass, edge, usedEnvelopes);
        }
        if (EDGE_OF_TWO_AA) {
            for (GapEdge[] edges : GAP_EDGES) {
                for (Edge edge : edges) {
                    findEdges(spectrum, reversedSpectrum, path, bestScore, peaks, parentMassCorrection, currentMass, edge, usedEnvelopes);
                }
            }
        }
        usedEnvelopes[arrayIndex] = false;
        peaks.remove(peaks.size() - 1);
    }

    private void findEdges(Spectrum spectrum, Spectrum reversedSpectrum, Path path, TreeMap<Path, Double> bestScore, ArrayList<Double> peaks, Double parentMassCorrection, double currentMass, Edge edge, boolean[] usedEnvelopes) {
        double needMass = currentMass + edge.getMass();
        for (int next = spectrum.getFirstMatchingEnvelopeIndex(currentMass, needMass, edge.getMass()); next < spectrum.envelopes.length && MassComparator.edgeMatches(currentMass, spectrum.envelopes[next].getMass(parentMassCorrection), edge.getMass()); ++next) {
            Envelope nextEnvelope = spectrum.envelopes[next];
            if (usedEnvelopes[next]) {
                continue;
            }
            Path newPath = path.append(edge, nextEnvelope.score);
            addTags(spectrum, reversedSpectrum, next, newPath, bestScore, peaks, parentMassCorrection, usedEnvelopes);
        }
        if (DOUBLE_MASSES) {
            for (int next = reversedSpectrum.getFirstMatchingEnvelopeIndex(currentMass, needMass, edge.getMass(), parentMassCorrection); next < reversedSpectrum.envelopes.length && MassComparator.edgeMatches(currentMass, reversedSpectrum.envelopes[next].getMass(parentMassCorrection), edge.getMass(), spectrum.parentMass, parentMassCorrection); ++next) {
                Envelope nextEnvelope = reversedSpectrum.envelopes[next];
                if (usedEnvelopes[spectrum.envelopes.length - 1 - next]) {
                    continue;
                }
                Path newPath = path.append(edge, nextEnvelope.score);
                Double newMassCorrection = parentMassCorrection;
                if (newMassCorrection == null) {
                    newMassCorrection = currentMass + edge.getMass() - spectrum.parentMass + nextEnvelope.getMass(parentMassCorrection);
                }
                addTags(spectrum, reversedSpectrum, -next - 1, newPath, bestScore, peaks, newMassCorrection, usedEnvelopes);
            }
        }
    }

    public String getOutputFilename() {
        return (EDGE_OF_TWO_AA ? "2" : "1") + "_" + ((int) (MassComparator.ERROR_THRESHOLD * 1e6)) + "ppm_" + (SCORE_BY_LENGTH ? "len" : "score") + "_" + (ConfigReader.getBooleanProperty("ALIGN") ? "virt" : "exp") + ".html";
    }
}