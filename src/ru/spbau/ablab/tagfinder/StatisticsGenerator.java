package ru.spbau.ablab.tagfinder;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import ru.spbau.ablab.tagfinder.path.AAEdge;
import ru.spbau.ablab.tagfinder.path.Edge;
import ru.spbau.ablab.tagfinder.path.GapEdge;
import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.FastScanner;
import ru.spbau.ablab.tagfinder.util.HtmlWriter;
import ru.spbau.ablab.tagfinder.util.MassComparator;

public class StatisticsGenerator implements Runnable {
    public static final boolean USE_DEFAULT_FILENAME = ConfigReader.getBooleanProperty("USE_DEFAULT_FILENAME");
    public static final boolean GET_ENTIRE_PROTEIN = ConfigReader.getBooleanProperty("GET_ENTIRE_PROTEIN");
    public static final boolean SCORE_BY_LENGTH = ConfigReader.getBooleanProperty("SCORE_BY_LENGTH");
    public static final boolean EDGE_OF_TWO_AA = ConfigReader.getBooleanProperty("EDGE_OF_TWO_AA");
    public static final int MAX_TAG_LENGTH = ConfigReader.getIntProperty("MAX_TAG_LENGTH");
    public static final double MASS_EPS = ConfigReader.getDoubleProperty("MASS_EPS");
    public static final int MAX_PATHS = ConfigReader.getIntProperty("MAX_PATHS");
    public static final boolean DOUBLE_MASSES = /*!ConfigReader.getBooleanProperty("ALIGN") && */ConfigReader.getBooleanProperty("DOUBLE_MASSES");

    public static final String ENVELOPES_DIR = ConfigReader.getProperty("ENVELOPES_DIR");
    public static final String SPECTRUM_FILE_SUFFIX = ConfigReader.getProperty("SPECTRUM_FILE_SUFFIX");
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


    private static final String MATCHES_PATH = ConfigReader.getProperty("MATCHES_PATH");
    private static final String PROTEIN_DB_PATH = ConfigReader.getProperty("PROTEIN_DB_PATH");
    private static final String SPECTRUM_FILE_PREFIX = ConfigReader.getProperty("SPECTRUM_FILE_PREFIX");
    private static final String ALIGN_RESULT_FILE = ConfigReader.getProperty("ALIGN_RESULT_FILE");
    private static final String ALIGN_SPECTRA_FILE = ConfigReader.getProperty("ALIGN_SPECTRA_FILE");


    private String[][] matchesData;
    private Map<Integer, Integer> scanToRow;
    private Map<Integer, Spectrum> scanToSpectrum;
    private Map<Integer, Spectrum> scanToVirtualSpectrum;
    private Map<String, Protein> proteinDB;

    public static void main(String[] args) {
        new StatisticsGenerator().run();
    }

    @Override
    public void run() {
        runFor(null);
    }

    public void runFor(ArrayList<Integer> scanIds) {
        try {
            initDB();
            if (ConfigReader.getBooleanProperty("ALIGN")) {
                scanToSpectrum = scanToVirtualSpectrum;
            }
            scanIds = filter(scanIds == null ? scanToSpectrum.keySet() : scanIds);
            HtmlWriter writer = new HtmlWriter(USE_DEFAULT_FILENAME ? OUTPUT_FILE : getOutputFilename());
            printStatistics(writer, scanToSpectrum, scanIds);
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(13);
        }
        System.out.println("done");
    }

    private ArrayList<Integer> filter(Collection<Integer> scanToSpectrum) {
        ArrayList<Integer> ans = new ArrayList<Integer>();
        HashSet<String> differentNames = new HashSet<String>();
        for (int id : scanToSpectrum) {
            if (scanToRow.containsKey(id)) {
                String proteinName = getProteinName(id);
                assert proteinName != null;
                if (differentNames.add(proteinName)) {
                    ans.add(id);
                }
            }
        }
        proteinDB.keySet().retainAll(differentNames);
        Collections.sort(ans);
        return ans;
    }

    private void printStatistics(HtmlWriter writer, Map<Integer, Spectrum> scanToSpectrum, ArrayList<Integer> scanIds) {
        writer.printOpenTag("table", "cellpadding=0 cellspacing=20");
        writer.printHeader();
        int[] found = new int[MAX_PATHS];
        int[][] count = new int[MAX_PATHS][MAX_TAG_LENGTH + 1];
        for (int id : scanIds) {
            processScan(id, writer, count, found, scanToSpectrum);
        }
        writer.printFrequences(count);
        writer.printRatio(found, scanIds.size());
        writer.printCloseTag("table");
    }

    private void initDB() throws FileNotFoundException {
        proteinDB = getProteinDB();
        matchesData = getMatches();
        scanToSpectrum = getExperimentalSpectra();
        scanToVirtualSpectrum = getVirtualSpectra();
    }

    private TreeMap<Integer, Spectrum> getVirtualSpectra() throws FileNotFoundException {
        assert scanToSpectrum != null;
        HashMap<Integer, Integer> idToScan = getIdToScanMap();
        FastScanner scanner = new FastScanner(new File(ALIGN_RESULT_FILE));
        TreeMap<Integer, Spectrum> ans = new TreeMap<Integer, Spectrum>();
        while (scanner.hasNextLine()) {
            scanner.skipLine("BEGIN PRSM");
            int id = idToScan.get(scanner.getNextIntProperty("SPECTRUM_ID"));
            double parentMass = 0;
            if (scanToSpectrum.containsKey(id)) {
                parentMass = scanToSpectrum.get(id).parentMass;
            }
            double eValue = scanner.getNextDoubleProperty("E_VALUE");
            scanner.skipLine("BEGIN MATCH_PAIR");
            ArrayList<Envelope> envelopes = new ArrayList<Envelope>();
            for (String s; !(s = scanner.nextLine().trim()).equals("END MATCH_PAIR"); ) {
                String[] strings = s.split("	+");
                double mass = Double.parseDouble(strings[3]);
                if (!scanToSpectrum.containsKey(id)) {
                    continue;
                }
                Envelope closestExp = scanToSpectrum.get(id).getClosest(mass);
                mass = strings[4].equals("B") ? mass : (parentMass - mass);
                envelopes.add(new Envelope(mass, closestExp.score, closestExp.intensity));
            }
            ans.put(id, new Spectrum(id, envelopes.toArray(new Envelope[envelopes.size()]), parentMass, eValue));
            scanner.skipLine("END PRSM");
        }
        scanner.close();
        return ans;
    }

    private HashMap<Integer, Integer> getIdToScanMap() throws FileNotFoundException {
        FastScanner scanner = new FastScanner(new File(ALIGN_SPECTRA_FILE));
        HashMap<Integer, Integer> ans = new HashMap<Integer, Integer>();
        while (scanner.skipLine("BEGIN IONS")) {
            ans.put(scanner.getNextIntProperty("ID"), scanner.getNextIntProperty("SCANS"));
        }
        scanner.close();
        return ans;
    }

    private HashMap<String, Protein> getProteinDB() throws FileNotFoundException {
        FastScanner scanner = new FastScanner(new File(PROTEIN_DB_PATH));
        HashMap<String, Protein> ans = new HashMap<String, Protein>();
        String peptideName = null;
        StringBuilder protein = null;
        for (String line; (line = scanner.nextLine()) != null; ) {
            if (line.charAt(0) == '>') {
                if (peptideName != null) {
                    ans.put(peptideName, new Protein(protein.toString().trim().replace('I', 'L')));
                }
                peptideName = line.substring(1);
                protein = new StringBuilder();
            } else {
                assert protein != null;
                protein.append(line);
            }
        }
        return ans;
    }

    private void processScan(int id, HtmlWriter writer, int[][] count, int[] found, Map<Integer, Spectrum> scanToSpectrum) {
        writer.printOpenTag("tr");
        writer.printThTaggedValue(id);
        writer.printThTaggedValue(scanToSpectrum.get(id).envelopes.length);
        Protein protein = getProtein(id);
        TreeSet<Path> paths = getAllPaths(id, scanToSpectrum);
        writer.printThTaggedValue(paths.size());
        printProteinDBMatches(writer, paths);
        writer.printOpenTh();
        writer.printf("%.2E", scanToVirtualSpectrum.get(id).eValue);
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
        int pathN = 0;
        HashSet<String> matchedProteins = new HashSet<String>();
        for (Path path : paths) {
            for (Map.Entry<String, Protein> entry : proteinDB.entrySet()) {
                if (entry.getValue().contains(path)) {
                    matchedProteins.add(entry.getKey());
                }
            }
            ++pathN;
            if (pathN == MAX_PATHS) {
                break;
            }
        }
        writer.printThTaggedValue(matchedProteins.size());
    }

    private Protein getProtein(int id) {
        if (GET_ENTIRE_PROTEIN) {
            return proteinDB.get(getProteinName(id));
        }
        return new Protein(matchesData[scanToRow.get(id)][12]);
    }

    private double getEValue(int id) {
        return Double.parseDouble(matchesData[scanToRow.get(id)][13]);
    }

    private String getProteinName(int id) {
        return matchesData[scanToRow.get(id)][6];
    }

    private TreeSet<Path> getAllPaths(int id, Map<Integer, Spectrum> scanToSpectrum) {
        TreeMap<Path, Double> bestScore = new TreeMap<Path, Double>(Path.LENGTH_FIRST_COMPARATOR);
        Spectrum spectrum = scanToSpectrum.get(id);
        double parentMass = spectrum.parentMass;
        Envelope[] envelopes = spectrum.envelopes;
        Envelope[] reversedEnvelopes = new Envelope[envelopes.length];
        for (int i = 0; i < envelopes.length; ++i) {
            reversedEnvelopes[i] = envelopes[envelopes.length - 1 - i].getReversed(parentMass);
        }
        Spectrum reversedSpectrum = new Spectrum(id, reversedEnvelopes, parentMass, spectrum.eValue);
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

    private TreeMap<Integer, Spectrum> getExperimentalSpectra() throws FileNotFoundException {
        TreeMap<Integer, Spectrum> ans = new TreeMap<Integer, Spectrum>();
        File envDir = new File(ENVELOPES_DIR);
        for (File file : envDir.listFiles()) {
            String name = file.getName();
            if (!name.toLowerCase().startsWith(SPECTRUM_FILE_PREFIX) || !name.toLowerCase().endsWith(SPECTRUM_FILE_SUFFIX)) {
                continue;
            }
            String stringId = name.substring(SPECTRUM_FILE_PREFIX.length(), name.length() - SPECTRUM_FILE_SUFFIX.length());
            int id = Integer.parseInt(stringId);
            Integer rowId = scanToRow.get(id);
            if (rowId == null) {
                continue;
            }
            Spectrum spectrum = new Spectrum(id, file, getEValue(id));
            ans.put(id, spectrum);
        }
        return ans;
    }

    private String[][] getMatches() throws FileNotFoundException {
        ArrayList<String[]> ans = new ArrayList<String[]>();
        scanToRow = new TreeMap<Integer, Integer>();
        FastScanner matchesScanner = new FastScanner(new File(MATCHES_PATH));
        while (matchesScanner.hasNextLine()) {
            String[] current = new String[14];
            for (int i = 0; i < current.length; ++i) {
                String s = matchesScanner.nextToken();
                if (s.startsWith("gi")) {
                    while (!s.endsWith("]")) {
                        s += " " + matchesScanner.nextToken();
                    }
                }
                current[i] = s;
            }
            current[12] = current[12].replace("I", "L");
            current[12] = null; //Can't use match string due to unknown format
            int scanId = Integer.parseInt(current[1]);
            scanToRow.put(scanId, ans.size());
            ans.add(current);
        }
        matchesScanner.close();
        return ans.toArray(new String[ans.size()][]);
    }

    public String getOutputFilename() {
        return (EDGE_OF_TWO_AA ? "2" : "1") + "_" + ((int) (MassComparator.ERROR_THRESHOLD * 1e6)) + "ppm_" + (SCORE_BY_LENGTH ? "len" : "score") + "_" + (ConfigReader.getBooleanProperty("ALIGN") ? "virt" : "exp") + ".html";
    }
}