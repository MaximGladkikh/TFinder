package ru.spbau.ablab.tagfinder.util;

import edu.ucsd.msalign.spec.id.IdCompEValue;
import ru.spbau.ablab.tagfinder.Protein;
import ru.spbau.ablab.tagfinder.StatisticsGenerator;
import ru.spbau.ablab.tagfinder.TagGenerator;
import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.io.FastScanner;
import ru.spbau.ablab.tagfinder.util.pairs.Pair;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

public class Database {
    public static final double WATER_MASS = 18.010565;

    public static final String ENVELOPES_DIR = ConfigReader.getProperty("ENVELOPES_DIR");
    public static final String SPECTRUM_FILE_SUFFIX = ConfigReader.getProperty("SPECTRUM_FILE_SUFFIX");
    private static final String PROTEIN_DB_PATH = ConfigReader.getProperty("PROTEIN_DB_PATH");
    private static final String SPECTRUM_FILE_PREFIX = ConfigReader.getProperty("SPECTRUM_FILE_PREFIX");
    private static final String ALIGN_RESULT_FILE = ConfigReader.getProperty("ALIGN_RESULT_FILE");
    private static final String ALIGN_SPECTRA_FILE = ConfigReader.getProperty("ALIGN_SPECTRA_FILE");
    public static final boolean ALIGN = ConfigReader.getBooleanProperty("ALIGN");

    private Map<Integer, Spectrum> scanToSpectrum;
    private Map<String, Protein> nameToProtein;
    private Map<Integer, Protein> sequenceIdToProtein;
    private Map<Integer, Double> idToEValue;
    private Map<Integer, Protein> idToMatch;
    private Map<Integer, Protein> scanToAlignedProtein;
    private static final double UNIDENTIFIED_THRESHOLD = ConfigReader.getDoubleProperty("UNINDENTIFIED_THRESHOLD");
    private static final Database INSTANCE;
    private static final String ALIGN_RESULT_TABLE_FILE = ConfigReader.getProperty("ALIGN_RESULT_TABLE_FILE");

    static {
        Database database = null;
        try {
            database = new Database();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        INSTANCE = database;
    }

    public static Database getInstance() {
        return INSTANCE;
    }

    private Database() throws FileNotFoundException {
        long startTime = System.currentTimeMillis();
        nameToProtein = getNameToProtein();
        scanToSpectrum = getExperimentalSpectra();
        Map<Integer, Spectrum> scanToVirtualSpectrum = getVirtualSpectra();
        if (ALIGN) {
            scanToSpectrum = scanToVirtualSpectrum;
        }
        scanToAlignedProtein = getAlignedProteins();
        System.out.printf("Database loaded in %.3fms\n", 1e-3 * (System.currentTimeMillis() - startTime));
    }

    private Map<Integer, Protein> getAlignedProteins() {
        HashMap<Integer, Protein> ans = new HashMap<Integer, Protein>();
        try {
            FastScanner scanner = new FastScanner(new File(ALIGN_RESULT_TABLE_FILE));
            for(String s; (s = scanner.nextLine()) != null;) {
                String[] ss = s.split("\t+");
                int id = Integer.parseInt(ss[5].split(" +")[2]);
                String proteinString = ss[ss.length - 5].replace('I', 'L');
                String proteinName = ss[9];
                Protein protein = new Protein(proteinString, proteinName);
                ans.put(id, protein);
//                System.err.println((protein.parentMass - scanToSpectrum.get(id).parentMass));
            }
        } catch (FileNotFoundException e) {
            throw  new RuntimeException(e);
        }
        return ans;
    }

    public Protein getAlignedProtein(int id) {
        return scanToAlignedProtein.get(id);
    }

    public void setSpectrum(int id, Spectrum spectrum) {
        assert spectrum.id == id;
        scanToSpectrum.put(id, spectrum);
    }

    public ArrayList<Integer> getUnmatchedIds() {
        ArrayList<Integer> ans = new ArrayList<Integer>();
        for (int id : scanToSpectrum.keySet()) {
            Double eValue = idToEValue.get(id);
            if (eValue != null && eValue >= UNIDENTIFIED_THRESHOLD) {
                ans.add(id);
            }
        }
        return ans;
    }

    private HashMap<Integer, Protein> bestMatches = new HashMap<Integer, Protein>();

    public String getProteinShortName(int id) {
        Protein protein = getProtein(id);
        return protein == null ? null : protein.getShortName();
    }

    public Protein getProtein(int id) {
        Protein protein;
        if (!bestMatches.containsKey(id)) {
            bestMatches.put(id, protein = getBestMatch(id));
        } else {
            protein = bestMatches.get(id);
        }
        return protein;
    }

    private PrintWriter alignWriter;

    {
        try {
            alignWriter = new PrintWriter(new File("pair_list"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }


        private IdCompEValue comp;
    private HashMap<Integer, Pair<Double, Protein>> alignResults;

    public Pair<Double, Protein> getBestFromAlign(int id, Collection<Path> paths) {
//        if (alignResults == null) {
//            try {
//                FastScanner scanner = new FastScanner(new File("result_list"));
//                alignResults = new HashMap<Integer, Pair<Double, Protein>>();
//                for (String s; (s = scanner.nextLine()) != null && s.trim().length() > 0; ) {
//                    String[] ss = s.split("\t+");
//                    int scanId = Integer.parseInt(ss[5]);
//                    Protein protein = nameToProtein.get(ss[9]);
//                    double eValue = Double.parseDouble(ss[ss.length - 3]);
//                    Double prev = alignResults.get(scanId) == null ? Double.POSITIVE_INFINITY : alignResults.get(scanId).a;
//                    if (eValue < prev) {
//                        alignResults.put(scanId, new Pair<Double, Protein>(eValue, protein));
//                    }
//                }
//            } catch (FileNotFoundException e) {
//                throw new RuntimeException(e);
//            }
//        }
//        return alignResults.get(id);
        if (comp == null) {
            try {
                comp = new IdCompEValue(PROTEIN_DB_PATH, ALIGN_SPECTRA_FILE, 15);
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
        double eValue = Double.POSITIVE_INFINITY;
        Protein best = null;
        ArrayList<Protein> proteins = new ArrayList<Protein>();
//        proteins.add(getProteinFromTable(id));
        AKAutomaton automaton = new AKAutomaton(paths);
        for (Map.Entry<String, Protein> entry : nameToProtein.entrySet()) {
            if (automaton.acceptsString(entry.getValue().getString())) {
                proteins.add(entry.getValue());
            }
        }
        for (Protein protein : proteins) {
            try {
                System.err.println("starting for " + id + " " + protein.getShortName());
                comp.compEValue(id, protein.getShortName());
                alignWriter.println(id + " " + protein.getShortName());
                alignWriter.flush();
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        double value = comp.getEValue(i, j);
                        if (value > 1e-50 && value < eValue) {
                            eValue = value;
                            best = protein;
                        }
                    }
                }
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
        System.err.println("id = " + id + " evalue = " + eValue);
//        return null;
        return new Pair<Double, Protein>(eValue, best);
    }

    public Protein getProteinFromTable(int id) {
        return idToMatch.get(id);
    }

    public double getEValue(int id) {
        Double d;
        return (d = idToEValue.get(id)) == null ? 0 : d;
    }

    public String getProteinName(int id) {
        Protein protein = getProtein(id);
        return protein == null ? null : protein.getName();
    }

    public ArrayList<Integer> filter(boolean matched) {
        Collection<Integer> scanIds = ConfigReader.getIntListProperty("SCANS_INVEST");
        if (scanIds.isEmpty()) {
            scanIds = scanToSpectrum.keySet();
        }
        ArrayList<Integer> ans = new ArrayList<Integer>();
        for (int id : scanIds) {
            Double eValue = idToEValue.get(id);
            if (eValue == null) {
                continue;
            }
            if (eValue <= UNIDENTIFIED_THRESHOLD == matched) {
                ans.add(id);
            }
        }
        Collections.sort(ans);
        return ans;
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
            Spectrum spectrum = new Spectrum(id, file);
            ans.put(id, spectrum);
        }
        return ans;
    }

    private HashSet<String> lastMatchedProteins;

    public Set<String> getLastMatchedProteins() {
        return lastMatchedProteins;
    }

    public int getMatchedProteinsNumber(Collection<Path> paths) {
        assert paths.size() <= StatisticsGenerator.MAX_PATHS;
        lastMatchedProteins = new HashSet<String>();
        AKAutomaton automaton = new AKAutomaton(paths);
        for (Map.Entry<String, Protein> entry : nameToProtein.entrySet()) {
            if (automaton.acceptsString(entry.getValue().getString())) {
                lastMatchedProteins.add(entry.getKey());
            }
        }
//        System.err.println(lastMatchedProteins);
        return lastMatchedProteins.size();
    }

    private Protein getBestMatch(int id) {
        Collection<Path> tags = TagGenerator.getTopTags(this, id, StatisticsGenerator.MAX_PATHS * 5);
        if (tags.isEmpty()) {
            return null;
        }
        Spectrum spectrum = getSpectrum(id);
        int bestScore = Integer.MIN_VALUE;
        Protein bestProtein = null;
        double minDiff = Double.POSITIVE_INFINITY;
        AKAutomaton automaton = new AKAutomaton(tags);
        for (Map.Entry<String, Protein> entry : nameToProtein.entrySet()) {
            Protein protein = entry.getValue();
            int score = automaton.getMatchesNumber(protein.getString(), protein.parentMass);
            double diff = Math.abs(protein.parentMass - spectrum.parentMass);
            if (score > bestScore || score == bestScore && minDiff > diff) {
                bestScore = score;
                bestProtein = protein;
                minDiff = diff;
            }
        }
        return bestProtein;
    }

    public Spectrum getSpectrum(int id) {
        return scanToSpectrum.get(id);
    }

    private TreeMap<Integer, Spectrum> getVirtualSpectra() throws FileNotFoundException {
        assert scanToSpectrum != null;
        idToEValue = new HashMap<Integer, Double>();
        idToMatch = new HashMap<Integer, Protein>();
        HashMap<Integer, Integer> idToScan = getIdToScanMap();
        FastScanner scanner = new FastScanner(new File(ALIGN_RESULT_FILE));
        TreeMap<Integer, Spectrum> ans = new TreeMap<Integer, Spectrum>();
        while (scanner.hasNextLine()) {
            scanner.skipLine("BEGIN PRSM");
            int id = idToScan.get(scanner.getNextIntProperty("SPECTRUM_ID"));
            int proteinId = scanner.getNextIntProperty("SEQUENCE_ID");
            idToMatch.put(id, sequenceIdToProtein.get(proteinId));
            double eValue = scanner.getNextDoubleProperty("E_VALUE");
            idToEValue.put(id, eValue);
            double parentMass = sequenceIdToProtein.get(proteinId).parentMass;
            scanner.skipLine("BEGIN MATCH_PAIR");
            ArrayList<Envelope> envelopes = new ArrayList<Envelope>();
            for (String s; !(s = scanner.nextLine().trim()).equals("END MATCH_PAIR"); ) {
                String[] strings = s.split("	+");
                double mass = Double.parseDouble(strings[3]);
                if (!scanToSpectrum.containsKey(id)) {
                    continue;
                }
                Envelope closestExp = scanToSpectrum.get(id).getClosest(mass);
                mass = strings[4].equals("B") ? mass : (MassUtil.convertIonsType(mass, parentMass));
                envelopes.add(new Envelope(mass, closestExp.score, closestExp.intensity));
            }
            ans.put(id, new Spectrum(id, envelopes.toArray(new Envelope[envelopes.size()]), parentMass));
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

    private HashMap<String, Protein> getNameToProtein() throws FileNotFoundException {
        FastScanner scanner = new FastScanner(new File(PROTEIN_DB_PATH));
        HashMap<String, Protein> ans = new HashMap<String, Protein>();
        String proteinName = null;
        StringBuilder proteinString = null;
        sequenceIdToProtein = new HashMap<Integer, Protein>();
        for (String line; (line = scanner.nextLine()) != null; ) {
            if (line.charAt(0) == '>') {
                if (proteinName != null) {
                    Protein protein = new Protein(proteinString.toString().trim().replace('I', 'L'), proteinName);
                    sequenceIdToProtein.put(sequenceIdToProtein.size(), protein);
                    ans.put(proteinName, protein);
                }
                proteinName = line.substring(1);
                proteinString = new StringBuilder();
            } else {
                assert proteinString != null;
                proteinString.append(line);
            }
        }
        if (proteinName != null) {
            ans.put(proteinName, new Protein(proteinString.toString().trim().replace('I', 'L'), proteinName));
        }
        return ans;
    }
}
