package ru.spbau.ablab.tagfinder.util;

import ru.spbau.ablab.tagfinder.Protein;
import ru.spbau.ablab.tagfinder.StatisticsGenerator;
import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.io.FastScanner;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class Database {
    public static final boolean GET_ENTIRE_PROTEIN = ConfigReader.getBooleanProperty("GET_ENTIRE_PROTEIN");
    public static final String ENVELOPES_DIR = ConfigReader.getProperty("ENVELOPES_DIR");
    public static final String SPECTRUM_FILE_SUFFIX = ConfigReader.getProperty("SPECTRUM_FILE_SUFFIX");
    private static final String MATCHES_PATH = ConfigReader.getProperty("MATCHES_PATH");
    private static final String PROTEIN_DB_PATH = ConfigReader.getProperty("PROTEIN_DB_PATH");
    private static final String SPECTRUM_FILE_PREFIX = ConfigReader.getProperty("SPECTRUM_FILE_PREFIX");
    private static final String ALIGN_RESULT_FILE = ConfigReader.getProperty("ALIGN_RESULT_FILE");
    private static final String ALIGN_SPECTRA_FILE = ConfigReader.getProperty("ALIGN_SPECTRA_FILE");

    private String[][] matchesData;
    private Map<Integer, Integer> scanToRow;
    private Map<Integer, Spectrum> scanToSpectrum;
    private Map<String, Protein> proteinDB;
    private Map<Integer, Double> idToEValue;
    private Map<Integer, Spectrum> scanToVirtualSpectrum;

    public Database(boolean skip) throws FileNotFoundException {
        proteinDB = getProteinDB();
        matchesData = getMatches();
        scanToSpectrum = getExperimentalSpectra(skip);
        scanToVirtualSpectrum = getVirtualSpectra();
        if (ConfigReader.getBooleanProperty("ALIGN")) {
            scanToSpectrum = scanToVirtualSpectrum;
        }
    }

    public Spectrum getVirtualSpectrum(int id) {
        return scanToVirtualSpectrum.get(id);
    }

    public void setSpectrum(int id, Spectrum spectrum) {
        scanToSpectrum.put(id, spectrum);
    }

    public ArrayList<Integer> getUnmatchedIds() {
        ArrayList<Integer> ans = new ArrayList<Integer>();
        for (int id : scanToSpectrum.keySet()) {
            if (!scanToRow.containsKey(id)) {
                ans.add(id);
            }
        }
        return ans;
    }

    public Protein getProtein(int id) {
        if (GET_ENTIRE_PROTEIN) {
            return proteinDB.get(getProteinName(id));
        }
        return new Protein(matchesData[scanToRow.get(id)][12]);
    }

    public double getEValue(int id) {
        Double d;
        return (d = idToEValue.get(id)) == null ? 0 : d;
    }

    public String getProteinName(int id) {
        return matchesData[scanToRow.get(id)][6];
    }

    public ArrayList<Integer> filter(Collection<Integer> scanIds) {
        if (scanIds == null) {
            scanIds = scanToSpectrum.keySet();
        }
        ArrayList<Integer> ans = new ArrayList<Integer>();
        HashSet<String> differentNames = new HashSet<String>();
        for (int id : scanIds) {
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

    private TreeMap<Integer, Spectrum> getExperimentalSpectra(boolean skipUnmatched) throws FileNotFoundException {
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
            if (skipUnmatched && rowId == null) {
                continue;
            }
            if (!skipUnmatched && rowId != null) {
                continue;
            }
            Spectrum spectrum = new Spectrum(id, file);
            ans.put(id, spectrum);
        }
        return ans;
    }

    private HashSet<String> lastMatchedProteins;

    public String getLastMatchedProtein() {
        return lastMatchedProteins.isEmpty() ? null : lastMatchedProteins.iterator().next();
    }

    public Set<String> getLastMatchedProteins() {
        return lastMatchedProteins;
    }

    public int getMatchedProteinsNumber(Collection<Path> paths) {
        int pathN = 0;
        lastMatchedProteins = new HashSet<String>();
        HashSet<String> matchedProteins = new HashSet<String>();
        for (Path path : paths) {
            if (path == null) {
                continue;
            }
            for (Map.Entry<String, Protein> entry : proteinDB.entrySet()) {
                if (entry.getValue().contains(path)) {
                    lastMatchedProteins.add(entry.getKey());
                    matchedProteins.add(entry.getKey());
                }
            }
            ++pathN;
            if (pathN == StatisticsGenerator.MAX_PATHS) {
                break;
            }
        }
        return matchedProteins.size();
    }

    public Spectrum getSpectrum(int id) {
        return scanToSpectrum.get(id);
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

    private TreeMap<Integer, Spectrum> getVirtualSpectra() throws FileNotFoundException {
        assert scanToSpectrum != null;
        idToEValue = new HashMap<Integer, Double>();
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
            idToEValue.put(id, eValue);
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

    private HashMap<String, Protein> getProteinDB() throws FileNotFoundException {
        FastScanner scanner = new FastScanner(new File(PROTEIN_DB_PATH));
        HashMap<String, Protein> ans = new HashMap<String, Protein>();
        String proteinName = null;
        StringBuilder protein = null;
        for (String line; (line = scanner.nextLine()) != null; ) {
            if (line.charAt(0) == '>') {
                if (proteinName != null) {
                    ans.put(proteinName, new Protein(protein.toString().trim().replace('I', 'L')));
                }
                proteinName = line.substring(1);
                protein = new StringBuilder();
            } else {
                assert protein != null;
                protein.append(line);
            }
        }
        if (proteinName != null) {
            ans.put(proteinName, new Protein(protein.toString().trim().replace('I', 'L')));
        }
        return ans;
    }
}
