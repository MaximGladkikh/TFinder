package ru.spbau.ablab.tagfinder.util;

import ru.spbau.ablab.tagfinder.Protein;
import ru.spbau.ablab.tagfinder.StatisticsGenerator;
import ru.spbau.ablab.tagfinder.TagGenerator;
import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.io.FastScanner;
import ru.spbau.ablab.tagfinder.util.pairs.Pair;
import ru.spbau.ablab.tagfinder.util.trie.HashTrie;
import ru.spbau.ablab.tagfinder.util.trie.Trie;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class Database {
    public static final String ENVELOPES_DIR = ConfigReader.getProperty("ENVELOPES_DIR");
    public static final String SPECTRUM_FILE_SUFFIX = ConfigReader.getProperty("SPECTRUM_FILE_SUFFIX");
    private static final String PROTEIN_DB_PATH = ConfigReader.getProperty("PROTEIN_DB_PATH");
    private static final String SPECTRUM_FILE_PREFIX = ConfigReader.getProperty("SPECTRUM_FILE_PREFIX");
    private static final String ALIGN_RESULT_FILE = ConfigReader.getProperty("ALIGN_RESULT_FILE");
    private static final String ALIGN_SPECTRA_FILE = ConfigReader.getProperty("ALIGN_SPECTRA_FILE");

    private Map<Integer, Spectrum> scanToSpectrum;
    private Map<String, Protein> nameToProtein;
    private Map<Integer, Protein> sequenceIdToProtein;
    private Map<Integer, Double> idToEValue;
    private Map<Integer, Spectrum> scanToVirtualSpectrum;
    private Map<Integer, Protein> idToMatch;
    private static final double MAX_DISTANCE = ConfigReader.getDoubleProperty("MAX_DISTANCE");
    private static final double UNIDENTIFIED_THRESHOLD = ConfigReader.getDoubleProperty("UNINDENTIFIED_THRESHOLD");

    public Database() throws FileNotFoundException {
        nameToProtein = getNameToProtein();
        scanToSpectrum = getExperimentalSpectra();
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
            Double eValue = idToEValue.get(id);
            if (eValue != null && eValue >= UNIDENTIFIED_THRESHOLD) {
                ans.add(id);
            }
        }
        return ans;
    }

    private HashMap<Integer, Protein> bestMatches = new HashMap<Integer, Protein>();

    public Protein getProtein(int id) {
        Protein protein;
        if (!bestMatches.containsKey(id)) {
            bestMatches.put(id, protein = getBestMatch(id));
        } else {
            protein = bestMatches.get(id);
        }
        return protein;
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
//        HashSet<String> differentNames = new HashSet<String>();
        for (int id : scanIds) {
            Double eValue = idToEValue.get(id);
            if (eValue == null) {
                continue;
            }
            if (eValue <= UNIDENTIFIED_THRESHOLD == matched) {
//                String proteinName = getProteinFromTable(id).getName();
//                assert proteinName != null;
//                differentNames.add(proteinName);
//                if (differentNames.add(proteinName)) {
                    ans.add(id);
//                }
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

    public String getLastMatchedProtein() {
        return lastMatchedProteins.isEmpty() ? null : lastMatchedProteins.iterator().next();
    }

    public Set<String> getLastMatchedProteins() {
        return lastMatchedProteins;
    }

    public int getMatchedProteinsNumber(Collection<Path> paths) {
        int pathN = 0;
        lastMatchedProteins = new HashSet<String>();
        for (Path path : paths) {
            if (path == null) {
                continue;
            }
            for (Pair<Protein, Double> pair : getMatchingProteins(path)) {
                Protein protein = pair.a;
                lastMatchedProteins.add(protein.getName());
            }
            ++pathN;
            if (pathN == StatisticsGenerator.MAX_PATHS) {
                break;
            }
        }
        return lastMatchedProteins.size();
    }

    private Trie trie;

    private Collection<Pair<Protein, Double>> getMatchingProteins(Path path) {
        if (TagGenerator.EDGE_OF_TWO_AA) {
            throw new AssertionError("not implemented");
        }
        if (trie == null) {
            fillTrie();
        }
        return trie.getMatching(path.toString());
    }

    private Protein getBestMatch2(int id) {
        Collection<Path> tags = TagGenerator.getTopTags(this, id, StatisticsGenerator.MAX_PATHS);
        if (tags.isEmpty()) {
            return null;
        }
        Spectrum spectrum = getSpectrum(id);
        ArrayList<Pair<Protein, Pair<Double, Double>>> proteins = new ArrayList<Pair<Protein, Pair<Double, Double>>>();
        for (Path tag : tags) {
            for (Pair<Protein, Double> tagMatch : getMatchingProteins(tag)) {
                Protein protein = tagMatch.a;
                double matchBegin = tagMatch.b;
                double position = tagMatch.b;
                if (Math.abs(position - tag.beginMass) > MAX_DISTANCE) {
                    continue;
                }
                double shift = matchBegin - tag.beginMass;
                proteins.add(new Pair<Protein, Pair<Double, Double>>(protein, new Pair<Double, Double>(1., shift)));
            }
            Path rev = tag.getReversed();
            double tagMass = tag.getMass();
            for (Pair<Protein, Double> tagMatch : getMatchingProteins(rev)) {
                Protein protein = tagMatch.a;
                Double matchBegin = tagMatch.b;
                double position = spectrum.parentMass - tagMatch.b - tagMass;
                if (Math.abs(position - tag.beginMass) > MAX_DISTANCE) {
                    continue;
                }
                double shift = matchBegin + tag.beginMass + tagMass;
//                System.err.println(tag + " " + shift + " " + (matchBegin + tagMass) + " " + (-tag.beginMass + shift) + " " + protein.masses[ArrayUtil.getClosestIndex(protein.masses, (-tag.beginMass + shift))]);
                assert Math.abs(tag.beginMass * (-1) + shift - matchBegin - tagMass) < 1e-1;
                proteins.add(new Pair<Protein, Pair<Double, Double>>(protein, new Pair<Double, Double>(-1., shift)));
            }
        }
        int bestScore = Integer.MIN_VALUE;
        Protein bestProtein = null;

        for (Pair<Protein, Pair<Double, Double>> pair : proteins) {
            Protein protein = pair.a;
            double mult = pair.b.a;
            double shift = pair.b.b;
            int score = protein.getMatchScore(spectrum, mult, shift);
            if (score > bestScore) {
                bestScore = score;
                bestProtein = protein;
            }
        }
        return bestProtein;
    }


    private Protein getBestMatch(int id) {
        Collection<Path> tags = TagGenerator.getTopTags(this, id, StatisticsGenerator.MAX_PATHS * 5);
        if (tags.isEmpty()) {
            return null;
        }
        Spectrum spectrum = getSpectrum(id);
//        Collection<Path> allPaths = TagGenerator.getAllPaths(this, id, 0.15);
        HashMap<Protein, Integer> proteins = new HashMap<Protein, Integer>();
        for (Path tag : tags) {
            for (Pair<Protein, Double> tagMatch : getMatchingProteins(tag)) {
                Protein protein = tagMatch.a;
                double position = tagMatch.b;
                if (Math.abs(position - tag.beginMass) > MAX_DISTANCE) {
                    continue;
                }
                Integer before = (before = proteins.get(protein)) == null ? 0 : before;
                proteins.put(protein, before + 1);
            }
            double mass = tag.getMass();
            for (Pair<Protein, Double> tagMatch : getMatchingProteins(tag.getReversed())) {
                Protein protein = tagMatch.a;
                double position = protein.parentMass - tagMatch.b - mass;
                if (Math.abs(position - tag.beginMass) > MAX_DISTANCE) {
                    continue;
                }
                Integer before = (before = proteins.get(protein)) == null ? 0 : before;
                proteins.put(protein, before + 1);
            }
        }
        int bestScore = Integer.MIN_VALUE;
        Protein bestProtein = null;
        double minDiff = Double.POSITIVE_INFINITY;
        for (Map.Entry<Protein, Integer> entry : proteins.entrySet()) {
            double diff = Math.abs(entry.getKey().parentMass - spectrum.parentMass);
            if (entry.getValue() > bestScore || entry.getValue().equals(bestScore) && minDiff > diff) {
                bestScore = entry.getValue();
                bestProtein = entry.getKey();
                minDiff = diff;
            }
        }
//        System.err.println("bestmatch id=" + id + " score = " + bestScore + " " +(bestProtein == null ? null : bestProtein.getName()));
//        System.err.println(bestProtein != null && bestProtein.getName().equals(getProteinName(id)));
        return bestProtein;
    }

    private void fillTrie() {
        trie = new HashTrie();
        for (Map.Entry<String, Protein> entry : nameToProtein.entrySet()) {
            char[] s = entry.getValue().getString().toCharArray();
            for (int length = TagGenerator.MIN_TAG_LENGTH; length <= TagGenerator.MAX_TAG_LENGTH; ++length) {
                for (int i = 0; i + length <= s.length; ++i) {
                    trie.add(s, i, i + length, entry.getValue(), entry.getValue().masses[i]);
                }
            }
        }
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
            double parentMass = 0;
            if (scanToSpectrum.containsKey(id)) {
                parentMass = scanToSpectrum.get(id).parentMass;
            }
            int proteinId = scanner.getNextIntProperty("SEQUENCE_ID");
            idToMatch.put(id, sequenceIdToProtein.get(proteinId));
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
