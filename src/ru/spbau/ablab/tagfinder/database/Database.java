package ru.spbau.ablab.tagfinder.database;

import ru.spbau.ablab.tagfinder.Protein;
import ru.spbau.ablab.tagfinder.StatisticsGenerator;
import ru.spbau.ablab.tagfinder.TagGenerator;
import ru.spbau.ablab.tagfinder.database.proteindb.ExperimentalProteinDb;
import ru.spbau.ablab.tagfinder.database.proteindb.ProteinDb;
import ru.spbau.ablab.tagfinder.database.proteindb.VirtualProteinDb;
import ru.spbau.ablab.tagfinder.database.spectradb.ExperimentalSpectraDb;
import ru.spbau.ablab.tagfinder.database.spectradb.SpectraDb;
import ru.spbau.ablab.tagfinder.database.spectradb.VirtualSpectraDb;
import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.AKAutomaton;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.MassUtil;
import ru.spbau.ablab.tagfinder.util.io.FastScanner;
import ru.spbau.ablab.tagfinder.util.pairs.ComparablePair;

import javax.xml.stream.events.EntityDeclaration;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class Database {
    public static final double WATER_MASS = 18.010565;

    public static final String ENVELOPES_DIR = ConfigReader.getProperty("ENVELOPES_DIR");
    public static final String SPECTRUM_FILE_SUFFIX = ConfigReader.getProperty("SPECTRUM_FILE_SUFFIX");
    public static final String SPECTRUM_FILE_PREFIX = ConfigReader.getProperty("SPECTRUM_FILE_PREFIX");
    public static final String ALIGN_RESULT_FILE = ConfigReader.getProperty("ALIGN_RESULT_FILE");
    public static final String ALIGN_SPECTRA_FILE = ConfigReader.getProperty("ALIGN_SPECTRA_FILE");
    public static final boolean ALIGN = ConfigReader.getBooleanProperty("ALIGN");
    private static final String MASS_LIST = ConfigReader.getProperty("MASS_LIST");
    public static final char[] AA_LET;
    public static final double[] AA_MONO_MASS;
    public static final int ALPHABET_SIZE = 19;
    public static final double[] AA_MASS_ARRAY = new double[256];
    public static final int[] AA_INDEX = new int[Character.MAX_VALUE];

    private Map<Integer, Double> idToEValue;
    private ProteinDb proteinDb;
    private VirtualProteinDb virtualProteinDb;
    private ExperimentalProteinDb experimentalProteinDb;
    private VirtualSpectraDb virtualSpectraDb;
    private ExperimentalSpectraDb experimentalSpectraDb;
    private SpectraDb spectraDb;

    private static final double UNIDENTIFIED_THRESHOLD = ConfigReader.getDoubleProperty("UNINDENTIFIED_THRESHOLD");
    private static Database INSTANCE;

    static {
        @SuppressWarnings("unchecked")
        ComparablePair<Double, Character>[] acids = new ComparablePair[Database.ALPHABET_SIZE];
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
        for (int i = 0; i < ALPHABET_SIZE; ++i) {
            AA_INDEX[AA_LET[i]] = i;
        }
        Database database = null;
        Arrays.fill(AA_MASS_ARRAY, Double.NaN);
        for (int i = 0; i < AA_LET.length; ++i) {
            AA_MASS_ARRAY[AA_LET[i]] = AA_MONO_MASS[i];
        }
        try {
            database = new Database();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static Database getInstance() {
        return INSTANCE;
    }

    private Database() throws FileNotFoundException {
        INSTANCE = this;
        long startTime = System.currentTimeMillis();
        experimentalSpectraDb = new ExperimentalSpectraDb();
        virtualSpectraDb = new VirtualSpectraDb();
        idToEValue = virtualSpectraDb.getIdToEValue();
        virtualProteinDb = new VirtualProteinDb();
        experimentalProteinDb = new ExperimentalProteinDb();
        if (ALIGN) {
            proteinDb = virtualProteinDb;
            spectraDb = virtualSpectraDb;
        } else {
            proteinDb = experimentalProteinDb;
            spectraDb = experimentalSpectraDb;
        }
        idToEValue = virtualSpectraDb.getIdToEValue();
        System.out.printf("Database loaded in %.3fms\n", 1e-3 * (System.currentTimeMillis() - startTime));
    }

    public ProteinDb getProteinDb() {
        return proteinDb;
    }

    public SpectraDb getSpectraDb() {
        return spectraDb;
    }

    public Envelope getClosestEnvelopeFromExperimental(int id, double mass) {
        Spectrum spectrum = experimentalSpectraDb.getSpectrum(id);
        if (spectrum == null) {
            return null;
        }
        return spectrum.getClosest(mass);
    }

    public Protein getBestMatch(int spectrumId) {
        Collection<Path> tags = TagGenerator.getTopTags(spectrumId, StatisticsGenerator.MAX_PATHS * 5);
        if (tags.isEmpty()) {
            return null;
        }
        Spectrum spectrum = spectraDb.getSpectrum(spectrumId);
        int bestScore = Integer.MIN_VALUE;
        Protein bestProtein = null;
        double minDiff = Double.POSITIVE_INFINITY;
        AKAutomaton automaton = new AKAutomaton(tags);
        for (Protein protein : proteinDb.getProteins()) {
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


    public ArrayList<Integer> getUnmatchedIds() {
        ArrayList<Integer> ans = new ArrayList<Integer>();
        for (int id : spectraDb.getSpectraIds()) {
            Double eValue = idToEValue.get(id);
            if (eValue != null && eValue >= UNIDENTIFIED_THRESHOLD) {
                ans.add(id);
            }
        }
        return ans;
    }

    public double getEValue(int id) {
        Double d;
        return (d = idToEValue.get(id)) == null ? 0 : d;
    }

    public ArrayList<Integer> filter(boolean matched) {
        Collection<Integer> scanIds = ConfigReader.getIntListProperty("SCANS_INVEST");
        if (scanIds.isEmpty()) {
            scanIds = spectraDb.getSpectraIds();
        }
        ArrayList<Integer> ans = new ArrayList<Integer>();
        for (int id : scanIds) {
            Double eValue = idToEValue.get(id);
//            System.err.println(id + " " + eValue);
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

    public Protein getProteinPredictedByAlign(int id) {
        String name = virtualProteinDb.getProteinByScanId(id).getName();
//        return virtualProteinDb.getProteinByScanId(id);
        return experimentalProteinDb.getProteinByName(name);
    }
}
