package ru.spbau.ablab.tagfinder.database.proteindb;

import edu.ucsd.msalign.spec.id.IdCompEValue;
import ru.spbau.ablab.tagfinder.Protein;
import ru.spbau.ablab.tagfinder.StatisticsGenerator;
import ru.spbau.ablab.tagfinder.database.Database;
import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.util.AKAutomaton;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.pairs.Pair;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

import static ru.spbau.ablab.tagfinder.database.Database.*;

public abstract class ProteinDb {
    protected static final String FASTA_PATH = ConfigReader.getProperty("FASTA_PATH");

    protected Map<Integer, Protein> scanIdToProtein;
    protected Map<String, Protein> fullnameToProtein;
    protected Map<String, Protein> nameToProtein;
    protected Map<Integer, Protein> idToProtein;
    protected HashSet<Protein> proteins;

    private HashSet<String> lastMatchedProteins;

    public HashSet<String> getLastMatchedProteins() {
        return lastMatchedProteins;
    }

    public HashSet<Protein> getProteins() {
        return proteins;
    }

    public Protein getProteinById(int id) {
        Protein protein = idToProtein.get(id);
        if (protein == null) {
            throw new RuntimeException("No protein id=" + id);
        }
        return protein;
    }

    public Protein getProteinByScanId(int id) {
        Protein protein = scanIdToProtein.get(id);
        if (protein == null) {
            throw new RuntimeException("No protein for scan id=" + id);
        }
        return protein;
    }

    public Protein getProteinByName(String name) {
        Protein protein = nameToProtein.get(name);
        if (protein == null) {
            throw new RuntimeException("No protein name=" + name);
        }
        return protein;
    }

    public Protein getProteinByFullName(String fullName) {
        Protein protein = nameToProtein.get(fullName);
        if (protein == null) {
            throw new RuntimeException("No protein name=" + fullName);
        }
        return protein;
    }

    public int getMatchedProteinsNumber(Collection<Path> paths) {
        assert paths.size() <= StatisticsGenerator.MAX_PATHS;
        lastMatchedProteins = new HashSet<String>();
        AKAutomaton automaton = new AKAutomaton(paths);
        for (Protein protein : proteins) {
            if (automaton.acceptsString(protein.getString())) {
                lastMatchedProteins.add(protein.getFullname());
            }
        }
        return lastMatchedProteins.size();
    }

    private HashMap<Integer, Protein> bestMatches = new HashMap<Integer, Protein>();

    public Protein getBestMatchedProtein(int id) {
        Protein protein;
        if (!bestMatches.containsKey(id)) {
            bestMatches.put(id, protein = Database.getInstance().getBestMatch(id));
        } else {
            protein = bestMatches.get(id);
        }
        return protein;
    }

    private IdCompEValue comp;
    private HashMap<Integer, Pair<Double, Protein>> alignResults;
    private PrintWriter alignWriter;
    {
        try {
            alignWriter = new PrintWriter(new File("pair_list"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public Pair<Double, Protein> getBestFromAlign(int id, Collection<Path> paths) {
//        if (alignResults == null) {
//            try {
//                FastScanner scanner = new FastScanner(new File("result_list"));
//                alignResults = new HashMap<Integer, Pair<Double, Protein>>();
//                for (String s; (s = scanner.nextLine()) != null && s.trim().length() > 0; ) {
//                    String[] ss = s.split("\t+");
//                    int scanId = Integer.parseInt(ss[5]);
//                    Protein protein = fullnameToProtein.get(ss[9]);
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
                comp = new IdCompEValue(FASTA_PATH, ALIGN_SPECTRA_FILE, 15);
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
        double eValue = Double.POSITIVE_INFINITY;
        Protein best = null;
        ArrayList<Protein> proteins = new ArrayList<Protein>();
        AKAutomaton automaton = new AKAutomaton(paths);
        for (Map.Entry<String, Protein> entry : nameToProtein.entrySet()) {
            if (automaton.acceptsString(entry.getValue().getString())) {
                proteins.add(entry.getValue());
            }
        }
        for (Protein protein : proteins) {
            try {
                System.out.println("starting for " + id + " " + protein.getName());
                comp.compEValue(id, protein.getName());
                alignWriter.println(id + " " + protein.getName());
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
//        System.err.println("id = " + id + " evalue = " + eValue);
        return new Pair<Double, Protein>(eValue, best);
    }
}
