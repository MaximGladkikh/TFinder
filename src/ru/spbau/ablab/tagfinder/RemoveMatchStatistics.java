package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.Database;
import ru.spbau.ablab.tagfinder.util.MassComparator;
import ru.spbau.ablab.tagfinder.util.io.HtmlWriter;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Set;

public class RemoveMatchStatistics extends StatisticsGenerator {
    static boolean FILTER_DB = true;

    public static void main(String[] args) {
        new RemoveMatchStatistics().run();
    }

    public void run() {
        try {
            database = Database.getInstance();
            ArrayList<Integer> ids = database.filter(false);
            HtmlWriter writer = new HtmlWriter("removedpeaks.html");
            int count = 0;
            int all = 0;
            for (int id : ids) {
                removeSharedPeaks(id);
                ++all;
                if (database.getProteinFromTable(id).getName().equals(database.getProteinName(id))) {
                    ++count;
                }
            }
            RemovedScanProcessor processor = new RemovedScanProcessor(false);
            printStatistics(writer, ids, processor);
            System.out.println("done " + processor.getProcessedScansNumber());
            System.out.println(count + " of " + all + " prsms matched");
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(13);
        }
    }

    public class RemovedScanProcessor extends MatchedScanProcessor {
        public RemovedScanProcessor(boolean skip) {
            super(skip, !FILTER_DB);
        }

        @Override
        protected void printMatchedProteinsStats(HtmlWriter writer, int id) {
            Set<String> lastMatchedProteins = database.getLastMatchedProteins();
            if (lastMatchedProteins == null) {
                return;
            }
            writer.printTaggedValue("td", lastMatchedProteins.contains(database.getProteinFromTable(id).getName()));
            super.printMatchedProteinsStats(writer, id);
        }
    }

    private void removeSharedPeaks(int id) {
        boolean success;
        do {
            Spectrum experimental = database.getSpectrum(id);
            boolean[] delete = new boolean[experimental.envelopes.length];
        //            Protein protein = database.getProtein(id);
//            System.err.println(Arrays.toString(experimental.envelopes));
//            System.err.println(Arrays.toString(protein.masses));
//            System.err.println(protein.getName());
            Protein protein = database.getProtein(id);
//        System.err.println(protein != null && protein.getName().equals(database.getProteinName(id)));
//        if (protein == null || !protein.getName().equals(database.getProteinFromTable(id).getName())) {
//            System.err.println(((protein == null) ? null: protein.getName()) + " " + database.getProteinFromTable(id).getName());
//            System.err.println(protein.getName() + " " + protein.getMatchScore(experimental));
//            System.err.println(database.getProteinName(id) + " " + database.getProtein(id).getMatchScore(experimental));
//            throw new AssertionError("" + id);
//        }
//        assert protein == null || protein.getName().equals(database.getProteinName(id));
//            removeTheoreticalSpectrum(experimental, delete, protein);
//            success = removeTags(experimental, delete, protein);
            ArrayList<Envelope> envelopes = new ArrayList<Envelope>();
            for (int i = 0; i < delete.length; ++i) {
                if (!delete[i]) {
                    envelopes.add(experimental.envelopes[i]);
                }
            }
            Spectrum newSpectrum = new Spectrum(id, envelopes.toArray(new Envelope[envelopes.size()]), experimental.parentMass);
            database.setSpectrum(id, newSpectrum);
        } while (false && success);
        System.err.println("cleaned " + id);
    }

    private boolean removeTags(Spectrum experimental, boolean[] delete, Protein protein) {
        if (protein == null) {
            return false;
        }
        boolean found = false;
        Set<Path> paths = TagGenerator.getAllPaths(database, experimental.id);
        String proteinString = protein.getString();
        for (Path path : paths) {
            found |= removeTag(experimental, delete, protein, path, proteinString);
//            found |= removeTag(experimental, delete, protein, path.getReversed(), proteinString);
        }
        return found;
    }

    private boolean removeTag(Spectrum experimental, boolean[] delete, Protein protein, Path path, String proteinString) {
        Edge[] edges = path.getEdges();
//        StringBuilder builder = new StringBuilder();
//        for (Edge e : edges) {
//            builder.append(e);
//        }
//        String s = builder.toString();
        if (!protein.contains(path)) {
            return false;
        }
        double offset = path.beginMass;
        boolean found = false;
        for (int i = 0; i <= edges.length; ++i) {
            int index = experimental.getClosestIndex(offset);
            if (MassComparator.sameForDeletion(experimental.envelopes[index].getMass(), offset)) {
                delete[index] = true;
                found = true;
            }
            if (i < edges.length) {
                offset += edges[i].getMass();
            }
        }
        return found;
    }

    private boolean removeTag(Spectrum experimental, boolean[] delete, Protein protein, Edge[] edges, String proteinString) {
        StringBuilder builder = new StringBuilder();
        for (Edge e : edges) {
            builder.append(e);
        }
        String s = builder.toString();
        boolean found = false;
        for (int pos = proteinString.indexOf(s); pos < delete.length && pos >= 0; pos = proteinString.indexOf(s, pos + 1)) {
            double offset = protein.masses[pos];
            for (int i = 0; i <= edges.length; ++i) {
                int index = experimental.getClosestIndex(offset);
                if (MassComparator.sameForDeletion(experimental.envelopes[index].getMass(), offset)) {
                    delete[index] = true;
                    found = true;
                }
                if (i < edges.length) {
                    offset += edges[i].getMass();
                }
            }
        }
        return found;
    }

    private void removeTheoreticalSpectrum(Spectrum experimental, boolean[] delete, Protein protein) {
        if (protein == null) {
            return;
        }
        for (int i = 1; i + 1 < protein.masses.length; ++i) {
            double mass = protein.masses[i];
            if (mass == 0 || MassComparator.compare(mass, experimental.parentMass) == 0) {
                continue;
            }
            removePeak(experimental, delete, mass);
            removePeak(experimental, delete, experimental.parentMass - mass);
        }
    }

    private void removePeak(Spectrum experimental, boolean[] delete, double mass) {
        int index = experimental.getClosestIndex(mass);
        Envelope closest = experimental.envelopes[index];
        if (MassComparator.compare(mass, closest.getMass()) == 0) {
            delete[index] = true;
        }
    }
}
