package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.Database;
import ru.spbau.ablab.tagfinder.util.MassUtil;
import ru.spbau.ablab.tagfinder.util.io.HtmlWriter;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Set;
import java.util.TreeSet;

public class RemoveMatchStatistics extends StatisticsGenerator {
    static boolean FILTER_DB = true;
    private static final int MAX_PATHS_FOR_SHIFTS = ConfigReader.getIntProperty("MAX_PATHS_FOR_SHIFTS");

    public static void main(String[] args) {
        new RemoveMatchStatistics().run();
    }

    public void run() {
        try {
            database = Database.getInstance();
            ArrayList<Integer> ids = database.filter(true);
            HtmlWriter writer = new HtmlWriter("removedpeaks.html");
            for (int id : ids) {
                removeSharedPeaks(id);
            }
            RemovedScanProcessor processor = new RemovedScanProcessor(false);
            printStatistics(writer, ids, processor);
            System.out.println("done " + processor.getProcessedScansNumber());
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(13);
        }
    }

    public class RemovedScanProcessor extends MatchedScanProcessor {
        public RemovedScanProcessor(boolean skip) {
            super(skip, !FILTER_DB);
            printShifts = true;
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
            success = false;
            Spectrum experimental = database.getSpectrum(id);
            boolean[] delete = new boolean[experimental.envelopes.length];
//        Protein protein = database.getProteinFromTable(id);
            Protein protein = database.getAlignedProtein(id);
//            Protein protein = database.getProtein(id);
//            removeTheoreticalSpectrum(experimental, delete, protein);
//            success = removeTags(experimental, delete, protein);
            Set<Path> paths = TagGenerator.getAllPaths(database, id);
            TreeSet<Double> shifts = new TreeSet<Double>(MassUtil.MASS_COMPARATOR);
            int count = 0;
            for (Path path : paths) {
                if (protein.contains(path)) {
                    shifts.add(protein.getBestLastAlignment());
                }
                if (++count == MAX_PATHS_FOR_SHIFTS) {
                    break;
                }
            }
            shifts.add(0.);
            System.err.println(shifts.size());
            for (double d : shifts) {
                success |= removeTheoreticalSpectrum(experimental, delete, protein, -d);
                success |= removeTheoreticalSpectrum(experimental, delete, protein, d);
            }
            ArrayList<Envelope> envelopes = new ArrayList<Envelope>();
            for (int i = 0; i < delete.length; ++i) {
                if (!delete[i]) {
                    envelopes.add(experimental.envelopes[i]);
                }
            }
            int removed = 0;
            for (boolean b : delete) {
                if (b) {
                    ++removed;
                }
            }
            System.err.println("removed " + removed + " peaks");
            Spectrum newSpectrum = new Spectrum(id, envelopes.toArray(new Envelope[envelopes.size()]), experimental.parentMass);
            database.setSpectrum(id, newSpectrum);
        } while (success);
        System.err.println("cleaned " + id);
    }


    private boolean removeTheoreticalSpectrum(Spectrum experimental, boolean[] delete, Protein protein, double d) {
        if (protein == null) {
            return false;
        }
        boolean deleted = false;
        for (int i = 0; i < protein.masses.length; ++i) {
            double mass = protein.masses[i] + d;
            if (mass == 0 || MassUtil.compare(mass, experimental.parentMass) == 0) {
                continue;
            }
            deleted |= removePeak(experimental, delete, mass);
            deleted |= removePeak(experimental, delete, experimental.parentMass - mass);
        }
        return deleted;
    }

    private boolean removePeak(Spectrum experimental, boolean[] delete, double mass) {
        int index = experimental.getClosestIndex(mass);
        if (experimental.envelopes.length == 0) {
            return  false;
        }
        Envelope closest = experimental.envelopes[index];
        return MassUtil.sameForDeletion(mass, closest.getMass()) && !delete[index] && (delete[index] = true);
    }
}
