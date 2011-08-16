package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.Database;
import ru.spbau.ablab.tagfinder.util.MassComparator;
import ru.spbau.ablab.tagfinder.util.io.HtmlWriter;

import java.io.FileNotFoundException;
import java.util.ArrayList;

public class DeleteMatchStatistics extends StatisticsGenerator{
    public static void main(String[] args) {
        new DeleteMatchStatistics().run();
    }
    public void run() {
        try {
            database = new Database(true);
            ArrayList<Integer> ids = database.filter(null);
            HtmlWriter writer = new HtmlWriter("removedpeaks.html");
            for (int id : ids) {
                removeSharedPeaks(id);
            }
            int processed = printStatistics(writer, ids, new RemovedScanProcessor(true));
            System.out.println("done " + processed);
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(13);
        }
    }

    public class RemovedScanProcessor extends MatchedScanProcessor {
        public RemovedScanProcessor(boolean skip) {
            super(skip);
        }

        @Override
        protected void printMatchedProteinsStats(HtmlWriter writer, int id) {
            writer.printTaggedValue("td", database.getLastMatchedProteins().contains(database.getProteinName(id)));
            super.printMatchedProteinsStats(writer, id);
        }
    }

    private void removeSharedPeaks(int id) {
        Spectrum experimental = database.getSpectrum(id);
        Spectrum virtual = database.getVirtualSpectrum(id);
        boolean[] delete = new boolean[experimental.envelopes.length];
        for (double mass : database.getProtein(id).masses) {
            if (mass == 0 || MassComparator.compare(mass, experimental.parentMass) == 0) {
                continue;
            }
            removePeak(experimental, delete, mass);
            removePeak(experimental, delete, experimental.parentMass - mass);
        }
        for (Envelope envelope : virtual.envelopes) {
            double mass = envelope.getMass();
            if (mass == 0 || MassComparator.compare(mass, experimental.parentMass) == 0) {
                continue;
            }
            removePeak(experimental, delete, mass);
            removePeak(experimental, delete, experimental.parentMass - mass);
        }
        ArrayList<Envelope> envelopes = new ArrayList<Envelope>();
        for (int i = 0 ; i< delete.length; ++i) {
            if (!delete[i]) {
                envelopes.add(experimental.envelopes[i]);
            }
        }
        Spectrum newSpectrum = new Spectrum(id, envelopes.toArray(new Envelope[envelopes.size()]), experimental.parentMass);
        database.setSpectrum(id, newSpectrum);
    }

    private void removePeak(Spectrum experimental, boolean[] delete, double mass) {
        int index = experimental.getClosestIndex(mass);
        Envelope closest = experimental.envelopes[index];
        if (MassComparator.compare(mass, closest.getMass()) == 0) {
            delete[index] = true;
        }
    }
}
