package ru.spbau.ablab.tagfinder.spectrum;

import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.io.FastScanner;
import ru.spbau.ablab.tagfinder.util.MassComparator;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;

import static ru.spbau.ablab.tagfinder.util.MassComparator.FIRST_BY_EDGE_ERROR;

public class Spectrum {
    public static final boolean REMOVE_DUPLICATES = ConfigReader.getBooleanProperty("REMOVE_DUPLICATES");
    public final Envelope[] envelopes;
    public final int id;
    public final double parentMass;
    public final double massEps;

    public Spectrum(int id, File file) throws FileNotFoundException {
        this.id = id;
        FastScanner scanner = new FastScanner(file);
        Envelope[] envelopes = new Envelope[scanner.getNextIntValue("ENVELOPE_NUMBER")];
        this.parentMass = scanner.getNextDoubleValue("MONOISOTOPIC_MASS");
        for (int i = 0; i < envelopes.length; ++i) {
            envelopes[i] = readEnvelope(scanner);
        }
        Arrays.sort(envelopes);
        this.envelopes = unique(envelopes);
        massEps = FIRST_BY_EDGE_ERROR * parentMass;
        scanner.close();
    }

    private Envelope[] unique(Envelope[] envelopes) {
        if (!REMOVE_DUPLICATES) {
            return envelopes;
        }
        ArrayList<Envelope> ans = new ArrayList<Envelope>();
        for (int i = 0, j; i < envelopes.length; i = j) {
            int best = i;
            for (j = i; j < envelopes.length && envelopes[i].compareTo(envelopes[j]) == 0; ++j) {
                if (envelopes[best].score < envelopes[j].score) {
                    best = j;
                }
            }
            ans.add(envelopes[best]);
        }
        return ans.toArray(new Envelope[ans.size()]);
    }

    public Spectrum(int id, Envelope[] envelopes, double parentMass) {
        this.id = id;
        this.envelopes = envelopes;
        this.parentMass = parentMass;
        massEps = FIRST_BY_EDGE_ERROR * parentMass;
        Arrays.sort(envelopes);
    }

    private Envelope readEnvelope(FastScanner scanner) {
        scanner.skipLine("BEGIN ENVELOPE");
        double score = scanner.getNextDoubleValue("SCORE");
        double mass = scanner.getNextDoubleValue("REAL_MONO_MASS");
        double intensity = scanner.getNextDoubleValue("REAL_INTE_SUM");
        Envelope envelope = new Envelope(mass, score, intensity);
        scanner.skipLine("END ENVELOPE");
        return envelope;
    }

    public int getFirstMatchingEnvelopeIndex(double firstMass, double needMass, double edgeMass, Double massCorrection) {
        if (massCorrection == null) {
            int index = getClosestIndex(needMass);
            while (index >= 0 && Math.abs(Math.abs(envelopes[index].getMass(massCorrection) - firstMass) - edgeMass) < massEps) {
                --index;
            }
            ++index;
            return index;
        }
        return getFirstMatchingCorrectedEnvelopeIndex(firstMass, needMass, edgeMass, massCorrection);
    }

    public int getFirstMatchingEnvelopeIndex(double firstMass, double needMass, double edgeMass) {
        return getFirstMatchingCorrectedEnvelopeIndex(firstMass, needMass, edgeMass, 0);
    }

    private int getFirstMatchingCorrectedEnvelopeIndex(double firstMass, double needMass, double edgeMass, double delta) {
        int index = getClosestIndex(needMass);
        while (index >= 0 && MassComparator.edgeMatches(firstMass, envelopes[index].getMass() + delta, edgeMass)) {
            --index;
        }
        ++index;
        if (index >= envelopes.length || !MassComparator.edgeMatches(firstMass, envelopes[index].getMass() + delta, edgeMass)) {
            return envelopes.length;
        }
        return index;
    }

    public Envelope getClosest(double mass) {
        return envelopes[getClosestIndex(mass)];
    }

    public int getClosestIndex(double mass) {
        int l = 0;
        int r = envelopes.length - 1;
        int closest = r;
        while (l <= r) {
            int med = (l + r) / 2;
            if (envelopes[med].getMass(null) >= mass) {
                closest = med;
                r = med - 1;
            } else {
                l = med + 1;
            }
        }
        if (closest == 0 && envelopes.length > 1 && Math.abs(envelopes[0].getMass() - mass) > Math.abs(envelopes[1].getMass() - mass)) {
            closest = 1;
        } else if (closest > 0 && Math.abs(envelopes[closest - 1].getMass() - mass) < Math.abs(envelopes[closest].getMass() - mass)) {
            --closest;
        }
        return closest;
    }
}