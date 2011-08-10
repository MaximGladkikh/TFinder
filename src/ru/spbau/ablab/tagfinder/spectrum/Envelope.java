package ru.spbau.ablab.tagfinder.spectrum;

import ru.spbau.ablab.tagfinder.util.MassComparator;

public class Envelope implements Comparable<Envelope> {
    public final double score;
    public final double intensity;
    public final boolean reversed;

    private final double mass;

    public Envelope(double mass, double score, double intensity) {
        this(mass, score, intensity, false);
    }

    public Envelope(double mass, double score, double intensity, boolean reversed) {
        this.mass = mass;
        this.score = score;
        this.intensity = intensity;
        this.reversed = reversed;
    }

    public double getMass() {
        return mass;
    }

    public double getMass(Double massCorrection) {
        return mass + (reversed && massCorrection != null ? massCorrection : 0);
    }

    public Envelope getReversed(double parentMass) {
        return new Envelope(parentMass - mass, score, intensity, true);
    }

    @Override
    public int compareTo(Envelope o) {
        return MassComparator.compare(mass, o.mass);
    }

    public String toString() {
        return "" + mass;
    }
}
