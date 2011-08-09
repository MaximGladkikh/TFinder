package ru.spbau.ablab.tagfinder.spectrum;

import ru.spbau.ablab.tagfinder.util.MassComparator;

public class Envelope implements Comparable<Envelope> {
	public final double mass;
	public final double score;
	public final double intensity;

	public Envelope(double mass, double score, double intensity) {
		super();
		this.mass = mass;
		this.score = score;
		this.intensity = intensity;
	}

	@Override
	public int compareTo(Envelope o) {
		return MassComparator.compare(mass, o.mass);
	}

	public String toString() {
		return "" + mass;
	}
}
