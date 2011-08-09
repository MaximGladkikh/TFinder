package ru.spbau.ablab.tagfinder.spectrum;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;

import ru.spbau.ablab.tagfinder.StatisticsGenerator;
import ru.spbau.ablab.tagfinder.util.FastScanner;

public class Spectrum {
	public final Envelope[] envelopes;
	public final int id;
	public final double parentMass;
    public final double eValue;

	public Spectrum(int id, File file, double eValue) throws FileNotFoundException {
		this.id = id;
		FastScanner scanner = new FastScanner(file);
		Envelope[] tenvelopes = new Envelope[scanner.getNextIntValue("ENVELOPE_NUMBER")];
		this.parentMass = scanner.getNextDoubleValue("MONOISOTOPIC_MASS");
		for (int i = 0; i < tenvelopes.length; ++i) {
			tenvelopes[i] = readEnvelope(scanner);
		}
		if (StatisticsGenerator.doubleMasses) {
			tenvelopes = doubleEnvelopes(tenvelopes);
		}
		envelopes = tenvelopes;
		Arrays.sort(envelopes);
        this.eValue = eValue;
		scanner.close();
	}

	private Envelope[] doubleEnvelopes(Envelope[] envelopes2) {
		Envelope[] envelopes = new Envelope[envelopes2.length * 2];
		for (int i = 0; i < envelopes2.length; ++i) {
			Envelope envelope = envelopes2[i];
			envelopes[i * 2] = envelope;
			envelopes[i * 2 + 1] = new Envelope(parentMass - envelope.mass, envelope.score, envelope.intensity);
		}
		return envelopes;
	}

	public Spectrum(int id, Envelope[] envelopes, double parentMass, double eValue) {
		this.id = id;
		this.envelopes = envelopes;
        this.eValue = eValue;
		this.parentMass = parentMass;
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

	public Envelope getClosest(double mass) {
		int l = 0;
		int r = envelopes.length - 1;
		int closest = r;
		while (l <= r) {
			int med = (l + r) / 2;
			if (envelopes[med].mass >= mass) {
				closest = med;
				r = med - 1;
			} else {
				l = med + 1;
			}
		}
		if (closest == 0 && envelopes.length > 1 && Math.abs(envelopes[0].mass - mass) > Math.abs(envelopes[1].mass - mass)) {
			closest = 1;
		} else if (closest > 0 && Math.abs(envelopes[closest - 1].mass - mass) < Math.abs(envelopes[closest].mass - mass)) {
			--closest;
		}
		return envelopes[closest];
	}
}