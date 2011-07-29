import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Locale;
import java.util.Scanner;

public class Spectrum {
	public static final double MASS_EPS = 1e-2;

	public Envelope[] envelopes;

	public Spectrum(int id, File file) throws FileNotFoundException {
		this.id = id;
		Scanner scanner = new Scanner(file);
		scanner.useLocale(Locale.US);
		while (!"ENVELOPE_NUMBER".equals(scanner.next()))
			;
		envelopes = new Envelope[scanner.nextInt()];
		
		for (int i = 0; i < envelopes.length; ++i) {
			envelopes[i] = readEnvelope(scanner);
		}
		
		Arrays.sort(envelopes);
		scanner.close();
	}

	private Envelope readEnvelope(Scanner scanner) {
		while (!"BEGIN ENVELOPE".equals(scanner.nextLine())) ;
		
		while (!"SCORE".equals(scanner.next())) ;
		double score = scanner.nextDouble();
		
		while (!"REAL_MONO_MASS".equals(scanner.next())) ;
		Envelope envelope = new Envelope(scanner.nextDouble(), score);
		
		while (!"END ENVELOPE".equals(scanner.nextLine())) ;
		return envelope;
	}

	static class Envelope implements Comparable<Envelope> {
		double mass;
		double score;
		
		public Envelope(double mass, double score) {
			super();
			this.mass = mass;
			this.score = score;
		}

		@Override
		public int compareTo(Envelope o) {
			if (Math.abs(mass - o.mass) < MASS_EPS) {
				return 0;
			}
			return mass < o.mass ? -1 : 1;
		}
	}

	public final int id;
}
