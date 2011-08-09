package ru.spbau.ablab.tagfinder;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.FastScanner;

public class StatisticsGeneratorVirtual extends StatisticsGeneratorExperimental {
	private static final String ALIGN_RESULT_FILE = ConfigReader.getProperty("ALIGN_RESULT_FILE");
	private static final String ALIGN_SPECTRA_FILE = ConfigReader.getProperty("ALIGN_SPECTRA_FILE");

	protected TreeMap<Integer, Spectrum> getSpectra() throws FileNotFoundException {
		StatisticsGeneratorExperimental.doubleMasses = false;
		scanToSpectrum = super.getSpectra();
		HashMap<Integer, Integer> idToScan = getIdToScanMap();
		FastScanner scanner = new FastScanner(new File(ALIGN_RESULT_FILE));
		TreeMap<Integer, Spectrum> ans = new TreeMap<Integer, Spectrum>();
		while (scanner.hasNextLine()) {
			scanner.skipLine("BEGIN PRSM");
			int id = idToScan.get(scanner.getNextIntProperty("SPECTRUM_ID"));
			double parentMass = 0;
			if (scanToSpectrum.containsKey(id)) {
				parentMass = scanToSpectrum.get(id).parentMass;
			}
			scanner.skipLine("BEGIN MATCH_PAIR");
			ArrayList<Envelope> envelopes = new ArrayList<Envelope>();
			for (String s; !(s = scanner.nextLine().trim()).equals("END MATCH_PAIR");) {
				String[] ss = s.split("	+");
				double mass = Double.parseDouble(ss[3]);
				if (!scanToSpectrum.containsKey(id)) {
					continue;
				}
				Envelope closestExp = scanToSpectrum.get(id).getClosest(mass);
				mass = ss[4].equals("B") ? mass : (parentMass - mass);
				envelopes.add(new Envelope(mass, closestExp.score, closestExp.intensity));
			}
			ans.put(id, new Spectrum(id, envelopes.toArray(new Envelope[envelopes.size()]), parentMass));
			scanner.skipLine("END PRSM");
		}
		scanner.close();
		return ans;
	}

	private HashMap<Integer, Integer> getIdToScanMap() throws FileNotFoundException {
		FastScanner scanner = new FastScanner(new File(ALIGN_SPECTRA_FILE));
		HashMap<Integer, Integer> ans = new HashMap<Integer, Integer>();
		while (scanner.skipLine("BEGIN IONS")) {
			ans.put(scanner.getNextIntProperty("ID"), scanner.getNextIntProperty("SCANS"));
		}
		scanner.close();
		return ans;
	}
}