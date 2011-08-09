package ru.spbau.ablab.tagfinder;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import ru.spbau.ablab.tagfinder.util.ComparablePair;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.FastScanner;

public class Multifiler {
	public static void main(String[] args) throws FileNotFoundException {
		FastScanner scanner = new FastScanner(new File(ConfigReader.getProperty("ALIGN_SPECTRA_FILE")));
		while (scanner.skipLine("BEGIN IONS")) {
			int id = scanner.getNextIntProperty("SCANS");
			scanner.getNextDoubleProperty("PRECURSOR_MASS");
			ArrayList<ComparablePair<Double, Double>> list = new ArrayList<ComparablePair<Double,Double>>();
			for (String s; !(s = scanner.nextLine()).equals("END IONS");) {
				String[] split = s.split("	+");
				list.add(new ComparablePair<Double, Double>(Double.parseDouble(split[0]), Double.parseDouble(split[1])));
			}
			Collections.sort(list);
			PrintWriter writer = new PrintWriter(new File(StatisticsGenerator.ENVELOPES_DIR +"/scan_" + id + ".vir"));
			writer.println(list.size());
			for (ComparablePair<?, ?> pair : list) {
				writer.println(pair.a);
			}
			writer.println(list.size());
			for (ComparablePair<?, ?> pair : list) {
				writer.println(pair.b);
			}
			writer.close();
		}
	}
}
