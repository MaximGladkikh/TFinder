package ru.spbau.ablab.tagfinder;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.StringTokenizer;

import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.FastScanner;

public class StatForScanGenerator {
	public static void main(String[] args) throws FileNotFoundException {
		String sscans = ConfigReader.getProperty("SCAN_INVEST");
		StringTokenizer tokenizer = new StringTokenizer(sscans);
		ArrayList<Integer> scans = new ArrayList<Integer>();
		while (tokenizer.hasMoreTokens()) {
			scans.add(Integer.parseInt(tokenizer.nextToken()));
		}
		if (ConfigReader.getBoolProperty("ALIGN")) {
			new StatisticsGeneratorVirtual().runFor(scans);
		} else {
			new StatisticsGeneratorExperimental().runFor(scans);
		}
//		mergeTables(scans);
	}

	public static void mergeTables(ArrayList<Integer> scans) throws FileNotFoundException {
		new StatisticsGeneratorExperimental().runFor(scans);
		ArrayList<String> strings = new ArrayList<String>();
		FastScanner scanner = new FastScanner(new File(StatisticsGeneratorExperimental.OUTPUT_FILE));
		for (String s; (s = scanner.nextLine()) != null; ) {
			strings.add(s);
		}
		scanner.close();
		new StatisticsGeneratorVirtual().runFor(scans);
		scanner = new FastScanner(new File(StatisticsGeneratorVirtual.OUTPUT_FILE));
		for (String s; (s = scanner.nextLine()) != null; ) {
			strings.add(s);
		}
		scanner.close();
		PrintWriter writer = new PrintWriter(StatisticsGeneratorVirtual.OUTPUT_FILE);
		for (String s : strings) {
			writer.println(s);
		}
		writer.close();
	}
}