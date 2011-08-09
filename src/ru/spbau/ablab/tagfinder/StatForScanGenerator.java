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
        new StatisticsGenerator().runFor(scans);
    }

    public static void mergeTables(ArrayList<Integer> scans) throws FileNotFoundException {
        ConfigReader.setProperty("ALIGN", "false");
        new StatisticsGenerator().runFor(scans);
        ArrayList<String> strings = new ArrayList<String>();
        FastScanner scanner = new FastScanner(new File(StatisticsGenerator.OUTPUT_FILE));
        for (String s; (s = scanner.nextLine()) != null; ) {
            strings.add(s);
        }
        scanner.close();
        ConfigReader.setProperty("ALIGN", "true");
        new StatisticsGenerator().runFor(scans);
        scanner = new FastScanner(new File(StatisticsGenerator.OUTPUT_FILE));
        for (String s; (s = scanner.nextLine()) != null; ) {
            strings.add(s);
        }
        scanner.close();
        PrintWriter writer = new PrintWriter(StatisticsGenerator.OUTPUT_FILE);
        for (String s : strings) {
            writer.println(s);
        }
        writer.close();
    }
}