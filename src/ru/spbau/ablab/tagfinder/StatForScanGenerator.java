package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.io.FastScanner;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class StatForScanGenerator {
    public static void main(String[] args) throws FileNotFoundException {
        if (true) {
            String s1 = "1_5ppm_len_exp.html 0.46601\n" +
                    "\t\n" +
                    "0.05825\n" +
                    "\t\n" +
                    "0.01941\n" +
                    "\t\n" +
                    "0.01941\n" +
                    "\t\n" +
                    "0.03883\n" +
                    "\t\n" +
                    "0.0\n" +
                    "\t\n" +
                    "0.00970\n" +
                    "\t\n" +
                    "0.00970\n" +
                    "\t\n" +
                    "0.0\n" +
                    "\t\n" +
                    "0.0 ";
            String s2 = "1_5ppm_len_virt.html 0.50485\n" +
                    "\t\n" +
                    "0.04854\n" +
                    "\t\n" +
                    "0.00970\n" +
                    "\t\n" +
                    "0.0\n" +
                    "\t\n" +
                    "0.02912\n" +
                    "\t\n" +
                    "0.0\n" +
                    "\t\n" +
                    "0.0\n" +
                    "\t\n" +
                    "0.0\n" +
                    "\t\n" +
                    "0.0\n" +
                    "\t\n" +
                    "0.0 ";
            String s3 = "2_5ppm_len_exp.html  0.38834 \n" +
                    "\t\n" +
                    "0.08737\n" +
                    "\t\n" +
                    "0.03883\n" +
                    "\t\n" +
                    "0.00970\n" +
                    "\t\n" +
                    "0.02912\n" +
                    "\t\n" +
                    "0.00970\n" +
                    "\t\n" +
                    "0.01941\n" +
                    "\t\n" +
                    "0.00970\n" +
                    "\t\n" +
                    "0.01941\n" +
                    "\t\n" +
                    "0.00970 ";
            String s4 = "2_5ppm_len_virt.html 0.60194\n" +
                    "\t\n" +
                    "0.06796\n" +
                    "\t\n" +
                    "0.00970\n" +
                    "\t\n" +
                    "0.01941\n" +
                    "\t\n" +
                    "0.00970\n" +
                    "\t\n" +
                    "0.0\n" +
                    "\t\n" +
                    "0.0\n" +
                    "\t\n" +
                    "0.00970\n" +
                    "\t\n" +
                    "0.0\n" +
                    "\t\n" +
                    "0.0 ";
            for (String s : new String[] {s1, s2, s3, s4}) {
                FastScanner scanner = new FastScanner(s);
                String ss = scanner.nextToken();
                double sum = 0;
                for (int i = 0; i < 10 ; ++i) {
                    sum += scanner.nextDouble();
                }
                System.out.println(ss + " " + sum);
            }
            return;
        }
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