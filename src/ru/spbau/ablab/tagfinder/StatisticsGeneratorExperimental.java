package ru.spbau.ablab.tagfinder;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import ru.spbau.ablab.tagfinder.path.AAEdge;
import ru.spbau.ablab.tagfinder.path.Edge;
import ru.spbau.ablab.tagfinder.path.GapEdge;
import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.FastScanner;
import ru.spbau.ablab.tagfinder.util.HtmlWriter;
import ru.spbau.ablab.tagfinder.util.MassComparator;

public class StatisticsGeneratorExperimental implements Runnable {
	public static final boolean GET_ENTIRE_PROTEIN = ConfigReader.getBoolProperty("GET_ENTIRE_PROTEIN");
	public static final boolean SCORE_BY_LENGH = ConfigReader.getBoolProperty("SCORE_BY_LENGTH");
	public static final boolean EDGE_OF_TWO_AA = ConfigReader.getBoolProperty("EDGE_OF_TWO_AA");
	public static final int MAX_TAG_LENGTH = ConfigReader.getIntProperty("MAX_TAG_LENGTH");
	public static final double MASS_EPS = ConfigReader.getDoubleProperty("MASS_EPS");
	public static final int MAX_PATHS = ConfigReader.getIntProperty("MAX_PATHS");

	public static final String ENVELOPES_DIR = ConfigReader.getProperty("ENVELOPES_DIR");
	public static final String SPECTRUM_FILE_SUFFIX = ConfigReader.getProperty("SPECTRUM_FILE_SUFFIX");
	public static final String OUTPUT_FILE = ConfigReader.getProperty("OUTPUT_FILE");

	public static final char[] AA_LET;
	public static final double[] AA_MONO_MASS;
	public static final double[] AA_AVG_MASS;

	public static boolean doubleMasses = ConfigReader.getBoolProperty("DOUBLE_MASSES");

	private static final String MATCHES_PATH = ConfigReader.getProperty("MATCHES_PATH");
	private static final String PROTEIN_DB_PATH = ConfigReader.getProperty("PROTEIN_DB_PATH");
	private static final String SPECTRUM_FILE_PREFIX = ConfigReader.getProperty("SPECTRUM_FILE_PREFIX");

	private static final String MASS_LIST = ConfigReader.getProperty("MASS_LIST");
	static {
		AA_LET = new char[19];
		AA_MONO_MASS = new double[AA_LET.length];
		AA_AVG_MASS = new double[AA_LET.length];
		FastScanner scanner = new FastScanner(MASS_LIST);
		for (int i = 0; i < AA_LET.length; ++i) {
			AA_LET[i] = scanner.nextToken().charAt(0);
			AA_MONO_MASS[i] = scanner.nextDouble();
			AA_AVG_MASS[i] = scanner.nextDouble();
		}
	}

	private String[][] matchesData;
	private Map<Integer, Integer> scanToRow;
	protected Map<Integer, Spectrum> scanToSpectrum;
	private Map<String, Protein> proteinDB;

	@Override
	public void run() {
		runFor(null);
	}

	public void runFor(ArrayList<Integer> scanIds) {
		try {
			initDB();
			scanToSpectrum = getSpectra();
			scanIds = intersect(scanIds == null ? scanToSpectrum.keySet() : scanIds);
			HtmlWriter writer = new HtmlWriter(OUTPUT_FILE);
			printStatistics(writer, scanToSpectrum, scanIds);
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(13);
		}
		System.out.println("done");
	}

	public ArrayList<Integer> intersect(Collection<Integer> scanToSpectrum) {
		ArrayList<Integer> ans = new ArrayList<Integer>();
		HashSet<String> differentNames = new HashSet<String>();
		for (int id : scanToSpectrum) {
			if (scanToRow.containsKey(id)) {
				if (differentNames.add(getProteinName(id))) {
					ans.add(id);
				}
			}
		}
		ArrayList<String> toRemove = new ArrayList<String>();
		for (String s : proteinDB.keySet()) {
			if (!differentNames.contains(s)) {
				toRemove.add(s);
			}
		}
		for (String s : toRemove) {
			proteinDB.remove(s);
		}
		Collections.sort(ans);
		return ans;
	}

	public void printStatistics(HtmlWriter writer, Map<Integer, Spectrum> scanToSpectrum, ArrayList<Integer> scanIds) {
		writer.printOpenTag("table", "cellpadding=0 cellspacing=20");
		writer.printHeader();
		int[] found = new int[MAX_PATHS];
		int[][] count = new int[MAX_PATHS][MAX_TAG_LENGTH + 1];
		for (int id : scanIds) {
			processScan(id, writer, count, found, scanToSpectrum);
		}
		writer.printFrequences(count);
		writer.printRatio(found, scanIds.size());
		writer.printCloseTag("table");
	}

	public void initDB() throws FileNotFoundException {
		proteinDB = getProteinDB();
		matchesData = getMatches();
	}

	private HashMap<String, Protein> getProteinDB() throws FileNotFoundException {
		FastScanner scanner = new FastScanner(new File(PROTEIN_DB_PATH));
		HashMap<String, Protein> ans = new HashMap<String, Protein>();
		String peptideName = null;
		StringBuilder protein = null;
		for (String line; (line = scanner.nextLine()) != null;) {
			if (line.charAt(0) == '>') {
				if (peptideName != null) {
					ans.put(peptideName, new Protein(protein.toString().trim()));
				}
				peptideName = line.substring(1);
				protein = new StringBuilder();
			} else {
				protein.append(line);
			}
		}
		return ans;
	}

	private void processScan(int id, HtmlWriter writer, int[][] count, int[] found, Map<Integer, Spectrum> scanToSpectrum) {
		writer.printOpenTag("tr");
		writer.printThTaggedValue(id);
		writer.printThTaggedValue(scanToSpectrum.get(id).envelopes.length);
		Protein protein = getProtein(id);
		TreeSet<Path> paths = getAllPaths(id, protein, scanToSpectrum);
		writer.printThTaggedValue(paths.size());
		printMatches(writer, paths);
		writer.printOpenTh();
		writer.printf("%.2E", getEValue(id));
		writer.printCloseTh();
		printTags(writer, count, found, protein, paths);
		writer.printCloseTag("tr");
		writer.flush();
		System.out.println(id + " ok");
	}

	private void printTags(HtmlWriter writer, int[][] count, int[] found, Protein protein, TreeSet<Path> paths) {
		int pathN = 0;
		boolean foundFirst = false;
		for (Path path : paths) {
			if (pathN == MAX_PATHS) {
				break;
			}
			writer.printOpenTag("td");
			writer.printOpenTag("div", "align=center");
			writer.printTagPrefix("span");
			boolean notPrint = false;
			if (protein.contains(path)) {
				if (!foundFirst) {
					writer.printTagSuffix("style=\"color:red;font-weight:bold\"");
					++found[pathN];
					foundFirst = true;
				} else {
					writer.printTagSuffix("style=\"color:magenta\"");
				}
			} else {
				Path subPath = path.subPath(1, path.edges.length - 1);
				if (subPath.length() > 2 && protein.contains(subPath)) {
					writer.printTagSuffix("style=\"color:blue\"");
					writer.println(path.edges[0]);
					writer.printCloseTag("span");
					writer.printTaggedValue("span", subPath, "style=\"color:red\"");
					writer.printTaggedValue("span", path.edges[path.edges.length - 1] + " " + protein.getMaxMatch(subPath), "style=\"color:blue\"");
					notPrint = true;
				} else {
					writer.printTagSuffix("style=\"color:blue\"");
				}
			}
			++count[pathN][path.length()];
			if (!notPrint) {
				writer.println(path + " " + protein.getMaxMatch(path));
				writer.printCloseTag("span");
			}
			writer.printCloseTag("div");
			writer.printCloseTag("td");
			++pathN;
		}
	}

	private void printMatches(HtmlWriter writer, TreeSet<Path> paths) {
		int pathN = 0;
		HashSet<String> matchedProteins = new HashSet<String>();
		for (Path path : paths) {
			for (Map.Entry<String, Protein> entry : proteinDB.entrySet()) {
				if (entry.getValue().contains(path)) {
					matchedProteins.add(entry.getKey());
				}
			}
			++pathN;
			if (pathN == MAX_PATHS) {
				break;
			}
		}
		writer.printThTaggedValue(matchedProteins.size());
	}

	private Protein getProtein(int id) {
		if (GET_ENTIRE_PROTEIN) {
			return proteinDB.get(getProteinName(id));
		}
		return new Protein(matchesData[scanToRow.get(id)][12]);
	}

	private double getEValue(int id) {
		return Double.parseDouble(matchesData[scanToRow.get(id)][13]);
	}

	private String getProteinName(int id) {
		return matchesData[scanToRow.get(id)][6];
	}

	private TreeSet<Path> getAllPaths(int id, Protein protein, Map<Integer, Spectrum> scanToSpectrum) {
		TreeMap<Path, Double> bestScore = new TreeMap<Path, Double>(Path.LENGTH_FIRST_COMPARATOR);

		Envelope[] envelopes = scanToSpectrum.get(id).envelopes;
		for (int i = 0; i < envelopes.length; ++i) {
			ArrayList<Double> list = new ArrayList<Double>();
			list.add(envelopes[i].mass);
			addTags(envelopes, i, new Path(new Edge[0], envelopes[i].score), bestScore, protein, list);
		}

		return new TreeSet<Path>(bestScore.keySet());
	}

	private void addTags(Envelope[] envelopes, int v, Path path, TreeMap<Path, Double> bestScore, Protein protein, ArrayList<Double> peaks) {
		if (path.length() > 2) {
			Double d = (d = bestScore.get(path)) == null ? Double.NEGATIVE_INFINITY : d;
			bestScore.put(path, Math.max(d, path.score));
		}
		if (path.length() >= MAX_TAG_LENGTH) {
			return;
		}
		for (int i = 0; i < AA_LET.length; ++i) {
			double needMass = envelopes[v].mass + AA_MONO_MASS[i];
			for (int j = v + 1; j < envelopes.length && MassComparator.compare(needMass, envelopes[j].mass) >= 0; ++j) {
                if (MassComparator.edgeMatches(envelopes[v].mass, envelopes[j].mass, AA_MONO_MASS[i])) {
//				if (MassComparator.same(needMass, envelopes[j].mass)) {
					Path newPath = path.append(new AAEdge(AA_LET[i]), envelopes[j].intensity);
					peaks.add(envelopes[j].mass);
					addTags(envelopes, j, newPath, bestScore, protein, peaks);
					peaks.remove(peaks.size() - 1);
				}
			}
		}
		if (EDGE_OF_TWO_AA) {
			for (int i = 0; i < AA_LET.length; ++i) {
				for (int j = 0; j < AA_LET.length; ++j) {
					double needMass = envelopes[v].mass + AA_MONO_MASS[i] + AA_MONO_MASS[j];
					for (int k = v + 1; k < envelopes.length && MassComparator.compare(needMass, envelopes[k].mass) >= 0; ++k) {
                        if (MassComparator.edgeMatches(envelopes[v].mass, envelopes[k].mass, AA_MONO_MASS[i] + AA_MONO_MASS[j])) {
//						if (MassComparator.same(needMass, envelopes[k].mass)) {
							Path newPath = path.append(new GapEdge(needMass - envelopes[v].mass), envelopes[k].intensity);
							peaks.add(envelopes[k].mass);
							addTags(envelopes, k, newPath, bestScore, protein, peaks);
							peaks.remove(peaks.size() - 1);
						}
					}
				}
			}
		}
	}

	protected TreeMap<Integer, Spectrum> getSpectra() throws FileNotFoundException {
		TreeMap<Integer, Spectrum> ans = new TreeMap<Integer, Spectrum>();
		File envDir = new File(ENVELOPES_DIR);
		for (File file : envDir.listFiles()) {
			String name = file.getName();
			if (!name.toLowerCase().startsWith(SPECTRUM_FILE_PREFIX) || !name.toLowerCase().endsWith(SPECTRUM_FILE_SUFFIX)) {
				continue;
			}
			String stringId = name.substring(SPECTRUM_FILE_PREFIX.length(), name.length() - SPECTRUM_FILE_SUFFIX.length());
			int id = Integer.parseInt(stringId);
			Integer rowId = scanToRow.get(id);
			if (rowId == null) {
				continue;
			}
			Spectrum spectrum = new Spectrum(id, file);
			ans.put(id, spectrum);
		}
		return ans;
	}

	private String[][] getMatches() throws FileNotFoundException {
		ArrayList<String[]> ans = new ArrayList<String[]>();
		scanToRow = new TreeMap<Integer, Integer>();
		FastScanner matchesScanner = new FastScanner(new File(MATCHES_PATH));
		while (matchesScanner.hasNextLine()) {
			String[] current = new String[14];
			for (int i = 0; i < current.length; ++i) {
				String s = matchesScanner.nextToken();
				if (s.startsWith("gi")) {
					while (!s.endsWith("]")) {
						s += " " + matchesScanner.nextToken();
					}
				}
				current[i] = s;
			}
			current[12] = current[12].replace("I", "L");
			int scanId = Integer.parseInt(current[1]);
			scanToRow.put(scanId, ans.size());
			ans.add(current);
		}
		matchesScanner.close();
		return ans.toArray(new String[ans.size()][]);
	}
}