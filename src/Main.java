import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeSet;

public class Main implements Runnable {
	private static final boolean EDGE_OF_TWO_AA = false;
	private static final int MAX_TAG_LENGTH = 6;
	private static final double MASS_EPS = 1e-2;
	private static final int MAX_PATHS = 10;
	private static final int MAX_LENGTH = 30;

	private static final String SPECTRUM_FILE_PREFIX = "scan_";
	private static final String ENVELOPES_DIR = "env_multiple_mass";
	private static final String OUTPUT_FILE = "output.html";
	private static final String SPECTRUM_FILE_SUFFIX = ".env";

	private static final String MASS_LIST = "A 	 71.03711 	71.0788 R 	156.10111 	156.1875 N 	114.04293 	114.1038 D 	115.02694 	115.0886 C 	103.00919 	103.1388 E 	129.04259 	129.1155 Q 	128.05858 	128.1307 G 	 57.02146 	57.0519 H 	137.05891 	137.1411 L 	113.08406 	113.1594 K 	128.09496 	128.1741 M 	131.04049 	131.1926 F 	 147.06841 	147.1766 P 	 97.05276 	97.1167 S 	87.03203 	87.0782 T 	101.04768 	101.1051 W 	186.07931 	186.2132 Y 	163.06333 	163.1760 V 	99.06841 	99.1326";
	private static final char[] AA_LET;
	private static final double[] AA_MONO_MASS;
	private static final double[] AA_AVG_MASS;

	static {
		AA_LET = new char[19];
		AA_MONO_MASS = new double[AA_LET.length];
		AA_AVG_MASS = new double[AA_LET.length];
		Scanner scanner = new Scanner(MASS_LIST);
		scanner.useLocale(Locale.US);
		for (int i = 0; i < AA_LET.length; ++i) {
			AA_LET[i] = scanner.next().charAt(0);
			AA_MONO_MASS[i] = scanner.nextDouble();
			AA_AVG_MASS[i] = scanner.nextDouble();
		}
	}

	private String[][] matchesData;
	private HashMap<Integer, Integer> scanToRow;
	private ArrayList<Integer> scanIds = new ArrayList<Integer>();
	private HashMap<Integer, Spectrum> scanToSpectrum;

	public static void main(String[] args) {
		new Main().run();
	}

	@Override
	public void run() {
		try {
			matchesData = getMatches();
			scanToSpectrum = getSpectrums();

			PrintWriter writer = new PrintWriter(OUTPUT_FILE);
			writer.println("<table>");

			int[] found = new int[MAX_PATHS];
			int[][] count = new int[MAX_PATHS][MAX_LENGTH];

			int scanCount = 0;
			for (int id : scanIds) {
				processScan(id, writer, count, found);
				++scanCount;
				System.err.println(scanCount);
			}

			printTags(writer, count);
			printLongerThan5(writer, count);
			printRatio(writer, found);

			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(13);
		}
	}

	private void printRatio(PrintWriter writer, int[] found) {
		writer.println("<tr>");
		writer.println("<th> % red </th>");
		for (int i = 0; i < MAX_PATHS; ++i) {
			writer.println("<td>" + 1. * found[i] / scanIds.size() + "</td>");
		}
		writer.println("</tr>");
		writer.println("</table>");
	}

	private void printLongerThan5(PrintWriter writer, int[][] count) {
		writer.println("<tr>");
		writer.println("<th> >= 5</th>");
		for (int i = 0; i < MAX_PATHS; ++i) {
			int count5 = 0;
			for (int j = 5; j < MAX_LENGTH; ++j) {
				count5 += count[i][j];
			}
			writer.println("<td>" + count5 + "</td>");
		}
		writer.println("</tr>");
	}

	private void printTags(PrintWriter writer, int[][] count) {
		for (int len = MAX_LENGTH - 1; len >= 0; --len) {
			writer.println("<tr>");
			writer.println("<th>" + len + "</th>");
			for (int i = 0; i < MAX_PATHS; ++i) {
				writer.println("<td>" + count[i][len] + "</td>");
			}
			writer.println("</tr>");
		}
	}

	private void processScan(int id, PrintWriter writer, int[][] count,
			int[] found) {
		writer.println("<tr>");
		writer.println("<th>" + id + "</th>");
		String protein = getProtein(id);
		TreeSet<Path> paths = getAllPaths(id, protein);

		int pathN = 0;

		HashSet<Integer> matchedProteins = new HashSet<Integer>();
		for (Path path : paths) {
			for (int id2 : scanIds) {
				if (getProtein(id2).contains(path.path)) {
					matchedProteins.add(id2);
				}
			}
			++pathN;
			if (pathN == MAX_PATHS) {
				break;
			}
		}

		writer.println("<th>" + (matchedProteins.size()) + "</th>");

		pathN = 0;
		for (Path path : paths) {
			if (pathN == MAX_PATHS) {
				break;
			}
			writer.println("<td>");
			writer.print("<span");
			boolean notPrint = false;
			if (protein.contains(path.path)) {
				writer.println(" style=\"color:red\">");
				++found[pathN];
			} else {
				String spath = path.path;
				if (spath.length() > 2
						&& protein.contains(spath.substring(1,
								spath.length() - 1))) {
					String subPath = spath.substring(1, spath.length() - 1);
					writer.println(" style=\"color:blue\">");
					writer.println(spath.charAt(0));
					writer.println("</div> <span style=\"color:red\">");
					writer.println(subPath);
					writer.println("</span>");
					writer.println("<span style=\"color:blue\">");
					writer.println(spath.charAt(spath.length() - 1));
					writer.println("</span>");
					notPrint = true;
				} else {
					writer.println(" style=\"color:blue\">");
				}
			}
			++count[pathN][path.path.length()];
			if (!notPrint) {
				writer.println(path.path);
				writer.println(" " + getMaxMatch(protein, path.path));
				writer.println("</span>");
			}
			writer.println("</td>");
			++pathN;
		}
		writer.println("</tr>");

	}

	private int getMaxMatch(String protein, String path) {
		int ans = 0;
		for (int i = 0; i + path.length() <= protein.length(); ++i) {
			int count = 0;
			for (int j = 0; j < path.length(); ++j) {
				if (protein.charAt(i + j) == path.charAt(j)) {
					++count;
				}
			}
			ans = Math.max(ans, count);
		}
		return ans;
	}

	private String getProtein(int id) {
		String protein = matchesData[scanToRow.get(id)][12];
		return protein;
	}

	private TreeSet<Path> getAllPaths(int id, String protein) {
		HashMap<String, Double> bestScore = new HashMap<String, Double>();

		Spectrum.Envelope[] envelopes = scanToSpectrum.get(id).envelopes;
		for (int i = 0; i < envelopes.length; ++i) {
			addTags(envelopes, i, "", bestScore, protein, envelopes[i].score);
		}
		
		TreeSet<Path> ans = new TreeSet<Main.Path>();
		for (Map.Entry<String, Double>  entry : bestScore.entrySet()) {
			ans.add(new Path(entry.getKey(), entry.getValue()));
		}

		return ans;
	}

	private void addTags(Spectrum.Envelope[] envelopes, int v, String string,
			HashMap<String, Double> bestScore, String protein, double score) {
		if (string.length() > 2) {
			Double d = bestScore.get(string);
			if (d == null) {
				d = Double.NEGATIVE_INFINITY;
			}
			bestScore.put(string, Math.max(d, score));
		}
		// if (!protein.contains(string)) {
		// return;
		// }
		if (string.length() >= MAX_TAG_LENGTH) {
			return;
		}
		// System.err.println(v + " " + string);
		for (int i = 0; i < AA_LET.length; ++i) {
			double needMass = envelopes[v].mass + AA_MONO_MASS[i];
			for (int j = v + 1; j < envelopes.length
					&& envelopes[j].mass < needMass + MASS_EPS; ++j) {
				if (Math.abs(needMass - envelopes[j].mass) < MASS_EPS) {
					addTags(envelopes, j, string + AA_LET[i], bestScore, protein,
							score + envelopes[j].score);
				}
			}
		}
		if (EDGE_OF_TWO_AA) {
			for (int i = 0; i < AA_LET.length; ++i) {
				for (int j = 0; j < AA_LET.length; ++j) {
					double needMass = envelopes[v].mass + AA_MONO_MASS[i]
							+ AA_MONO_MASS[j];
					for (int k = v + 1; k < envelopes.length
							&& envelopes[k].mass < needMass + MASS_EPS; ++k) {
						if (Math.abs(needMass - envelopes[k].mass) < MASS_EPS) {
							addTags(envelopes, k, string + AA_LET[i]
									+ AA_LET[j], bestScore, protein, score
									+ envelopes[k].score);
						}
					}
				}
			}
		}
	}

	private HashMap<Integer, Spectrum> getSpectrums()
			throws FileNotFoundException {
		HashMap<Integer, Spectrum> ans = new HashMap<Integer, Spectrum>();
		File envDir = new File(ENVELOPES_DIR);
		for (File file : envDir.listFiles()) {
			String name = file.getName();
			if (!name.toLowerCase().startsWith(SPECTRUM_FILE_PREFIX)) {
				continue;
			}
			String stringId = name.substring(SPECTRUM_FILE_PREFIX.length(),
					name.length() - SPECTRUM_FILE_SUFFIX.length());
			int id = Integer.parseInt(stringId);
			Integer rowId = scanToRow.get(id);
			if (rowId == null) {
				continue;
			}
			scanIds.add(id);
			Spectrum spectrum = new Spectrum(id, file);
			ans.put(id, spectrum);
		}
		Collections.sort(scanIds);
		return ans;
	}

	private String[][] getMatches() throws FileNotFoundException {
		ArrayList<String[]> ans = new ArrayList<String[]>();
		scanToRow = new HashMap<Integer, Integer>();
		HashSet<String> differentNames = new HashSet<String>();
		scanIds = new ArrayList<Integer>();
		Scanner matchesScanner = new Scanner(new File("matches.in"));
		while (matchesScanner.hasNextLine()) {
			String[] current = new String[14];
			for (int i = 0; i < current.length; ++i) {
				String s = matchesScanner.next();
				if (s.startsWith("gi")) {
					while (!s.endsWith("]")) {
						s += " " + matchesScanner.next();
					}
				}
				current[i] = s;
			}
			String name = current[6];
			if (!differentNames.add(name)) {
				continue;
			}
			current[12] = current[12].replace("I", "L");
			int scanId = Integer.parseInt(current[1]);
			scanToRow.put(scanId, ans.size());
			ans.add(current);
		}
		matchesScanner.close();
		return ans.toArray(new String[ans.size()][]);
	}

	static class Path implements Comparable<Path> {
		public String path;
		public double score;
		public static final double SCORE_EPS = 1e-5;

		public Path(String path, double score) {
			this.path = path;
			this.score = score;
		}

		@Override
		public int compareTo(Path o) {
			if (Math.abs(score - o.score) > SCORE_EPS) {
				return score > o.score ? -1 : 1;
			}
			if (path.length() == o.path.length()) {
				return path.compareTo(o.path);
			}
			return -path.length() + o.path.length();
		}
	}
}
