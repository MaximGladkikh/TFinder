package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.AAEdge;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.util.MassComparator;

import static ru.spbau.ablab.tagfinder.TagGenerator.AA_LET;
import static ru.spbau.ablab.tagfinder.TagGenerator.AA_MONO_MASS;

public class Protein {
	public static final double[] AA_MASS_ARRAY = new double[256];
	static {
		for (int i = 0; i < AA_LET.length; ++i) {
			AA_MASS_ARRAY[AA_LET[i]] = AA_MONO_MASS[i];
		}
	}

	private String protein;
	private String revProtein;
	double[] masses;

	public Protein(String s) {
		protein = s;
		masses = getMasses(s);
		revProtein = new StringBuilder(s).reverse().toString();
	}

	private double[] getMasses(String s) {
		double[] ans = new double[s.length() + 1];
		for (int i = 0; i < s.length(); ++i) {
			ans[i + 1] = ans[i] + AA_MASS_ARRAY[s.charAt(i)];
		}
		return ans;
	}

	public boolean contains(Path path) {
		boolean ans = getMaxMatch(path) == path.length();
		assert TagGenerator.EDGE_OF_TWO_AA || ans == (protein.contains(path.toString()) || revProtein.contains(path.toString()));
		return ans;
	}

	public int getMaxMatch(Path path) {
        Edge[] edges = path.getEdges();
        Edge[] rev = new Edge[edges.length];
		for (int i = 0; i < edges.length; ++i) {
			rev[edges.length - 1 - i] = edges[i];
		}
		return Math.max(getMaxMatch(edges, masses), getMaxMatch(rev, masses));
	}

	private int getMaxMatch(Edge[] edges, double[] masses) {
		int maxScore = 0;
		for (int i = 0; i < masses.length - 1; ++i) {
			int pos = i;
			int matchedEdges = 0;
			for (int j = 0; j < edges.length && pos < masses.length - 1; ++j) {
				if (edges[j] instanceof AAEdge) {
					int add = 0;
					if (protein.charAt(pos) == edges[j].getLetter()) {
						add = 1;
					}
					++pos;
					matchedEdges += add;
				} else {
					double needMass = masses[pos] + edges[j].getMass();
					int nextInd = masses.length - 1;
					int add = 0;
                    if (MassComparator.edgeMatches(masses[pos], masses[pos + 1], edges[j].getMass())) {
						nextInd = pos + 1;
						add = 1;
					} else if (masses.length - pos < 10) {
						for (int k = pos + 1; k < masses.length; ++k) {
							if (Math.abs(masses[k] - needMass) < Math.abs(masses[nextInd] - needMass)) {
								nextInd = k;
                                if (MassComparator.edgeMatches(masses[pos], masses[k], edges[j].getMass())) {
									add = 1;
									break;
								}
							}
						}
					} else {
						int l = pos + 1;
						int r = masses.length - 1;
						while (l <= r) {
							int med = (l + r) >> 1;
							if (masses[med] > needMass) {
								nextInd = med;
								r = med - 1;
							} else {
								l = med + 1;
							}
						}
						if (nextInd - 1 > pos && Math.abs(masses[nextInd - 1] - needMass) < Math.abs(masses[nextInd] - needMass)) {
							--nextInd;
						}
                        if (MassComparator.edgeMatches(masses[pos], masses[nextInd], edges[j].getMass())) {
							add = 1;
						}
					}
					pos = nextInd;
					matchedEdges += add;
				}
			}
			maxScore = Math.max(matchedEdges, maxScore);
		}
		return maxScore;
	}
}
