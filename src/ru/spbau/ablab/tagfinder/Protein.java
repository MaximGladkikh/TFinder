package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.AAEdge;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.ArrayUtil;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.MassComparator;
import ru.spbau.ablab.tagfinder.util.StringUtil;

import java.util.Arrays;

import static ru.spbau.ablab.tagfinder.TagGenerator.AA_LET;
import static ru.spbau.ablab.tagfinder.TagGenerator.AA_MONO_MASS;

public class Protein {
	public static final double[] AA_MASS_ARRAY = new double[256];
	static {
		for (int i = 0; i < AA_LET.length; ++i) {
			AA_MASS_ARRAY[AA_LET[i]] = AA_MONO_MASS[i];
		}
	}

    private final String name;
	private final String protein;
	private final String revProtein;
    private Integer hashCode;
	public double[] masses;
    public final double parentMass;
    private double bestAbs;
    private static final double MAX_ACCEPTED_DISTANCE = ConfigReader.getDoubleProperty("MAX_ACCEPTED_DISTANCE");

    public Protein(String s, String name) {
		protein = s;
		masses = getMasses(s);
        parentMass = masses[masses.length - 1];
        this.name = name;
		revProtein = new StringBuilder(s).reverse().toString();
	}

    public String getString() {
        return protein;
    }

    public String getName() {
        return name;
    }

    @Override
    public int hashCode() {
        return hashCode == null ? hashCode = (int)StringUtil.getHash(name) : hashCode;
    }

    private double[] getMasses(String s) {
		double[] ans = new double[s.length() * 2 + 1];
		for (int i = 0; i < s.length(); ++i) {
			ans[i + 1] = ans[i] + AA_MASS_ARRAY[s.charAt(i)];
		}
        double sum = 0;
        for (int i = 0; i < s.length(); ++i) {
            ans[s.length() + 1 + i] = sum += AA_MASS_ARRAY[s.charAt(s.length() - 1 - i)];
        }
        Arrays.sort(ans);
		return ans;
	}

	public boolean contains(Path path) {
		boolean ans = getMaxMatch(path) == path.length();
		assert TagGenerator.EDGE_OF_TWO_AA || ans == (protein.contains(path.toString()) || revProtein.contains(path.toString()));
		return ans;
	}

	public int getMaxMatch(Path path) {
        bestAbs = Double.POSITIVE_INFINITY;
		return Math.max(getMaxMatch(path.getEdges(), masses, path.beginMass, path.length()), getMaxMatch(path.getReversedEdges(), masses, path.spectrum.parentMass - path.beginMass - path.getMass(), path.length()));
	}

	private int getMaxMatch(Edge[] edges, double[] masses, double beginMass, int length) {
		int maxScore = 0;
		for (int i = 0; i < masses.length - 1; ++i) {
            double abs = Math.abs(beginMass - masses[i]);
            if (abs > MAX_ACCEPTED_DISTANCE) {
                continue;
            }
			int pos = i;
			int matchedEdges = 0;
			for (int j = 0; j < edges.length && pos < protein.length(); ++j) {
				if (edges[j] instanceof AAEdge) {
					int add = 0;
					if (protein.charAt(pos) == edges[j].getLetter()) {
						add = 1;
					}
					++pos;
					matchedEdges += add;
				} else {
                    System.err.println(edges[j]);
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
            if (maxScore == length) {
                bestAbs = Math.min(bestAbs, abs);
            }
		}
		return maxScore;
	}

    public double getBestLastAlignment() {
        return bestAbs;
    }

    public int getMatchScore(Spectrum spectrum, double mult, double shift) {
        int ans = 0;
        for (Envelope envelope : spectrum.envelopes) {
            double needMass = envelope.getMass() * mult + shift;
            int index = ArrayUtil.getClosestIndex(masses, needMass);
            if (MassComparator.sameForDeletion(needMass, masses[index])) {
                ++ans;
            }
        }
        return ans;
    }
}
