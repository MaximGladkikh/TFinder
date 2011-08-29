package ru.spbau.ablab.tagfinder.path;

import ru.spbau.ablab.tagfinder.StatisticsGenerator;
import ru.spbau.ablab.tagfinder.TagGenerator;
import ru.spbau.ablab.tagfinder.path.edges.AAEdge;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.Database;

import java.util.Comparator;

public class Path implements Comparable<Path> {
	public static final double SCORE_EPS = 1e-5;
	public static final int GAP_LENGTH_IN_AACOUNT = 1;
	
	public static final Comparator<Path> LENGTH_FIRST_COMPARATOR = new LengthFirstComparator();

	public final double score;
    public final Edge edge;
    public final int edges;
    private final int length;
    private final Path parent;
    public final double beginMass;
    public final Spectrum spectrum;
    private Boolean monoTag = null;
    private final boolean reversed;

    public Path(Path parent, double score, Edge edge, double beginMass, Spectrum spectrum, boolean reversed) {
        this.reversed = reversed;
        this.spectrum = spectrum;
        if (beginMass < -3) {
            throw new AssertionError("bm = " + beginMass);
        }
        this.score = score;
        this.edge = edge;
        if (parent != null && parent.edge == null) {
            parent = null;
        }
        length = (parent == null ? 0 : parent.length()) + (edge instanceof AAEdge ? 1 : GAP_LENGTH_IN_AACOUNT) ;
        edges = 1 + (parent == null ? 0 : parent.edges);
        this.parent = parent;
        this.beginMass = beginMass;
    }

    public Path(Path parent, double score, Edge edge, double beginMass, Spectrum spectrum) {
        this(parent, score, edge, beginMass, spectrum, false);
    }

    private Path(Edge[] edges, int pos, double score, double beginMass, Spectrum spectrum, boolean reversed) {
        this(pos <= 0 ? null : new Path(edges, pos - 1, score, beginMass, spectrum, reversed), score, pos >= 0 ? edges[pos] : null, beginMass, spectrum, reversed);
    }

	public Path(Edge[] edges, double beginMass, Spectrum spectrum, boolean reversed) {
        this(edges, 0, beginMass, spectrum, reversed);
	}

	public Path(Edge[] edges, double score, double beginMass, Spectrum spectrum, boolean reversed) {
        this(edges, edges.length - 1, score, beginMass, spectrum, reversed);
	}

	public int length() {
		return length;
	}

    public Path getReversed() {
        return new Path(getReversedEdges(), score, spectrum.parentMass - beginMass - getMass() + Database.WATER_MASS, spectrum, !reversed);
    }

    public double getMass() {
        double mass = 0;
        Path path = this;
        while (path != null && path.edge != null) {
            mass += path.edge.getMass();
            path = path.parent;
        }
        return mass;
    }

    public boolean isMonoTag() {
        if (monoTag != null) {
            return monoTag;
        }
        Edge[] edges = getEdges();
        return monoTag = spectrumHasPeaks(beginMass, edges, 1.) || spectrumHasPeaks(spectrum.parentMass - beginMass + Database.WATER_MASS, edges, -1.);
    }

    private boolean spectrumHasPeaks(double mass, Edge[] edges, double d) {
        boolean matched = spectrum.hasPeak(mass);
        for (Edge edge : edges) {
            mass += edge.getMass() * d;
            matched &= spectrum.hasPeak(mass);
        }
        return matched;
    }

    public boolean isReversed() {
        return reversed;
    }

    private static final class LengthFirstComparator implements Comparator<Path>{
		@Override
		public int compare(Path o1, Path o2) {
			int lengths = o1.compareLengths(o2);
			if (lengths != 0) {
				return lengths;
			}
			return o1.lexicographicCompare(o2);
		}
		
	}

	@Override
	public int compareTo(Path o) {
		if (StatisticsGenerator.SCORE_BY_LENGTH) {
			return LENGTH_FIRST_COMPARATOR.compare(this, o);
		}
		if (Math.abs(score - o.score) > SCORE_EPS) {
			return score > o.score ? -1 : 1;
		}
		return LENGTH_FIRST_COMPARATOR.compare(this, o);
	}

	private int lexicographicCompare(Path p2) {
        Path p1 = this;
        while (p1 != null && p2 != null) {
            int d = p1.edge.compareTo(p2.edge);
            if (d != 0) {
                return d;
            }
            p1 = p1.parent;
            p2 = p2.parent;
        }
        if (p1 == null) {
            return p2 == null ? 0 : -1;
        }
		return 1;
	}

	private int compareLengths(Path o) {
		if (length > o.length) {
			return -1;
		} else if (length < o.length) {
			return 1;
		}
		return 0;
	}

	public Path append(Edge edge, double score) {
        return new Path(this, this.score + Math.log(score), edge, beginMass, spectrum);
	}

    private static final Edge[] edgesBuffer = new Edge[TagGenerator.MAX_TAG_LENGTH * 2];

    public Edge[] getEdges() {
        int size = putPath();
        Edge[] ans = new Edge[size];
        for (int i = 0; i < size; ++i) {
            ans[i] = edgesBuffer[size - i - 1];
        }
        return ans;
    }

    public Edge[] getReversedEdges() {
        int size = putPath();
        Edge[] ans = new Edge[size];
        System.arraycopy(edgesBuffer, 0, ans, 0, size);
        return ans;
    }

    private int putPath() {
        Path path = this;
        int size = 0;
        while (path != null && path.edge != null) {
            edgesBuffer[size++] = path.edge;
            path = path.parent;
        }
        return size;
    }

    @Override
	public String toString() {
        StringBuilder builder = new StringBuilder();
		for (Edge edge : getEdges()) {
			builder.append(edge);
		}
		return builder.toString();
	}

	public Path subPath(int start, int end) {
		if (start > end) {
			throw new IllegalArgumentException();
		}
        Edge[] currentEdges = getEdges();
		Edge[] edges = new Edge[end - start];
        System.arraycopy(currentEdges, start, edges, 0, end - start);
		return new Path(edges, Double.NaN, spectrum, reversed);
	}
}
