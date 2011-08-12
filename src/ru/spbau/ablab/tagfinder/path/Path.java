package ru.spbau.ablab.tagfinder.path;

import java.util.Comparator;

import ru.spbau.ablab.tagfinder.StatisticsGenerator;
import ru.spbau.ablab.tagfinder.path.edges.AAEdge;
import ru.spbau.ablab.tagfinder.path.edges.Edge;

public class Path implements Comparable<Path> {
	public static final double SCORE_EPS = 1e-5;
	public static final int GAP_LENGTH_IN_AACOUNT = 1;
	
	public static final Comparator<Path> LENGTH_FIRST_COMPARATOR = new LengthFirstComparator();

	public final double score;
    public final Edge edge;
    public final int edges;
    private final int length;
    private final Path parent;

    public Path(Path parent, double score, Edge edge) {
        this.score = score;
        this.edge = edge;
        if (parent != null && parent.edge == null) {
            parent = null;
        }
        length = (parent == null ? 0 : parent.length()) + (edge instanceof AAEdge ? 1 : GAP_LENGTH_IN_AACOUNT) ;
        edges = 1 + (parent == null ? 0 : parent.edges);
        this.parent = parent;
    }

    private Path(Edge[] edges, int pos, double score) {
        this(pos <= 0 ? null : new Path(edges, pos - 1, score), score, pos >= 0 ? edges[pos] : null);
    }

	public Path(Edge[] edges) {
        this(edges, 0);
	}

	public Path(Edge[] edges, double score) {
        this(edges, edges.length - 1, score);
	}

	public int length() {
		return length;
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
        return new Path(this, this.score + Math.log(score), edge);
	}

    private static final Edge[] edgesBuffer = new Edge[StatisticsGenerator.MAX_TAG_LENGTH];

    public Edge[] getEdges() {
        Path path = this;
        int i = 0;
        while (path != null) {
            edgesBuffer[i++] = path.edge;
            path = path.parent;
        }
        Edge[] ans = new Edge[i];
        for (int j = 0; j < i; ++j) {
            ans[j] = edgesBuffer[i - j - 1];
        }
        return ans;
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
		return new Path(edges);
	}
}
