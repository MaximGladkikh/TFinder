package ru.spbau.ablab.tagfinder.path;

import java.util.Comparator;

import ru.spbau.ablab.tagfinder.StatisticsGenerator;

public class Path implements Comparable<Path> {
	public static final double SCORE_EPS = 1e-5;
	public static final int GAP_LENGTH_IN_AACOUNT = 1;
	
	public static final Comparator<Path> LENGTH_FIRST_COMPARATOR = new LengthFirstComparator();

	public Edge[] edges;
	public double score;

	public Path(String path, double score) {
		this.edges = new Edge[path.length()];
		for (int i = 0; i < path.length(); ++i) {
			this.edges[i] = new AAEdge(path.charAt(i));
		}
		this.score = score;
	}

	public Path(Edge[] edges) {
		this.edges = edges;
	}

	public Path(Edge[] edges, double score) {
		this.edges = edges;
		this.score = score;
	}

	public int length() {
		int length = 0;
		for (Edge edge : edges) {
			if (edge instanceof AAEdge) {
				++length;
			} else {
				length += GAP_LENGTH_IN_AACOUNT;
			}
		}
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

	private int lexicographicCompare(Path o) {
		for (int i = 0; i < edges.length && i < o.edges.length; ++i) {
			int d = edges[i].compareTo(o.edges[i]);
			if (d != 0) {
				return d;
			}
		}
		return 0;
	}

	private int compareLengths(Path o) {
		if (edges.length > o.edges.length) {
			return -1;
		} else if (edges.length < o.edges.length) {
			return 1;
		}
		return 0;
	}

	public Path append(Edge edge, double score) {
		Edge[] edges = new Edge[this.edges.length + 1];
		System.arraycopy(this.edges, 0, edges, 0, this.edges.length);
		edges[this.edges.length] = edge;
//		return new Path(edges, Math.min(this.score, Math.cos(Math.log(score))));
		return new Path(edges, this.score + Math.log(score));
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		for (Edge edge : edges) {
			builder.append(edge);
		}
		return builder.toString();
	}

	public Path subPath(int start, int end) {
		if (start > end) {
			throw new IllegalArgumentException();
		}
		Edge[] edges = new Edge[end - start];
        System.arraycopy(this.edges, start, edges, 0, end - start);
		return new Path(edges);
	}

	// private int getMaxMatch(String protein, String path) {
	// int ans = 0;
	// for (int i = 0; i + path.length() <= protein.length(); ++i) {
	// int count = 0;
	// for (int j = 0; j < path.length(); ++j) {
	// if (protein.charAt(i + j) == path.charAt(j)) {
	// ++count;
	// }
	// }
	// ans = Math.max(ans, count);
	// }
	// return ans;
	// }
}
