package ru.spbau.ablab.tagfinder.path.edges;

import static ru.spbau.ablab.tagfinder.Protein.AA_MASS_ARRAY;
import ru.spbau.ablab.tagfinder.util.MassComparator;

public class AAEdge implements Edge {
	private char c;

	public AAEdge(char c) {
		this.c = c;
	}

	@Override
	public double getMass() {
		return AA_MASS_ARRAY[c];
	}

	@Override
	public Character getLetter() {
		return c;
	}

	@Override
	public int compareTo(Edge o) {
		if (o instanceof AAEdge) {
			return c - ((AAEdge) o).c;
		}
		return MassComparator.compare(getMass(), o.getMass());
	}

	@Override
	public String toString() {
		return "" + c;
	}
}
