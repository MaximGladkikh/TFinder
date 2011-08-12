package ru.spbau.ablab.tagfinder.util.pairs;

import ru.spbau.ablab.tagfinder.util.pairs.Pair;

public class ComparablePair<T1 extends Comparable<T1>, T2 extends Comparable<T2>> extends Pair<T1, T2> implements Comparable<Pair<T1, T2>> {
	public ComparablePair(T1 a, T2 b) {
		super(a, b);
	}

	@Override
	public int compareTo(Pair<T1, T2> o) {
		int d = a.compareTo(o.a);
		if (d == 0) {
			return b.compareTo(o.b);
		}
		return d;
	}
}
