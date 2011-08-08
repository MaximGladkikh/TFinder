package ru.spbau.ablab.tagfinder.path;

import ru.spbau.ablab.tagfinder.util.MassComparator;
import ru.spbau.ablab.tagfinder.util.StringUtil;


public class GapEdge implements Edge {
	private final double mass;

	public GapEdge(double mass) {
		this.mass = mass;
	}

	@Override
	public double getMass() {
		return mass;
	}

	@Override
	public Character getLetter() {
		return null;
	}

	@Override
	public int compareTo(Edge o) {
		return MassComparator.compare(mass, o.getMass());
	}
	
	@Override
	public String toString() {
		return "[" + StringUtil.toStringPrecision(mass, 2) + "]";
	}
}
