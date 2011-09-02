package ru.spbau.ablab.tagfinder.path.edges;

import ru.spbau.ablab.tagfinder.util.MassComparator;

import java.util.ArrayList;

public class GapEdge implements Edge {
	private final double mass;

    private final ArrayList<String> decodingsList = new ArrayList<String>();
    private char[][] decodings;

	public GapEdge(double mass, String defaultDecoding) {
		this.mass = mass;
        decodingsList.add(defaultDecoding);
	}

    public void addDecoding(String decoding) {
        decodingsList.add(decoding);
    }

    public char[][] getDecodings() {
        if (decodings == null) {
            decodings = new char[decodingsList.size()][];
            for (int i = 0; i < decodingsList.size(); ++i) {
                decodings[i] = decodingsList.get(i).toCharArray();
            }
        }
        return decodings;
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
		return "[" + String.format("%.2f", mass) + "]";
	}
}
