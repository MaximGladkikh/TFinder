package ru.spbau.ablab.tagfinder.path.edges;

import ru.spbau.ablab.tagfinder.util.MassComparator;

import static ru.spbau.ablab.tagfinder.Protein.AA_MASS_ARRAY;

public class AAEdge implements Edge {
    private final char c;
    private final char[][] decodings;

    public AAEdge(char c) {
        this.c = c;
        decodings = new char[][]{{c}};
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
    public char[][] getDecodings() {
        return decodings;
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
