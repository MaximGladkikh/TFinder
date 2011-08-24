package ru.spbau.ablab.tagfinder.util.trie;

import ru.spbau.ablab.tagfinder.Protein;
import ru.spbau.ablab.tagfinder.util.pairs.Pair;

import java.util.Collection;

public interface Trie {
    void add(char[] s, Protein p, double start);
    void add(char[] s, int l, int r, Protein p, double start);
    Collection<Pair<Protein, Double>> getMatching(char[] s);
    Collection<Pair<Protein, Double>> getMatching(String s);
}
