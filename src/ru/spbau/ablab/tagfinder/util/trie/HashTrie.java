package ru.spbau.ablab.tagfinder.util.trie;

import ru.spbau.ablab.tagfinder.Protein;
import ru.spbau.ablab.tagfinder.util.StringUtil;
import ru.spbau.ablab.tagfinder.util.pairs.Pair;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

public class HashTrie implements Trie {
    private HashMap<Long, HashSet<Pair<Protein, Double>>> matches = new HashMap<Long, HashSet<Pair<Protein, Double>>>();

    @Override
    public void add(char [] s, int l, int r, Protein p, double startMass) {
        long key = StringUtil.getHash(s, l, r);
        HashSet<Pair<Protein, Double>> set = matches.get(key);
        if (set == null) {
            matches.put(key, set = new HashSet<Pair<Protein, Double>>());
        }
        set.add(new Pair<Protein, Double>(p, startMass));
    }

    @Override
    public void add(char[] s, Protein p,double startMass) {
        add(s, 0, s.length, p, startMass);
    }

    public Collection<Pair<Protein, Double>> getMatching(String s) {
        return getMatching(s.toCharArray());
    }

    @Override
    public Collection<Pair<Protein,Double>> getMatching(char[] s) {
        Collection<Pair<Protein, Double>> match = matches.get(StringUtil.getHash(s));
        if (match == null) {
            match = new ArrayList<Pair<Protein,Double>>(0);
        }
        Collection<Pair<Protein, Double>> revMatch = matches.get(StringUtil.getRevHash(s));
        if (revMatch != null) {
            for (Pair<Protein, Double> pair : revMatch) {
                match.add(pair);
            }
        }
        return match;
    }
}
