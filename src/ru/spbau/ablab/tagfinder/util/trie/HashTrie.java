package ru.spbau.ablab.tagfinder.util.trie;

import ru.spbau.ablab.tagfinder.Protein;
import ru.spbau.ablab.tagfinder.TagGenerator;
import ru.spbau.ablab.tagfinder.util.StringUtil;
import ru.spbau.ablab.tagfinder.util.pairs.Pair;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

public class HashTrie implements Trie {
    private HashMap<Long, ListEntry> matches = new HashMap<Long, ListEntry>();

    private static class ListEntry {
        private final ListEntry previous;
        private final Protein protein;
        private final double startMass;

        private ListEntry(ListEntry previous, Protein protein, double startMass) {
            this.previous = previous;
            this.protein = protein;
            this.startMass = startMass;
        }
    }

    public void add(long hashKey, Protein p, double startMass) {
        matches.put(hashKey, new ListEntry(matches.get(hashKey), p, startMass));
    }

    @Override
    public void add(char[] s, int l, int r, Protein p, double startMass) {
        long key = StringUtil.getHash(s, l, r);
        add(key, p, startMass);
    }

    @Override
    public void add(char[] s, Protein p, double startMass) {
        add(s, 0, s.length, p, startMass);
    }

    public Collection<Pair<Protein, Double>> getMatching(String s) {
        return getMatches(s.toCharArray());
    }

    @Override
    public Collection<Pair<Protein, Double>> getMatches(char[] s) {
        ArrayList<Pair<Protein, Double>> ans = new ArrayList<Pair<Protein, Double>>();
        for (ListEntry list : new ListEntry[]{matches.get(StringUtil.getHash(s)), matches.get(StringUtil.getRevHash(s))}) {
            if (list != null) {
                for (ListEntry entry = list; entry != null; entry = entry.previous) {
                    ans.add(new Pair<Protein, Double>(entry.protein, entry.startMass));
                }
            }
        }
        return ans;
    }
}
