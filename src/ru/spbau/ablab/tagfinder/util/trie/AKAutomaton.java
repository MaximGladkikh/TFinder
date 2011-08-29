package ru.spbau.ablab.tagfinder.util.trie;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.Database;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;

import static ru.spbau.ablab.tagfinder.TagGenerator.*;

public class AKAutomaton {
    private static final double MAX_DISTANCE = ConfigReader.getDoubleProperty("MAX_DISTANCE");
    private static final short ROOT = 0;
    private static final int[] AA_INDEX = new int[Character.MAX_VALUE];

    static {
        for (int i = 0; i < ALPHABET_SIZE; ++i) {
            AA_INDEX[AA_LET[i]] = i;
        }
    }

    private int[][] next;
    private Path[] accepted;
    private int[] sufLinks;
    private int verticesCount = 1;

    private AKAutomaton() {
    }

    public AKAutomaton(Collection<Path> paths) {
        next = new int[1][AA_LET.length];
        Arrays.fill(next[ROOT], -1);
        accepted = new Path[1];
        for (Path path : paths) {
            addPath(path);
            addPath(path.getReversed());
        }
        calcLinks();
    }

    public int getMatchesNumber(String string, double parentMass) {
        int vertex = 0;
        double prefixMass = 0;
        HashSet<Path> paths = new HashSet<Path>();
        for (int i = 0; i < string.length(); ++i) {
            int index = AA_INDEX[string.charAt(i)];
            prefixMass += AA_MONO_MASS[index];
            vertex = next[vertex][index];
            for (int link = vertex; accepted[link] != null; link = sufLinks[link]) {
                Path tag = accepted[link];
                if (!tag.isReversed()) {
                    double position = prefixMass - tag.getMass();
                    if (Math.abs(position - tag.beginMass) <= MAX_DISTANCE) {
                        paths.add(tag);
                    }
                } else {
                    double position = parentMass - prefixMass + Database.WATER_MASS;
                    if (Math.abs(position - tag.beginMass) <= MAX_DISTANCE) {
                        paths.add(tag);
                    }
                }
            }
        }
        return paths.size();
    }

    public boolean acceptsString(String s) {
        int vertex = 0;
        for (int i = 0; i < s.length(); ++i) {
            int c = AA_INDEX[s.charAt(i)];
            vertex = next[vertex][c];
            if (accepted[vertex] != null) {
                return true;
            }
        }
        return false;
    }

    private void calcLinks() {
        sufLinks = new int[next.length];
        calcLinks(ROOT, new StringBuilder());
    }

    private int getBestAlignmentVertexIndex(StringBuilder stringBuilder) {
        for (int l = 1; l < stringBuilder.length(); ++l) {
            int index = getIndex(stringBuilder.substring(l, stringBuilder.length()).toCharArray());
            if (index >= 0) {
                return index;
            }
        }
        return 0;
    }

    private void calcLinks(int vertex, StringBuilder stringBuilder) {
        sufLinks[vertex] = getBestAlignmentVertexIndex(stringBuilder);
        if (accepted[vertex] == null) {
            accepted[vertex] = accepted[sufLinks[vertex]];
        }
        for (int i = 0; i < ALPHABET_SIZE; ++i) {
            stringBuilder.append(AA_LET[i]);
            if (next[vertex][i] >= 0) {
                calcLinks(next[vertex][i], stringBuilder);
            } else {
                next[vertex][i] = getBestAlignmentVertexIndex(stringBuilder);
            }
            stringBuilder.replace(stringBuilder.length() - 1, stringBuilder.length(), "");
        }
    }

    private int getIndex(char[] s) {
        int vertex = ROOT;
        for (char c : s) {
            vertex = next[vertex][AA_INDEX[c]];
            if (vertex < 0) {
                return -1;
            }
        }
        return vertex;
    }

    private void addPath(Path path) {
        addPath(ROOT, path, path.getEdges(), 0);
    }

    private void addPath(int vertex, Path path, Edge[] edges, int pos) {
        if (edges.length == pos) {
            accepted[vertex] = path;
            return;
        }
        for (char[] decoding : edges[pos].getDecodings()) {
            int nextVertex = vertex;
            for (char c : decoding) {
                nextVertex = getNextVertex(nextVertex, c);
            }
            addPath(nextVertex, path, edges, pos + 1);
        }
    }

    private int getNextVertex(int vertex, int c) {
        c = AA_INDEX[c];
        if (next[vertex][c] < 0) {
            next[vertex][c] = getNewVertex();
        }
        return next[vertex][c];
    }

    public int getNewVertex() {
        if (verticesCount == next.length) {
            int[][] next2 = new int[next.length * 2][];
            System.arraycopy(next, 0, next2, 0, next.length);
            next = next2;
            Path[] accepted2 = new Path[accepted.length * 2];
            System.arraycopy(accepted, 0, accepted2, 0, accepted.length);
            accepted = accepted2;
        }
        next[verticesCount] = new int[ALPHABET_SIZE];
        Arrays.fill(next[verticesCount], -1);
        return verticesCount++;
    }
}
