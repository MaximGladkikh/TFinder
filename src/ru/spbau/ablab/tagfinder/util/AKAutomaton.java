package ru.spbau.ablab.tagfinder.util;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import static ru.spbau.ablab.tagfinder.database.Database.*;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;

public class AKAutomaton {
    private static final double MAX_DISTANCE = ConfigReader.getDoubleProperty("MAX_DISTANCE");
    private static final short ROOT = 0;

    private int[][] next;
    private Path[] acceptedPath;
    private int[] sufLinks;
    private int[] previousAcceptedLink;
    private int verticesCount = 1;

    private AKAutomaton() {
    }

    public AKAutomaton(Collection<Path> paths) {
        next = new int[1][AA_LET.length];
        Arrays.fill(next[ROOT], -1);
        acceptedPath = new Path[1];
        for (Path path : paths) {
            addPath(path);
            addPath(path.getReversed());
        }
        previousAcceptedLink = new int[next.length];
        Arrays.fill(previousAcceptedLink, -1);
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
            for (int link = acceptedPath[vertex] == null ? previousAcceptedLink[vertex] : vertex; link > 0; link = previousAcceptedLink[link]) {
                Path tag = acceptedPath[link];
                if (!tag.isReversed()) {
                    double position = prefixMass - tag.getMass();
                    if (Math.abs(position - tag.beginMass) <= MAX_DISTANCE) {
                        paths.add(tag);
                    }
                } else {
                    double position = MassUtil.convertIonsType(prefixMass, parentMass);
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
            if (acceptedPath[vertex] != null || previousAcceptedLink[vertex] >= 0) {
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
        int link = sufLinks[vertex] = getBestAlignmentVertexIndex(stringBuilder);
        if (link >= 0) {
            previousAcceptedLink[vertex] = acceptedPath[link] == null ? previousAcceptedLink[link] : link;
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
            acceptedPath[vertex] = path;
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
            Path[] accepted2 = new Path[acceptedPath.length * 2];
            System.arraycopy(acceptedPath, 0, accepted2, 0, acceptedPath.length);
            acceptedPath = accepted2;
        }
        next[verticesCount] = new int[ALPHABET_SIZE];
        Arrays.fill(next[verticesCount], -1);
        return verticesCount++;
    }
}

