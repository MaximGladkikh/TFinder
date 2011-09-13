package ru.spbau.ablab.tagfinder.path;

import ru.spbau.ablab.tagfinder.TagGenerator;
import ru.spbau.ablab.tagfinder.path.edges.AAEdge;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.MassUtil;

import java.util.ArrayList;
import java.util.Comparator;

public class Path implements Comparable<Path> {
    public static final int GAP_LENGTH_IN_AA_COUNT = ConfigReader.getIntProperty("GAP_LENGTH_IN_AA_COUNT");
    public static final Comparator<Path> LENGTH_FIRST_COMPARATOR = new LengthFirstComparator();

    private static final double SCORE_EPS = 1e-5;

    public final double score;
    public final Edge edge;
    public final int edges;
    private final int length;
    private final Path parent;
    public final double beginMass;
    public final Spectrum spectrum;
    private Boolean isMonoTag = null;
    private final boolean reversed;

    public Path(Path parent, double score, Edge edge, double beginMass, Spectrum spectrum, boolean reversed) {
        this.reversed = reversed;
        this.spectrum = spectrum;
        this.score = score;
        this.edge = edge;
        if (parent != null && parent.edge == null) {
            parent = null;
        }
        length = (parent == null ? 0 : parent.length()) + (edge instanceof AAEdge ? 1 : GAP_LENGTH_IN_AA_COUNT);
        edges = 1 + (parent == null ? 0 : parent.edges);
        this.parent = parent;
        this.beginMass = beginMass;
//        assert Database.ALIGN || beginMass + getMass() <= spectrum.parentMass;
    }

    public Path(Path parent, double score, Edge edge, double beginMass, Spectrum spectrum) {
        this(parent, score, edge, beginMass, spectrum, false);
    }

    private Path(Edge[] edges, int pos, double score, double beginMass, Spectrum spectrum, boolean reversed) {
        this(pos <= 0 ? null : new Path(edges, pos - 1, score, beginMass, spectrum, reversed), score, pos >= 0 ? edges[pos] : null, beginMass, spectrum, reversed);
    }

    public Path(Edge[] edges, double beginMass, Spectrum spectrum, boolean reversed) {
        this(edges, 0, beginMass, spectrum, reversed);
    }

    public Path(Edge[] edges, double score, double beginMass, Spectrum spectrum, boolean reversed) {
        this(edges, edges.length - 1, score, beginMass, spectrum, reversed);
    }

    public int length() {
        return length;
    }

    public Path getReversed() {
        return new Path(getReversedEdges(), score, MassUtil.convertIonsType(beginMass + getMass(), spectrum.parentMass), spectrum, !reversed);
    }

    public double getMass() {
        double mass = 0;
        for (Path path = this; path != null && path.edge != null; path = path.parent) {
            mass += path.edge.getMass();
        }
        return mass;
    }

    public ArrayList<Double> getAllPossibleShifts() {
        ArrayList<Double> ans = new ArrayList<Double>();
        for (Envelope envelope : spectrum.getEnvelopes()) {
//            if (toString().equals("VFS") && Math.abs(envelope.getMass() - 4808) > 5) {
//                continue;
//            }
            double mass = envelope.getMass();
            if (checkDirection(spectrum, mass, getEdges(), 1)) {
                ans.add(envelope.getMass());
            }
            if (checkDirection(spectrum, mass, getReversedEdges(), -1)) {
                ans.add(-envelope.getMass());
            }
        }
        return  ans;
    }

    private boolean checkDirection(Spectrum spectrum, double mass, Edge[] edges, double direction) {
        for (Edge edge : edges)  {
            double needMass = mass + edge.getMass() * direction;
            Envelope nextEnvelope = spectrum.getClosest(needMass);
            if (MassUtil.compare(needMass, nextEnvelope.getMass(), 40e-6 * needMass) != 0) {
//                System.err.println(edge + " : " + mass + " " + needMass + " " + nextEnvelope.getMass());
                return false;
            }
            mass = nextEnvelope.getMass();
        }
        return true;
    }

    public boolean isMonoTag() {
        if (isMonoTag != null) {
            return isMonoTag;
        }
        Edge[] edges = getEdges();
        Edge[] reversedEdges = getReversedEdges();
        isMonoTag = spectrumHasPeaks(beginMass, edges, 1.) || spectrumHasPeaks(MassUtil.convertIonsType(beginMass + getMass(), spectrum.parentMass), edges, -1.);
        isMonoTag |= spectrumHasPeaks(beginMass, reversedEdges, 1.) || spectrumHasPeaks(MassUtil.convertIonsType(beginMass + getMass(), spectrum.parentMass), reversedEdges, -1.);
//        assert length() <= 1 || (TagGenerator.DOUBLE_MASSES || isMonoTag);
        return isMonoTag;
    }

    private boolean spectrumHasPeaks(double mass, Edge[] edges, double direction) {
        Envelope envelope = spectrum.getClosest(mass);
        if (MassUtil.compare(envelope.getMass(), mass, MassUtil.ERROR_THRESHOLD * 4) != 0) {
            return false;
        }
        for (Edge edge : edges) {
            mass += edge.getMass() * direction;
            envelope = spectrum.getClosest(mass);
            if (MassUtil.compare(mass, envelope.getMass(), MassUtil.ERROR_THRESHOLD * 4) != 0) {
                return false;
            }
            mass = envelope.getMass();
        }
        return true;
    }

    public boolean isReversed() {
        return reversed;
    }

    public boolean canBeReversedTo(Path path) {
        if (edges != path.edges) {
            return false;
        }
        Edge[] thisEdges = getEdges();
        Edge[] thatEdges = path.getEdges();
        boolean matches = true;
        boolean reversedMatches = true;
        for (int i = 0; i < edges; ++i) {
            matches &= thisEdges[i].compareTo(thatEdges[i]) == 0;
            reversedMatches &= thisEdges[i].compareTo(thatEdges[edges - i - 1]) == 0;
        }
        return matches | reversedMatches;
    }

    private static final class LengthFirstComparator implements Comparator<Path> {
        @Override
        public int compare(Path o1, Path o2) {
            int lengths = o1.compareLengths(o2);
            if (lengths != 0) {
                return -lengths;
            }
            return o1.lexicographicCompare(o2);
        }

    }

    @Override
    public int compareTo(Path o) {
        if (TagGenerator.SCORE_BY_LENGTH) {
            return LENGTH_FIRST_COMPARATOR.compare(this, o);
        }
        if (Math.abs(score - o.score) > SCORE_EPS) {
            return score > o.score ? -1 : 1;
        }
        return LENGTH_FIRST_COMPARATOR.compare(this, o);
    }

    private int lexicographicCompare(Path p2) {
        Path p1 = this;
        while (p1 != null && p2 != null) {
            int d = p1.edge.compareTo(p2.edge);
            if (d != 0) {
                return d;
            }
            p1 = p1.parent;
            p2 = p2.parent;
        }
        if (p1 == null) {
            return p2 == null ? 0 : -1;
        }
        return 1;
    }

    private int compareLengths(Path o) {
        if (length > o.length) {
            return -1;
        } else if (length < o.length) {
            return 1;
        }
        return 0;
    }

    public Path append(Edge edge, double score) {
        return new Path(this, Math.min(this.score, score), edge, beginMass, spectrum);
    }

    private static final Edge[] edgesBuffer = new Edge[TagGenerator.MAX_TAG_LENGTH * 2];

    public Edge[] getEdges() {
        int size = putPath();
        Edge[] ans = new Edge[size];
        for (int i = 0; i < size; ++i) {
            ans[i] = edgesBuffer[size - i - 1];
        }
        return ans;
    }

    public Edge[] getReversedEdges() {
        int size = putPath();
        Edge[] ans = new Edge[size];
        System.arraycopy(edgesBuffer, 0, ans, 0, size);
        return ans;
    }

    private int putPath() {
        Path path = this;
        int size = 0;
        while (path != null && path.edge != null) {
            edgesBuffer[size++] = path.edge;
            path = path.parent;
        }
        return size;
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        for (Edge edge : getEdges()) {
            builder.append(edge);
        }
        return builder.toString();
    }

    public Path subPath(int start, int end) {
        if (start > end) {
            throw new IllegalArgumentException();
        }
        Edge[] currentEdges = getEdges();
        Edge[] edges = new Edge[end - start];
        System.arraycopy(currentEdges, start, edges, 0, end - start);
        return new Path(edges, Double.NaN, spectrum, reversed);
    }
}