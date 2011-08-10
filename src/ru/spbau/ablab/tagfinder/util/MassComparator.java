package ru.spbau.ablab.tagfinder.util;

import java.util.Comparator;

import ru.spbau.ablab.tagfinder.StatisticsGenerator;

public class MassComparator implements Comparator<Double> {
    public static final MassComparator MASS_COMPARATOR = new MassComparator();
    public static final double FIRST_BY_EDGE_ERROR = ConfigReader.getDoubleProperty("FIRST_BY_EDGE_ERROR");
    public static final double ERROR_THRESHOLD = ConfigReader.getDoubleProperty("ERROR_THRESHOLD");
    private static final boolean RELATIVE_COMPARE = ConfigReader.getBooleanProperty("RELATIVE_COMPARE");
    private static final boolean REL_COMPARE_WITHOUT_MASS = ConfigReader.getBooleanProperty("REL_COMPARE_WITHOUT_MASS");

    private MassComparator() {
    }

    @Override
    public int compare(Double o1, Double o2) {
        return compare(o1.doubleValue(), o2.doubleValue());
    }

    public static int compare(double d1, double d2) {
        if (RELATIVE_COMPARE && (d1 + d2) * ERROR_THRESHOLD > Math.abs(d1 - d2)) {
            return 0;
        } else if (!RELATIVE_COMPARE && StatisticsGenerator.MASS_EPS > Math.abs(d1 - d2)) {
            return 0;
        }
        return d1 < d2 ? -1 : 1;
    }

    public static boolean edgeMatches(double m1, double m2, double edgeMass, Double parentMass) {
        if (parentMass == null) {
            return edgeMatches(m1, m2, edgeMass);
        }
        return Math.abs(Math.abs(m1 - m2) - edgeMass) <= FIRST_BY_EDGE_ERROR * parentMass;
    }

    public static boolean edgeMatches(double m1, double m2, double edgeMass) {
        if (REL_COMPARE_WITHOUT_MASS) {
            return Math.abs(Math.abs(m1 - m2) - edgeMass) <= edgeMass * ERROR_THRESHOLD;
        }
        return same(m1 + edgeMass, m2);
    }

    public static boolean same(double d1, double d2) {
        return compare(d1, d2) == 0;
    }
}
