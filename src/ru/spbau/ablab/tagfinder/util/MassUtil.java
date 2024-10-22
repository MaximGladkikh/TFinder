package ru.spbau.ablab.tagfinder.util;

import java.util.Comparator;

public class MassUtil implements Comparator<Double> {
    public static final double FIRST_BY_EDGE_ERROR = ConfigReader.getDoubleProperty("FIRST_BY_EDGE_ERROR");
    public static final double ERROR_THRESHOLD = ConfigReader.getDoubleProperty("ERROR_THRESHOLD");
    public static final Comparator<Double> MASS_COMPARATOR = new MassUtil();
    private static final boolean RELATIVE_COMPARE = ConfigReader.getBooleanProperty("RELATIVE_COMPARE");
    private static final boolean REL_COMPARE_WITHOUT_MASS = ConfigReader.getBooleanProperty("REL_COMPARE_WITHOUT_MASS");
    public static final double MASS_EPS = ConfigReader.getDoubleProperty("MASS_EPS");

    private MassUtil() {
    }

    @Override
    public int compare(Double o1, Double o2) {
        return compare(o1.doubleValue(), o2.doubleValue());
    }

    public static int compare(double d1, double d2, double errorThreshold) {
        if (RELATIVE_COMPARE && (d1 + d2) * errorThreshold > Math.abs(d1 - d2)) {
            return 0;
        } else if (!RELATIVE_COMPARE && MASS_EPS > Math.abs(d1 - d2)) {
            return 0;
        }
        return d1 < d2 ? -1 : 1;

    }

    public static boolean same(double d1, double d2) {
        return compare(d1, d2, ERROR_THRESHOLD) == 0;
    }

    public static int compare(double d1, double d2) {
        return compare(d1, d2, ERROR_THRESHOLD);
    }

    public static boolean edgeMatches(double m1, double m2, double edgeMass, double parentMass, Double massCorrection) {
        if (massCorrection == null) {
            return Math.abs(Math.abs(m1 - m2) - edgeMass) <= FIRST_BY_EDGE_ERROR * parentMass;
        }
        return edgeMatches(m1, m2, edgeMass);
    }

    public static boolean edgeMatches(double m1, double m2, double edgeMass) {
        double difference = Math.abs(Math.abs(m1 - m2) - edgeMass);
        if (REL_COMPARE_WITHOUT_MASS) {
            return difference <= edgeMass * ERROR_THRESHOLD;
        }
        return difference <= ERROR_THRESHOLD * (m1 + m2);
    }

    public static double convertIonsType(double mass, double parentMass) {
        return parentMass - mass;
    }
}
