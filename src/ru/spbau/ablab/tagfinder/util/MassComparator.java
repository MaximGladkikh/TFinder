package ru.spbau.ablab.tagfinder.util;

import java.util.Comparator;

import ru.spbau.ablab.tagfinder.StatisticsGeneratorExperimental;

public class MassComparator implements Comparator<Double> {
	public static final MassComparator MASS_COMPARATOR = new MassComparator();
	public static final double ERROR_THRESHOLD = ConfigReader.getDoubleProperty("ERROR_THRESHOLD");
	private static final boolean RELATIVE_COMPARE = ConfigReader.getBoolProperty("RELATIVE_COMPARE");
	
	private MassComparator() {
	}
	
	@Override
	public int compare(Double o1, Double o2) {
		return compare(o1.doubleValue(), o2.doubleValue());
	}
	
	public static int compare(double d1, double d2) {
		if (RELATIVE_COMPARE && d1 * ERROR_THRESHOLD > Math.abs(d1 - d2)) {
			return 0;
		} else if (!RELATIVE_COMPARE && StatisticsGeneratorExperimental.MASS_EPS > Math.abs(d1 - d2)) {
			return 0;
		}
		return d1 < d2 ? -1 : 1;
	}
	
	public static boolean same(double d1, double d2) {
		return compare(d1, d2) == 0;
	}
}
