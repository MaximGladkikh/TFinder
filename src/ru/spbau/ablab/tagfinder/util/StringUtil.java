package ru.spbau.ablab.tagfinder.util;

public class StringUtil {
	private StringUtil() {
	}
	public static String toStringPrecision(double x, int afterMark) {
		String ans = Double.toString(x);
		int pos = ans.indexOf('.');
		if (pos >= 0) {
			ans = ans.substring(0, Math.min(ans.length(), pos + 1 + afterMark));
		}
		return ans;
	}
}
