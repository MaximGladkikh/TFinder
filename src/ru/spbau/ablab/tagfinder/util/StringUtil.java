package ru.spbau.ablab.tagfinder.util;

import java.io.PrintWriter;
import java.io.StringWriter;

public class StringUtil {
    private StringUtil() {
    }

    private static StringWriter stringWriter = new StringWriter();
    private static PrintWriter printWriter = new PrintWriter(stringWriter);


    public static String toStringScientific(double x, int afterMark) {
        printWriter.printf("%." + afterMark + "E", x);
        try {
            return stringWriter.toString();
        } finally {
            printWriter = new PrintWriter((stringWriter = new StringWriter()));
        }
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
