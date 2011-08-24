package ru.spbau.ablab.tagfinder.util;

import java.io.PrintWriter;
import java.io.StringWriter;

public class StringUtil {
    public static final int MUL = (int)1e9 + 7;

    private StringUtil() {
    }

    public static long getHash(String s) {
        return getHash(s.toCharArray());
    }

    public static long getHash(char[] s) {
        long hashKey = 0;
        for (char c : s) {
            hashKey = hashKey * MUL + c;
        }
        return hashKey;
    }

    public static long getHash(char[] s, int l, int r) {
        long hashKey = 0;
        for (int i = l; i < r; ++i) {
            hashKey = hashKey * MUL + s[i];
        }
        return hashKey;
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

    public static long getRevHash(char[] s) {
        long hashKey = 0;
        for (int i = s.length - 1; i >= 0; --i) {
            hashKey = hashKey * MUL + s[i];
        }
        return hashKey;
    }
}
