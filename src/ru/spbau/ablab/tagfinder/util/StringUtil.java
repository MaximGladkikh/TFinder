package ru.spbau.ablab.tagfinder.util;

import java.util.ArrayList;
import java.util.StringTokenizer;

public class StringUtil {
    public static final int MUL = (int)1e9 + 7;

    private StringUtil() {
    }

    public static long getHash(String s) {
        return getHash(s.toCharArray());
    }

    public static String[] getTokenArray(String s) {
        StringTokenizer tokenizer = new StringTokenizer(s);
        ArrayList<String> ans = new ArrayList<String>();
        while (tokenizer.hasMoreTokens()) {
            ans.add(tokenizer.nextToken());
        }
        return ans.toArray(new String[ans.size()]);
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

    public static long getRevHash(char[] s) {
        long hashKey = 0;
        for (int i = s.length - 1; i >= 0; --i) {
            hashKey = hashKey * MUL + s[i];
        }
        return hashKey;
    }
}
