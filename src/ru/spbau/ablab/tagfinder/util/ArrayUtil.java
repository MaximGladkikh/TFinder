package ru.spbau.ablab.tagfinder.util;

import java.util.ArrayList;
import java.util.List;

public class ArrayUtil {
    private ArrayUtil() {
    }

    public static <T> List<T> asList(T... o) {
        ArrayList<T> ans = new ArrayList<T>();
        for (T t : o.clone()) {
            if (t != null) {
                ans.add(t);
            }
        }
        return ans;
    }

    public static int getClosestIndex(double[] a, double key) {
        int l = 0;
        int r = a.length - 1;
        int ans = r;
        while (l <= r) {
            int med = (l + r) / 2;
            if (key <= a[med]) {
                r = med - 1;
                ans = med;
            } else {
                l = med + 1;
            }
        }
        if (ans > 0 && Math.abs(a[ans] - key) > Math.abs(a[ans - 1] - key)) {
            --ans;
        }
        if (ans + 1 < a.length && Math.abs(a[ans] - key) > Math.abs(a[ans + 1] - key)) {
            ++ans;
        }
        return ans;
    }
}
