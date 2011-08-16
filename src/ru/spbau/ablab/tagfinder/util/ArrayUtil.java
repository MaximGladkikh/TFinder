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
}
