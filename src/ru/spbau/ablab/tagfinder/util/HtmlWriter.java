package ru.spbau.ablab.tagfinder.util;

import static ru.spbau.ablab.tagfinder.StatisticsGenerator.MAX_PATHS;
import static ru.spbau.ablab.tagfinder.StatisticsGenerator.MAX_TAG_LENGTH;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

public class HtmlWriter extends PrintWriter {
    public HtmlWriter(String fileName) throws FileNotFoundException {
        super(fileName);
    }

    public void printHeader() {
        printThTaggedValue("scan id");
        printThTaggedValue("#peaks");
        printThTaggedValue("#tags");
        printOpenTh();
        printTaggedValue("div", "#proteins");
        println(" matched set");
        printCloseTh();
        printThTaggedValue("E-value");
        for (int i = 1; i <= MAX_PATHS; ++i) {
            printThTaggedValue("tag# " + i);
        }
    }

    public void printRatio(int[] found, int scansProcessed) {
        printOpenTag("tr");
        printThTaggedValue("% red");
        for (int i = 0; i < 4; ++i) {
            printEmptyTag("th");
        }
        for (int i = 0; i < MAX_PATHS; ++i) {
            printOpenTag("td");
            printTaggedValue("div", StringUtil.toStringPrecision(1. * found[i] / scansProcessed, 5), "align=center");
            printCloseTag("td");
        }
        printCloseTag("tr");
    }

    public void printFrequences(int[][] count) {
        for (int len = MAX_TAG_LENGTH; len >= 0; --len) {
            printOpenTag("tr");
            printThTaggedValue("length");
            printThTaggedValue("");
            printThTaggedValue("=");
            printThTaggedValue("");
            printThTaggedValue(len);
            for (int i = 0; i < MAX_PATHS; ++i) {
                printOpenTag("td");
                printTaggedValue("div", count[i][len], "align=center");
                printCloseTag("td");
            }
            printCloseTag("tr");
        }
    }

    public void printTagPrefix(String tag) {
        printTagPrefix(tag, "");
    }

    public void printTagPrefix(String tag, String attributes) {
        print("<" + tag + " " + attributes);
    }

    public void printEmptyTag(String tag) {
        printTagPrefix(tag);
        println("/>");
    }

    public void printTaggedValue(String tag, Object value) {
        printTaggedValue(tag, value, "");
    }

    public void printTaggedValue(String tag, Object value, String attributes) {
        printOpenTag(tag, attributes);
        println(value);
        printCloseTag(tag);
    }

    public void printOpenTag(String tag, String attributes) {
        print('<');
        printTagSuffix(tag, attributes);
    }

    private void printTagSuffix(String tag, String attributes) {
        print(tag);
        println(" nowrap " + attributes + ">");
    }

    public void printThTaggedValue(Object value) {
        printTaggedValue("th", value);
    }

    public void printOpenTh() {
        printOpenTag("th");
    }

    public void printOpenTag(String tag) {
        printOpenTag(tag, "");
    }

    public void printTagSuffix(String string) {
        printTagSuffix(string, "");
    }

    public void printCloseTh() {
        printCloseTag("th");
    }

    public void printCloseTag(String string) {
        print("</");
        printTagSuffix(string);
    }
}
