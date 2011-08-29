package ru.spbau.ablab.tagfinder.util.io;

import ru.spbau.ablab.tagfinder.util.StringUtil;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

import static ru.spbau.ablab.tagfinder.StatisticsGenerator.MAX_PATHS;
import static ru.spbau.ablab.tagfinder.TagGenerator.MIN_TAG_LENGTH;

public class HtmlWriter extends PrintWriter {
    private final String name;

    public HtmlWriter(String fileName) throws FileNotFoundException {
        super(fileName);
        if (fileName.contains("/")) {
            fileName = fileName.substring(fileName.lastIndexOf("/"));
        }
        name = fileName.substring(0, fileName.indexOf('.')).replace('_', ' ');
        printOpenTag("html");
        printHeadTag();
        printOpenTag("body");
    }

    public void printHeader() {
        printOpenTag("tr");
        printThTaggedValue("scan id");
        printThTaggedValue("#peaks");
        printOpenTh();
        printTaggedValue("div", "%peaks");
        printTaggedValue("div", "used");
        printCloseTh();
        printOpenTh();
        printTaggedValue("div", "#proteins");
        printTaggedValue("div", "matched");
        printTaggedValue("div", "set");
        printCloseTh();
        printThTaggedValue("E-value");
        for (int i = 1; i <= MAX_PATHS; ++i) {
            printThTaggedValue("tag# " + i);
        }
        printOpenTh();
        printTaggedValue("div", "predicted");
        printTaggedValue("div", "protein");
        printTaggedValue("div", "matches set");
        printCloseTh();
        printThTaggedValue("matched protein");
        printCloseTag("tr");
    }

    public void close() {
        printCloseTag("body");
        printCloseTag("html");
        super.close();
    }

    private void printHeadTag() {
        printOpenTag("head");
        printTaggedValue("title", name);
        printCloseTag("head");
    }

    public void printRatio(int[] found, int scansProcessed, int predictedProteinFoundNumber, int matchedMsAlignNumber) {
        printOpenTag("tr");
        printThTaggedValue("% red");
        for (int i = 0; i < 3; ++i) {
            printEmptyTag("th");
        }
        double sum = 0;
        for (int i = 0; i < MAX_PATHS; ++i) {
            sum += 1. * found[i] / scansProcessed;
        }
        printThTaggedValue(String.format("%.2f", sum));
        for (int i = 0; i < MAX_PATHS; ++i) {
            printOpenTag("td");
            printTaggedValue("div", StringUtil.toStringPrecision(1. * found[i] / scansProcessed, 5), "align=center");
            printCloseTag("td");
        }
        printThTaggedValue(String.format("%.2f", predictedProteinFoundNumber * 1. / scansProcessed));
        printThTaggedValue(String.format("%.2f", matchedMsAlignNumber * 1. / scansProcessed));
        printCloseTag("tr");
    }

    public void printFrequencies(int[][] count) {
        int len = count[0].length - 1;
        while (len >= MIN_TAG_LENGTH) {
            boolean found = false;
            for (int i = 0; i < MAX_PATHS; ++i) {
                if (count[i][len] > 0) {
                    found = true;
                }
            }
            if (found) {
                break;
            }
            --len;
        }
        for (; len >= MIN_TAG_LENGTH; --len) {
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

    public void printMonoTags(int[] nTags, int[] monoTags) {
        printOpenTag("tr");
        printThTaggedValue("percentage");
        printThTaggedValue("of");
        printThTaggedValue("mono");
        printThTaggedValue("tags");
        int all = 0;
        int mono = 0;
        assert nTags.length == monoTags.length;
        for (int i = 0; i < nTags.length; ++i) {
            all += nTags[i];
            mono += monoTags[i];
        }
        printThTaggedValue(String.format("%.2f", 1. * mono / all));
        for (int i = 0; i < nTags.length; ++i) {
            printOpenTag("td");
            printTaggedValue("div", String.format("%.2f", 1. * monoTags[i] / nTags[i]), "align=center");
            printCloseTag("td");
        }
        printCloseTag("tr");
    }
}
