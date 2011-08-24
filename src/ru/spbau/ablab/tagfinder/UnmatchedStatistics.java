package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.Database;
import ru.spbau.ablab.tagfinder.util.io.HtmlWriter;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;

public class UnmatchedStatistics extends StatisticsGenerator {
    public static void main(String[] args) {
        new UnmatchedStatistics().run();
    }

    @Override
    public void run() {
        try {
            database = new Database();
            database.filter(false);
            String fileName = "unmatched.html";
            if (!ConfigReader.getBooleanProperty("USE_DEFAULT_FILENAME")) {
                fileName = "unmatched" + TagGenerator.MIN_TAG_LENGTH + "-" + TagGenerator.MAX_TAG_LENGTH + ".html";
            }
            HtmlWriter writer = new HtmlWriter(fileName);
            ArrayList<Integer> ids = database.getUnmatchedIds();
            UnmatchedScanProcessor processor = new UnmatchedScanProcessor();
            printStatistics(writer, ids, processor);
            writer.close();
            System.out.println("done " + processor.getProcessedScansNumber());
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    protected class UnmatchedScanProcessor extends MatchedScanProcessor {
        private UnmatchedScanProcessor() {
            super(true);
        }

        @Override
        protected void printTags(HtmlWriter writer, int[][] count, int[] found, Protein protein, Collection<Path> paths) {
            int pathN = 0;
            for (Path path : paths) {
                if (pathN == StatisticsGenerator.MAX_PATHS) {
                    break;
                }
                writer.printOpenTag("td");
                writer.printTaggedValue("div", path, "align=center");
                ++count[pathN][path.length()];
                writer.printCloseTag("span");
                writer.printCloseTag("div");
                writer.printCloseTag("td");
                ++pathN;
            }
            for (; pathN < MAX_PATHS; ++pathN) {
                writer.printEmptyTag("td");
            }
        }
    }
}
