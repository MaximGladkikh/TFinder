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
            database = Database.getInstance();
            database.filter(false);
            String fileName = "unmatched.html";
            if (!ConfigReader.getBooleanProperty("USE_DEFAULT_FILENAME")) {
                fileName = "unmatched" + TagGenerator.MIN_TAG_LENGTH + "-" + TagGenerator.MAX_TAG_LENGTH + ".html";
            }
            HtmlWriter writer = new HtmlWriter(fileName);
            ArrayList<Integer> ids = database.getUnmatchedIds();
//            ids.retainAll(ArrayUtil.asList(671));
            UnmatchedScanProcessor processor = new UnmatchedScanProcessor();
            printStatistics(writer, ids, processor);
            writer.close();
            System.out.println("written successfully to " + fileName);
            System.out.println("done " + processor.getProcessedScansNumber());
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    protected class UnmatchedScanProcessor extends MatchedScanProcessor {
        private UnmatchedScanProcessor() {
            super(true);
        }

        protected void preProcess(int id, Collection<Path> paths) {
            bestFromAlign = database.getBestFromAlign(id, paths);
            if (bestFromAlign == null) {
                throw new AssertionError("");
            }
        }
    }
}
