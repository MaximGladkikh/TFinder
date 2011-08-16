package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.util.Database;
import ru.spbau.ablab.tagfinder.util.StringUtil;
import ru.spbau.ablab.tagfinder.util.io.HtmlWriter;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class UnmatchedStatistics extends StatisticsGenerator {
    public static void main(String[] args) {
        new UnmatchedStatistics().run();
    }

    @Override
    public void run() {
        try {
            database = new Database(false);
            HtmlWriter writer = new HtmlWriter("unmatched.html");
            ArrayList<Integer> ids = database.getUnmatchedIds();
            int processed = printStatistics(writer, ids, new UnmatchedScanProcessor());
            writer.close();
            System.out.println("done " + processed);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    protected class UnmatchedScanProcessor implements ScanProcessor {
        public boolean processScan(int id, HtmlWriter writer, int[][] count, int[] found) {
            List<Path> paths = TagGenerator.getTopTags(database, id, MAX_PATHS);
            if (paths.isEmpty()) {
                return false;
            }
            writer.printOpenTag("tr");
            writer.printThTaggedValue(id);
            writer.printThTaggedValue(database.getSpectrum(id).envelopes.length);
            writer.printThTaggedValue(TagGenerator.getLastRatio());
            printProteinDBMatches(writer, paths);
            writer.printThTaggedValue(StringUtil.toStringScientific(database.getEValue(id), 2));
            printTags(writer, count, paths);
            writer.printTaggedValue("td", database.getLastMatchedProtein());
            writer.printCloseTag("tr");
            writer.flush();
            System.out.println(id + " ok");
            return true;
        }

        private void printTags(HtmlWriter writer, int[][] count, Collection<Path> paths) {
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
            for (;pathN < MAX_PATHS; ++pathN) {
                writer.printEmptyTag("td");
            }
        }
    }
}
