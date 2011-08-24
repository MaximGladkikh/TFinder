package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.Database;
import ru.spbau.ablab.tagfinder.util.StringUtil;
import ru.spbau.ablab.tagfinder.util.io.HtmlWriter;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class StatisticsGenerator implements Runnable {
    public static final boolean USE_DEFAULT_FILENAME = ConfigReader.getBooleanProperty("USE_DEFAULT_FILENAME");
    public static final boolean SCORE_BY_LENGTH = ConfigReader.getBooleanProperty("SCORE_BY_LENGTH");
    public static final double MASS_EPS = ConfigReader.getDoubleProperty("MASS_EPS");
    public static final int MAX_PATHS = ConfigReader.getIntProperty("MAX_PATHS");
    public static final String OUTPUT_FILE = ConfigReader.getProperty("OUTPUT_FILE");

    protected Database database;

    public static void main(String[] args) {
        new StatisticsGenerator().run();
    }

    @Override
    public void run() {
        ScanProcessor processor = new MatchedScanProcessor();
        try {
            database = new Database();
            ArrayList<Integer> scanIds = database.filter(true);
            HtmlWriter writer = new HtmlWriter(USE_DEFAULT_FILENAME ? OUTPUT_FILE : getOutputFilename());
            printStatistics(writer, scanIds, processor);
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(13);
        }
        System.out.println("done " + processor.getProcessedScansNumber());
    }

    protected void printStatistics(HtmlWriter writer, ArrayList<Integer> scanIds, ScanProcessor processor) {
        writer.printOpenTag("table", "cellpadding=0 cellspacing=20");
        writer.printHeader();
        int[] found = new int[MAX_PATHS];
        int[][] count = new int[MAX_PATHS][TagGenerator.MAX_TAG_LENGTH + 1];
        for (int id : scanIds) {
            processor.processScan(id, writer, count, found);
        }
        writer.printFrequencies(count);
        writer.printRatio(found, processor.getProcessedScansNumber(), processor.getPredictedProteinFoundNumber(), processor.getMatchedMsAlignNumber());
        writer.printCloseTag("table");
        System.out.println(processor.getMatchedMsAlignNumber() + " proteins matched MS-ALIGN+");
    }

    protected interface ScanProcessor {
        boolean processScan(int id, HtmlWriter writer, int[][] count, int[] found);
        int getMatchedMsAlignNumber();
        int getPredictedProteinFoundNumber();
        int getProcessedScansNumber();
    }

    protected class MatchedScanProcessor implements ScanProcessor {
        private final boolean skip;
        private final boolean dontMatch;
        public int matchedAlign = 0;
        public int predictedProteinFound = 0;
        public int processed = 0;

        public MatchedScanProcessor(boolean skip, boolean dontMatch) {
            this.skip = skip;
            this.dontMatch = dontMatch;
        }

        public MatchedScanProcessor(boolean skip) {
            this(skip, false);
        }

        public MatchedScanProcessor() {
            this(false);
        }

        public boolean processScan(int id, HtmlWriter writer, int[][] count, int[] found) {
            List<Path> paths = TagGenerator.getTopTags(database, id, MAX_PATHS);
            if (skip && paths.isEmpty()) {
                return false;
            }
            ++processed;
            writer.printOpenTag("tr");
            writer.printThTaggedValue(id);
            writer.printThTaggedValue(database.getSpectrum(id).envelopes.length);
            Protein protein = database.getProteinFromTable(id);
//            Protein protein = database.getProtein(id);
            writer.printThTaggedValue(TagGenerator.getLastRatio());
            if (dontMatch) {
                writer.printTaggedValue("td", "?");
            } else {
                printProteinDBMatches(writer, paths);
            }
            writer.printThTaggedValue(StringUtil.toStringScientific(database.getEValue(id), 2));
            printTags(writer, count, found, protein, paths);
            printMatchedProteinsStats(writer, id);
            writer.printCloseTag("tr");
            writer.flush();
            System.out.println(id + " ok");
            return true;
        }

        @Override
        public int getMatchedMsAlignNumber() {
            return matchedAlign;
        }

        @Override
        public int getPredictedProteinFoundNumber() {
            return predictedProteinFound;
        }

        @Override
        public int getProcessedScansNumber() {
            return processed;
        }

        protected void printTags(HtmlWriter writer, int[][] count, int[] found, Protein protein, Collection<Path> paths) {
            int pathN = 0;
            boolean foundFirst = false;
            for (Path path : paths) {
                if (pathN == MAX_PATHS) {
                    break;
                }
                writer.printOpenTag("td");
                writer.printOpenTag("div", "align=center");
                writer.printTagPrefix("span");
                boolean notPrint = false;
                Double bestAbs = null;
                if (protein != null && protein.contains(path)) {
                    bestAbs = protein.getBestLastAlignment();
                    if (!foundFirst) {
                        writer.printTagSuffix("style=\"color:red;font-weight:bold\"");
                        ++found[pathN];
                        foundFirst = true;
                    } else {
                        writer.printTagSuffix("style=\"color:magenta\"");
                    }
                } else if (protein != null) {
                    for (int mask = 1; mask < 4; ++mask) {
                        Path subPath = path.subPath(mask & 1, path.edges - mask / 2);
                        Edge[] edges = path.getEdges();
                        if (subPath.length() > 0 && protein.contains(subPath)) {
                            writer.printTagSuffix("style=\"color:blue\"");
                            if (mask % 2 == 1) {
                                writer.println(edges[0]);
                            }
                            writer.printCloseTag("span");
                            writer.printTaggedValue("span", subPath, "style=\"color:red\"");
                            if (mask / 2 == 1) {
                                writer.printTaggedValue("span", edges[path.edges - 1], "style=\"color:blue\"");
                            }
                            writer.printTaggedValue("span", protein.getMaxMatch(subPath), "style=\"color:blue\"");
                            notPrint = true;
                            break;
                        }
                    }
                    if (!notPrint) {
                        writer.printTagSuffix("style=\"color:blue\"");
                    }
                }
                ++count[pathN][path.length()];
                if (!notPrint) {
                    writer.println(path + " " + (protein == null ? 0 : protein.getMaxMatch(path)) + " " + (bestAbs == null ? "" : String.format("%.2f", bestAbs)));
                    writer.printCloseTag("span");
                }
                writer.printCloseTag("div");
                writer.printCloseTag("td");
                ++pathN;
            }
            for (; pathN < MAX_PATHS; ++pathN) {
                writer.printEmptyTag("td");
            }
        }

        protected void printMatchedProteinsStats(HtmlWriter writer, int id) {
            boolean predictedFound = database.getLastMatchedProteins().contains(database.getProteinFromTable(id).getName());
            writer.printTaggedValue("td", predictedFound ? "+" : "-");
            if (predictedFound) {
                ++predictedProteinFound;
            }
            boolean matchesAlign = database.getProteinFromTable(id).equals(database.getProtein(id));
            if (matchesAlign) {
                ++matchedAlign;
            }
            writer.printOpenTag("td");
            writer.printTaggedValue("div", database.getProteinName(id), "align=center" + (matchesAlign ? " style=\"color:red\"" : ""));
            writer.printCloseTag("td");
        }
    }


    protected void printProteinDBMatches(HtmlWriter writer, Collection<Path> paths) {
        writer.printThTaggedValue(database.getMatchedProteinsNumber(paths));
    }

    private String getOutputFilename() {
        return (ConfigReader.getBooleanProperty("ALIGN") ? "virt" : "exp") + "_" + TagGenerator.MIN_TAG_LENGTH + "-" + TagGenerator.MAX_TAG_LENGTH + ".html";
//        return (TagGenerator.EDGE_OF_TWO_AA ? "2" : "1") + "_" + (TagGenerator.DOUBLE_MASSES ? "doub" : "for") + "_" + ((int) (MassComparator.ERROR_THRESHOLD * 2e6)) + "ppm_" + (SCORE_BY_LENGTH ? "len" : "score") + "_" + (ConfigReader.getBooleanProperty("ALIGN") ? "virt" : "exp") + ".html";
    }
}