package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.database.Database;
import ru.spbau.ablab.tagfinder.path.Path;
import ru.spbau.ablab.tagfinder.path.edges.Edge;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.MassUtil;
import ru.spbau.ablab.tagfinder.util.io.HtmlWriter;
import ru.spbau.ablab.tagfinder.util.pairs.Pair;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import static ru.spbau.ablab.tagfinder.TagGenerator.*;

public class StatisticsGenerator implements Runnable {
    public static final boolean USE_DEFAULT_OUTPUT_FILENAME = ConfigReader.getBooleanProperty("USE_DEFAULT_OUTPUT_FILENAME");
    public static final int MAX_TAGS_IN_SET = ConfigReader.getIntProperty("MAX_TAGS_IN_SET");
    public static final String OUTPUT_FILE = ConfigReader.getProperty("OUTPUT_FILE");
    public static final boolean DISABLE_DB_TAG_SEARCH = ConfigReader.getBooleanProperty("DISABLE_DB_TAG_SEARCH");

    protected Database database;

    public static void main(String[] args) {
        new StatisticsGenerator().run();
    }

    @Override
    public void run() {
        ScanProcessor processor = new MatchedScanProcessor(false, DISABLE_DB_TAG_SEARCH);
        try {
            database = Database.getInstance();
            ArrayList<Integer> scanIds = database.filter(true);
            String fileName = USE_DEFAULT_OUTPUT_FILENAME ? OUTPUT_FILE : getOutputFilename();
            HtmlWriter writer = new HtmlWriter(fileName);
            printStatistics(writer, scanIds, processor);
            writer.close();
            System.out.println("Written successfully to " + fileName);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(13);
        }
        System.out.println("done " + processor.getProcessedScansNumber());
    }

    protected void printStatistics(HtmlWriter writer, ArrayList<Integer> scanIds, ScanProcessor processor) {
        long startTime = System.currentTimeMillis();
        writer.printOpenTag("table", "cellpadding=0 cellspacing=20");
        writer.printHeader();
        int[] found = new int[MAX_TAGS_IN_SET];
        int[][] count = new int[MAX_TAGS_IN_SET][TagGenerator.MAX_TAG_LENGTH * 2];
        int[] nTags = new int[MAX_TAGS_IN_SET];
        int[] monoTags = new int[MAX_TAGS_IN_SET];
        for (int id : scanIds) {
            processor.processScan(id, writer, count, found, nTags, monoTags);
        }
        writer.printFrequencies(count);
        writer.printMonoTags(nTags, monoTags);
        writer.printRatio(found, processor.getProcessedScansNumber(), processor.getPredictedProteinFoundNumber(), processor.getMatchedMsAlignNumber());
        writer.printCloseTag("table");
        System.out.println(processor.getMatchedMsAlignNumber() + " proteins matched MS-USE_VIRTUAL_SPECTRA+");
        System.out.printf("Statistics generated in %.3fms\n", 1e-3 * (System.currentTimeMillis() - startTime));
    }

    protected interface ScanProcessor {
        boolean processScan(int id, HtmlWriter writer, int[][] count, int[] found, int[] nTags, int[] monoTags);

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
        protected Pair<Double, Protein> bestFromAlign;
        protected boolean printShifts;
        private boolean dontSearchForBest = ConfigReader.getBooleanProperty("DISABLE_PROTEIN_IDENTIFICATION_BY_TAGS");

        public MatchedScanProcessor(boolean skip, boolean dontMatch) {
            this.skip = skip;
            this.dontMatch = dontMatch;
            printShifts = Database.USE_VIRTUAL_SPECTRA;
        }

        public MatchedScanProcessor(boolean skip) {
            this(skip, false);
        }

        public MatchedScanProcessor() {
            this(false);
        }

        protected void preProcess(int id, Collection<Path> paths) {
        }

        public boolean processScan(int id, HtmlWriter writer, int[][] count, int[] found, int[] nTags, int[] monoTags) {
            List<Path> paths = TagGenerator.getTopTags(id, MAX_TAGS_IN_SET);
            if (skip && paths.isEmpty()) {
                return false;
            }
            ++processed;
            preProcess(id, paths);
            writer.printOpenTag("tr");
            writer.printThTaggedValue(id);
            writer.printThTaggedValue(database.getSpectraDb().getSpectrum(id).envelopes.length);
            Protein protein = bestFromAlign == null ? database.getProteinPredictedByAlign(id) : bestFromAlign.b;
            writer.printThTaggedValue(TagGenerator.getLastRatio());
            if (dontMatch) {
                writer.printTaggedValue("td", "?");
            } else {
                printProteinDBMatches(writer, paths);
            }
            writer.printThTaggedValue(String.format("%.2E", bestFromAlign == null ? database.getEValue(id) : bestFromAlign.a));
            printTags(writer, count, found, nTags, monoTags, protein, paths);
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

        protected void printTags(HtmlWriter writer, int[][] count, int[] found, int[] nTags, int[] monoTags, Protein protein, Collection<Path> paths) {
            int pathN = 0;
            boolean foundFirst = false;
            for (Path path : paths) {
                if (pathN == MAX_TAGS_IN_SET) {
                    break;
                }
                ++nTags[pathN];
                if (path.isMonoTag()) {
                    ++monoTags[pathN];
                }
                writer.printOpenTag("td");
                writer.printOpenTag("div", "align=center");
                writer.printTagPrefix("span");
                boolean notPrint = false;
                Double bestDiff = null;
                if (protein != null && protein.contains(path)) {
                    if (printShifts) {
                        bestDiff = protein.getLastBestShift();
                    }
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
                    writer.println(path + " " + (protein == null ? 0 : protein.getMaxMatch(path)) + " " + (bestDiff == null ? "" : String.format("%.2f", bestDiff)) + " " + String.format("%.2f", path.beginMass));
                    writer.printCloseTag("span");
                }
                writer.printCloseTag("div");
                writer.printCloseTag("td");
                ++pathN;
            }
            for (; pathN < MAX_TAGS_IN_SET; ++pathN) {
                writer.printEmptyTag("td");
            }
        }

        protected void printMatchedProteinsStats(HtmlWriter writer, int id) {
            Protein proteinFromAlign = bestFromAlign == null ? database.getProteinPredictedByAlign(id) : bestFromAlign.b;
            boolean predictedFound = !dontMatch && database.getProteinDb().getLastMatchedProteins().contains(proteinFromAlign.getFullname());
            writer.printTaggedValue("td", predictedFound ? "+" : "-");
            if (predictedFound) {
                ++predictedProteinFound;
            }
            if (!dontSearchForBest) {
                Protein bestMatch = dontMatch ? null : database.getBestMatch(id);
                boolean matchesAlign = proteinFromAlign.getFullname().equals(bestMatch == null ? null : bestMatch.getFullname());
                if (matchesAlign) {
                    ++matchedAlign;
                }
                writer.printOpenTag("td");
                writer.printTaggedValue("div", bestMatch == null ? null : bestMatch.getName(), "align=center" + (matchesAlign ? " style=\"color:red\"" : ""));
            } else {
                writer.printOpenTag("td");
            }
            writer.printTaggedValue("div", proteinFromAlign.getName());
            writer.printCloseTag("td");
        }
    }


    protected void printProteinDBMatches(HtmlWriter writer, Collection<Path> paths) {
        writer.printThTaggedValue(database.getProteinDb().getMatchedProteinsNumber(paths));
    }

    private String getOutputFilename() {
//        return (ConfigReader.getBooleanProperty("USE_VIRTUAL_SPECTRA") ? "virt" : "exp") + "_" + TagGenerator.MIN_TAG_LENGTH + "-" + TagGenerator.MAX_TAG_LENGTH + ".html";
        return (EDGE_OF_THREE_AA ? "3" : EDGE_OF_TWO_AA ? "2" : "1") + "_" + (USE_RED_BLUE_GRAPH ? "doub" : "for") + "_" + ((int) (MassUtil.ERROR_THRESHOLD * 2e6)) + "ppm_" + (TagGenerator.SCORE_BY_LENGTH ? "len" : "score") + "_" + (ConfigReader.getBooleanProperty("USE_VIRTUAL_SPECTRA") ? "virt" : "exp") + ".html";
    }
}