import junit.framework.TestCase;
import ru.spbau.ablab.tagfinder.StatisticsGenerator;
import ru.spbau.ablab.tagfinder.UnmatchedStatistics;
import ru.spbau.ablab.tagfinder.util.ConfigReader;

public class TagFinderTest extends TestCase {
    public void testStatGenByRun() {
        StatisticsGenerator.main(null);
    }

    public void testUnmatchedStatGenByRun() {
        UnmatchedStatistics.main(null);
    }


    public void testMonoTags() {
        ConfigReader.setProperty("DOUBLE_MASSES","false");
        StatisticsGenerator.main(null);
    }
}