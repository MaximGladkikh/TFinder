package ru.spbau.ablab.tagfinder.database.proteindb;

import ru.spbau.ablab.tagfinder.Protein;
import ru.spbau.ablab.tagfinder.util.ConfigReader;
import ru.spbau.ablab.tagfinder.util.io.FastScanner;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.HashSet;

public class VirtualProteinDb extends ProteinDb {
    private static final String ALIGN_RESULT_TABLE_FILE = ConfigReader.getProperty("ALIGN_RESULT_TABLE_FILE");

    public VirtualProteinDb() throws FileNotFoundException {
        FastScanner scanner = new FastScanner(new File(ALIGN_RESULT_TABLE_FILE));
        scanIdToProtein = new HashMap<Integer, Protein>();
        proteins = new HashSet<Protein>();
        for (String s; (s = scanner.nextLine()) != null; ) {
            String[] ss = s.split("\t+");
            int id = Integer.parseInt(ss[5].split(" +")[2]);
            String proteinString = ss[ss.length - 5].replace('I', 'L');
            String proteinName = ss[9];
            Protein protein = new Protein(proteinString, proteinName);
            proteins.add(protein);
            scanIdToProtein.put(id, protein);
        }
    }
}
