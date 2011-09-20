package ru.spbau.ablab.tagfinder.database.proteindb;

import ru.spbau.ablab.tagfinder.Protein;
import ru.spbau.ablab.tagfinder.util.io.FastScanner;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

public class ExperimentalProteinDb extends ProteinDb {
    public ExperimentalProteinDb() throws FileNotFoundException {
        FastScanner scanner = new FastScanner(new File(FASTA_PATH));
        fullnameToProtein = new HashMap<String, Protein>();
        String proteinName = null;
        StringBuilder proteinString = null;
        scanIdToProtein = new HashMap<Integer, Protein>();
        for (String line; (line = scanner.nextLine()) != null; ) {
            if (line.charAt(0) == '>') {
                if (proteinName != null) {
                    Protein protein = new Protein(proteinString.toString().trim().replace('I', 'L'), proteinName);
                    scanIdToProtein.put(scanIdToProtein.size(), protein);
                    fullnameToProtein.put(proteinName, protein);
                }
                proteinName = line.substring(1);
                proteinString = new StringBuilder();
            } else {
                assert proteinString != null;
                proteinString.append(line);
            }
        }
        if (proteinName != null) {
            fullnameToProtein.put(proteinName, new Protein(proteinString.toString().trim().replace('I', 'L'), proteinName));
        }
        nameToProtein = new HashMap<String, Protein>();
        proteins = new HashSet<Protein>();
        for (Map.Entry<String, Protein> entry : fullnameToProtein.entrySet()) {
            Protein protein = entry.getValue();
            nameToProtein.put(protein.getName(), protein);
            proteins.add(protein);
        }
    }
}
