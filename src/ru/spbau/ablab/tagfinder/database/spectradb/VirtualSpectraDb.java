package ru.spbau.ablab.tagfinder.database.spectradb;

import ru.spbau.ablab.tagfinder.database.Database;
import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.MassUtil;
import ru.spbau.ablab.tagfinder.util.StringUtil;
import ru.spbau.ablab.tagfinder.util.io.FastScanner;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import static ru.spbau.ablab.tagfinder.database.Database.ALIGN_RESULT_FILE;
import static ru.spbau.ablab.tagfinder.database.Database.ALIGN_SPECTRA_FILE;

public class VirtualSpectraDb extends SpectraDb {
    private Map<Integer, Double> idToEValue;

    public Map<Integer, Double> getIdToEValue() {
        return idToEValue;
    }

    public VirtualSpectraDb() throws FileNotFoundException {
        assert scanToSpectrum != null;
        idToEValue = new HashMap<Integer, Double>();
        HashMap<Integer, Integer> idToScan = getIdToScanMap();
        FastScanner scanner = new FastScanner(new File(ALIGN_RESULT_FILE));
        scanToSpectrum = new TreeMap<Integer, Spectrum>();
        Database database = Database.getInstance();
        while (scanner.hasNextLine()) {
            scanner.skipLine("BEGIN PRSM");
            int spectrumId = scanner.getNextIntProperty("SPECTRUM_ID");
            boolean toAdd = idToScan.containsKey(spectrumId);
            int id = toAdd ? idToScan.get(spectrumId) : -1;
            double eValue = scanner.getNextDoubleProperty("E_VALUE");
            idToEValue.put(id, eValue);
            double parentMass = toAdd ? idToPrecursorMass.get(id) : Double.NaN;
            scanner.skipLine("BEGIN MATCH_PAIR");
            ArrayList<Envelope> envelopes = new ArrayList<Envelope>();
            for (String s; !(s = scanner.nextLine().trim()).equals("END MATCH_PAIR"); ) {
                String[] strings = StringUtil.getTokenArray(s);
                double mass = Double.parseDouble(strings[3]);
                mass = strings[4].equals("B") ? mass : (MassUtil.convertIonsType(mass, parentMass));
                if (toAdd) {
                    Envelope closestExp = database.getClosestEnvelopeFromExperimental(id, mass);
                    envelopes.add(new Envelope(mass, closestExp.score, closestExp.intensity));
                }
            }
            if (toAdd) {
                scanToSpectrum.put(id, new Spectrum(id, envelopes.toArray(new Envelope[envelopes.size()]), parentMass));
            }
            scanner.skipLine("END PRSM");
        }
        spectraIds = scanToSpectrum.keySet();
        scanner.close();
    }

    private HashMap<Integer, Double> idToPrecursorMass;

    private HashMap<Integer, Integer> getIdToScanMap() throws FileNotFoundException {
        FastScanner scanner = new FastScanner(new File(ALIGN_SPECTRA_FILE));
        HashMap<Integer, Integer> ans = new HashMap<Integer, Integer>();
        idToPrecursorMass = new HashMap<Integer, Double>();
        while (scanner.skipLine("BEGIN IONS")) {
            int id = scanner.getNextIntProperty("ID");
            int scan = scanner.getNextIntProperty("SCANS");
            ans.put(id, scan);
            idToPrecursorMass.put(scan, scanner.getNextDoubleProperty("PRECURSOR_MASS"));
        }
        scanner.close();
        return ans;
    }
}
