package ru.spbau.ablab.tagfinder.database.spectradb;

import ru.spbau.ablab.tagfinder.spectrum.Spectrum;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.TreeMap;

import static ru.spbau.ablab.tagfinder.database.Database.*;

public class ExperimentalSpectraDb extends SpectraDb {
    public ExperimentalSpectraDb() throws FileNotFoundException {
        scanToSpectrum = new TreeMap<>();
        File envDir = new File(ENVELOPES_DIR);
        for (File file : envDir.listFiles()) {
            String name = file.getName();
            if (!name.toLowerCase().startsWith(SPECTRUM_FILE_PREFIX) || !name.toLowerCase().endsWith(SPECTRUM_FILE_SUFFIX)) {
                continue;
            }
            String stringId = name.substring(SPECTRUM_FILE_PREFIX.length(), name.length() - SPECTRUM_FILE_SUFFIX.length());
            int id = Integer.parseInt(stringId);
            Spectrum spectrum = new Spectrum(id, file);
            scanToSpectrum.put(id, spectrum);
        }
        spectraIds = scanToSpectrum.keySet();
    }
}
