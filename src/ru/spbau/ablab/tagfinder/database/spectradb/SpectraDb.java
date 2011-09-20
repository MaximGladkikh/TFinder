package ru.spbau.ablab.tagfinder.database.spectradb;

import ru.spbau.ablab.tagfinder.spectrum.Spectrum;

import java.util.Map;
import java.util.Set;

public abstract class SpectraDb {
    protected Map<Integer, Spectrum> scanToSpectrum;
    protected Set<Integer> spectraIds;

    public Set<Integer> getSpectraIds() {
        return spectraIds;
    }

    public Spectrum getSpectrum(int id) {
        Spectrum spectrum = scanToSpectrum.get(id);
        if (spectrum == null) {
            throw new RuntimeException("No spectrum id=" + id);
        }
        return spectrum;
    }
}
