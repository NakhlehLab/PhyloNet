package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.datatype.DataType;

import java.util.Collection;
import java.util.List;

/**
 * Created by wendingqiao on 5/3/16.
 */
public class Sequence {

    private int _totalCount; // number of states or the number of lineages for this species
    private String _taxon;
    private String _sequenceData;

    public Sequence(String taxon, String sequenceData) {
        _totalCount = -1;
        _taxon = taxon;
        _sequenceData = sequenceData;
    }

    public List<Integer> getSequence(DataType dataType) {
        List<Integer> sequence = dataType.stringToState(_sequenceData);
        _totalCount = dataType.getStateCount();
        return sequence;
    }

    /**
     * @return the taxon of this sequence as a string.
     */
    public final String getTaxon() {
        return _taxon;
    }

    /**
     * @return the data of this sequence as a string.
     */
    public final String getData() {
        return _sequenceData;
    }

    /**
     * @return return the state count
     */
    public final int getStateCount() {
        return _totalCount;
    }

    /**
     * @return the sequence in the collection with the given taxon, or null if its not in the collection.
     */
    public static Sequence getSequenceByTaxon(String taxon, Collection<Sequence> sequences) {
        for (Sequence seq : sequences) {
            if (seq.getTaxon().equals(taxon)) return seq;
        }
        return null;
    }

}