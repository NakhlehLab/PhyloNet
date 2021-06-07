package edu.rice.cs.bioinfo.programs.soranus.viewModels;

import java.util.LinkedList;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/31/13
 * Time: 6:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class SequencesAndSourceVM
{
    private class Record
    {
        public final Object Source;

        public final String Sequence;

        private Record(Object source, String sequence) {
            Source = source;
            Sequence = sequence;
        }
    }

    private LinkedList<Record> _sequences = new LinkedList<Record>();

    public void addSequence(String sequence, Object source)
    {
        _sequences.add(new Record(source, sequence));
    }
}
