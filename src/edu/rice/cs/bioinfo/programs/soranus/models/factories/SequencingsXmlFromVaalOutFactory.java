package edu.rice.cs.bioinfo.programs.soranus.models.factories;

import edu.rice.cs.bioinfo.library.language.vaal._1_1.out.reading.ir.LineDifferenceRecord;
import edu.rice.cs.bioinfo.library.language.vaal._1_1.out.reading.parsers.VAALOutStringSplitParser;
import edu.rice.cs.bioinfo.library.language.vaal.out.reading.OutReadResult;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.text.ParseException;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 3/29/13
 * Time: 4:17 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SequencingsXmlFromVaalOutFactory<T>
{
    class NucleotideLocation implements Comparable<NucleotideLocation>
    {
        public final int Contig;

        public final int Position;

        public NucleotideLocation(int contig, int position) {
            Contig = contig;
            Position = position;
        }

        @Override
        public boolean equals(Object obj)
        {
            try
            {
                return equals((NucleotideLocation)obj);
            }
            catch(ClassCastException e)
            {
                return false;
            }
        }

        @Override
        public int hashCode()
        {
            return Position;
        }

        public boolean equals(NucleotideLocation other)
        {
            return other.Contig == this.Contig && other.Position == this.Position;
        }

        public int compareTo(NucleotideLocation o) {
            int contigCompare = this.Contig - o.Contig;

            if(contigCompare != 0)
                return contigCompare;

            return this.Position - o.Position;
        }

        public String toString()
        {
            return Position + "";
        }
    }


    public String makeXML(Map<T,String> sequencingEventToVaalOutputs) throws ParseException
    {
        Map<NucleotideLocation, Character> referenceGenome = new HashMap<NucleotideLocation, Character>();
        Map<T, Map<NucleotideLocation,Character>> sequencingEventToDifferences = new HashMap<T, Map<NucleotideLocation,Character>>();

        VAALOutStringSplitParser outFileParser = new VAALOutStringSplitParser();

        for(T sequencingEvent : sequencingEventToVaalOutputs.keySet())
        {

            String vaalOutput = sequencingEventToVaalOutputs.get(sequencingEvent);
            OutReadResult<List<LineDifferenceRecord>> outContents;
            try
            {
                outContents = outFileParser.parse(new ByteArrayInputStream(vaalOutput.getBytes()));
            }
            catch (IOException e)
            {
                throw new RuntimeException(e);
            }

            HashMap<NucleotideLocation, Character> differencesForSequencingEvent = new HashMap<NucleotideLocation, Character>();
            sequencingEventToDifferences.put(sequencingEvent, differencesForSequencingEvent);
            for(LineDifferenceRecord difference : outContents.getDifferences())
            {
                if(difference.Sample.length() > 1)
                    throw new RuntimeException("non SNP"); // TODO do something better

                NucleotideLocation differencePosition = new NucleotideLocation(difference.Contig, difference.Position);
                if(!referenceGenome.containsKey(differencePosition))
                {
                    referenceGenome.put(differencePosition, difference.Ref.charAt(0));
                }
                else if(!referenceGenome.get(differencePosition).equals(difference.Ref.charAt(0)))
                {
                    throw new RuntimeException(); // TODO: promote to better exception
                }

                differencesForSequencingEvent.put(differencePosition, difference.Sample.charAt(0));
            }
        }

        List<NucleotideLocation> genomeWideDifferenceLocations = new LinkedList<NucleotideLocation>(referenceGenome.keySet());
        Collections.sort(genomeWideDifferenceLocations);

        StringBuffer sequencingsText = new StringBuffer("<Sequencings>\n");
        for(T sequencingEvent : sequencingEventToVaalOutputs.keySet())
        {
            Map<NucleotideLocation,Character> positionToDifference = sequencingEventToDifferences.get(sequencingEvent);

            StringBuffer sequenceAttributeValue = new StringBuffer();
            for(NucleotideLocation differenceLocation : genomeWideDifferenceLocations)
            {
                sequenceAttributeValue.append(positionToDifference.containsKey(differenceLocation) ?
                        positionToDifference.get(differenceLocation) :
                        referenceGenome.get(differenceLocation));
            }

            String sourceId = getSourceId(sequencingEvent);
            sequencingsText.append("\t<Sequencing sourceId=\"" + sourceId +"\" ");
            sequencingsText.append("sequence=\"" + sequenceAttributeValue +"\"/>\n");
        }
        sequencingsText.append("</Sequencings>");

        return sequencingsText.toString();
    }

    protected abstract String getSourceId(T sequencingEvent);
}


