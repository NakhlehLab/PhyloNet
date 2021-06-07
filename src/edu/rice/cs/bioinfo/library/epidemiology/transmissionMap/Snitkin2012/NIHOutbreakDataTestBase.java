package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import org.joda.time.LocalDate;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/22/13
 * Time: 7:40 PM
 * To change this template use File | Settings | File Templates.
 */
public class NIHOutbreakDataTestBase
{
    protected class Sequencing
    {
        public final String Sequence;

        public Sequencing(String sequence)
        {
            Sequence = sequence;
        }

        public int getGeneticDistance(Sequencing other)
        {
            if(Sequence.length() != other.Sequence.length())
            {
                throw new IllegalArgumentException("Sequences must have same length");
            }

            int variantsAccum = 0;
            for(int i = 0; i<Sequence.length(); i++)
            {
                if(Sequence.charAt(i) != other.Sequence.charAt(i))
                    variantsAccum++;
            }

            return variantsAccum;
        }
    }

    // rows are days, columns are patients, integer entries are locations, .5 suffix is "fist test positive day" for the patient
    protected final String[] patientTraces = new String[] {
            "0	0	0	0	0	0	0	0	0	5	4	4	0	0	0	1	0	0",
            "0	0	0	0	0	0	0	0	0	0	4	4	8	0	0	1	0	0",
            "0	0	0	8	0	0	0	0	0	0	4	4	8	0	0	1	0	1",
            "0	0	0	8	0	0	0	0	0	0	4	4	8	0	0	1	0	1.5",
            "0	0	0	8	0	0	0	0	0	0	4	4	0	0	0	1	0	8",
            "0	0	0	8	0	0	1	0	0	0	4	4	0	0	0	1	0	8",
            "0	0	0	8	0	0	1	0	0	0	4	4	0	0	0	1	0	8",
            "0	0	0	8	0	0	5	0	0	0	4	4	0	0	0	1	0	8",
            "0	0	0	8	0	0	5	0	0	0	4	4	0	0	0	1	0	8",
            "0	0	0	8	0	0	5	0	0	0	4	4	0	0	0	1	0	8",
            "0	0	0	0	0	0	5	0	0	0	4	4	0	0	0	1	0	8",
            "0	0	0	0	0	0	5	0	0	0	4	4	0	0	6	1	0	8",
            "0	0	0	0	0	0	5	0	0	0	4	4	0	0	6	1	0	8",
            "0	0	0	0	0	0	0	0	0	0	4	4	0	0	6	1	0	8",
            "0	0	0	0	7	0	5	0	0	8	4	4	0	0	6	1	0	8",
            "0	0	0	0	7	0	5	0	0	8	4	4	0	0	6	1	0	8",
            "0	0	0	0	7	0	0	0	0	8	4	4	0	0	6	1	0	8",
            "0	0	0	0	0	0	0	0	0	8	0	4	0	0	6	1	0	1",
            "0	0	0	0	0	0	0	0	0	8	0	4	0	0	6	1	0	1",
            "0	0	0	0	0	0	0	0	0	0	0	4	0	0	6	1	0	8",
            "0	0	0	0	0	0	0	0	0	0	0	4	0	0	6	1	0	8",
            "0	0	0	0	0	0	0	0	0	0	0	4	0	0	6	1	0	8",
            "0	0	0	0	0	0	0	0	0	0	0	4	0	0	6	1	0	8",
            "0	0	0	0	0	0	0	4	0	0	0	4	0	0	6	1	0	8",
            "0	0	0	0	0	0	1	4	0	0	4	0	0	0	6	1	0	8",
            "0	0	0	0	7	0	1	4	0	0	4	0	0	0	6	1	0	8",
            "0	0	0	0	7	0	1	4	0	0	4	0	0	0	6	1	0	8",
            "0	0	0	0	7	0	1	4	0	0	4	0	0	0	6	1	0	8",
            "0	0	0	0	7	0	1	4	0	0	4	0	0	0	6	1	0	8",
            "0	4	0	0	7	0	1	4	0	0	4	0	0	0	6	1	5	8",
            "0	4	0	0	7	0	1	4	0	0	0	0	0	5	6	1	5	8",
            "0	4	0	0	7	0	1	4	0	0	0	0	0	5	6	1	5	8",
            "0	4	0	0	7	0	1	0	0	0	0	0	0	5	6	1	5	8",
            "0	4	0	0	0	0	1	0	0	0	0	0	0	5	6	1	5	8",
            "0	4	0	0	0	0	1	0	0	0	0	0	0	5	6	1	5	0",
            "0	4	0	0	0	0	1	0	0	0	0	0	0	5	6	1	5	0",
            "0	4	0	0	0	0	1	0	0	0	0	0	0	5	6	1	5	0",
            "0	4	0	0	0	0	1	0	0	0	0	0	0	5	6	1	5	0",
            "0	4	0	0	0	0	1	4	0	0	0	0	0	5	6	1	1	0",
            "0	4	0	0	0	0	1	4	0	0	0	0	0	5	6	1	1	0",
            "0	4	0	0	0	0	1	1	0	0	0	0	0	5	6	1	1	0",
            "0	4	0	0	0	0	1	1	0	0	0	0	0	5	6	1	1	0",
            "0	4	0	0	0	0	1	1	0	0	0	0	0	5	6	1	1	0",
            "0	4	0	0	0	0	1	1	0	0	0	0	0	1	6	1	1	0",
            "0	4	0	0	0	0	3	1	0	0	0	0	0	1	6	4	1	0",
            "0	4	0	0	0	0	1	1	0	0	0	4	0	1	6	4	1	0",
            "0	4	0	0	0	0	1	4	0	0	0	4	0	1	6	4	1	0",
            "0	4	0	0	0	0	1	4	0	0	0	4	0	1	6	4	1	0",
            "0	4	0	0	0	0	1	4	0	0	0	4	0	1	6	4	1	0",
            "0	4	0	0	0	0	1	4	0	0	0	4	0	1	6	4	1	0",
            "0	4	0	0	0	0	1	4	0	0	0	4	0	1	6	4	1	0",
            "0	4	0	0	7	0	1	4	0	8	0	4	0	1	6	4	1	0",
            "0	4	0	0	7	0	1	4	0	8	0	4	0	1	6	4	1	0",
            "0	4	0	0	7	0	1	4	0	8	0	4	0	1	6	4	1	0",
            "0	4	0	0	7	0	1	4	0	8	0	4	7	1	6	4	1	0",
            "0	4	0	0	7	0	1	4	0	8	0	4	7	1	6	4	1	0",
            "0	4	0	0	7	0	1	4	0	8	0	4	7	1	6	4	1.5	0",
            "0	4	0	0	7	0	1	4	0	8	0	0	4	1	6	4	1	0",
            "0	4	0	0	7	0	1	4	5	8	0	0	4	1	6	4	1	0",
            "0	4	0	0	0	0	1	4	5	8	0	0	4	1	6	4	1	0",
            "0	0	0	0	0	0	1	4	5	8	0	0	4	1	6	1	1	0",
            "0	0	0	0	0	0	1	4	5	8	0	4	4	1	1	1	1	0",
            "0	0	0	0	0	0	1	4	5	0	0	4	4	1	1	1	1	0",
            "0	0	0	0	0	0	1	4	5	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	0	0	1	4	5	0	0	4	4	1	6	1.5	1	0",
            "0	0	0	0	0	0	1	4	5	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	0	0	1	4	5	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	0	0	1	0	5	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	0	0	1	0	5	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	0	0	1	0	1	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	0	0	1	0	1	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	0	0	1	4	1	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	7	0	1	4	1	0	0	4	4	1	6.5	1	1	0",
            "0	0	0	0	8	0	1	4	1	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	8	0	1	4	1	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	8	0	1	4	1	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	8	0	1	4	1	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	8	0	1	4	1	0	0	4	4	1	6	1	1	0",
            "0	0	0	0	8	0	1	4	1	0	0	4	4	1.5	6	1	1	0",
            "0	0	0	8	8	0	1	4	1	0	0	4	4	1	1	1	1	0",
            "0	0	0	8	8	0	1	4	1	0	0	4	4	1	1	1	1	0",
            "0	0	0	8	8	0	1	4	1	0	0	4	0	1	1	1	1	0",
            "0	0	0	8	7	0	1	0	1	0	0	4	0	1	1	1	1	0",
            "0	0	0	8	7	0	1	0	1	0	0	4	4	2	2	2	2	0",
            "0	0	0	8	7	0	1	0	1	0	0	4	4	2	2	2	2	0",
            "0	0	0	8	7	0	1	0	1	0	0	4	4	2	2	2	2	0",
            "0	0	0	8	7	0	1	0	1	0	0	4	4	2	2	2	2	0",
            "0	0	0	8	7	0	1	0	1	0	0	4	4	2	2	2	2	0",
            "0	0	0	8	7	0	1	0	1	0	0	4	1	2	2	2	0	0",
            "0	4	0	8	7	0	1	0	1	0	0	4	1	2	2	2	0	0",
            "0	4	0	8	7	0	5	0	1	0	0	4	1	2	2	2	0	0",
            "0	4	0	8	7	0	5	4	1	0	0	4	1	2	2	2	0	0",
            "0	4	0	8	7	0	5	4	1	0	0	4	1	2	2	2	0	0",
            "0	0	0	8	7	0	5	4	1	0	0	4	1	2	2	2	0	0",
            "0	0	0	8	7	0	5	4	1	0	0	4	1	2	2	2	0	0",
            "0	0	0	8	7	0	5	4	1	0	0	4	1.5	2	2	2	0	0",
            "0	0	0	8	7	0	5	4	1	0	0	1	1	2	2	2	0	0",
            "0	0	0	8	7	0	5	4	1	0	0	1	1	2	2	2	0	0",
            "0	0	0	8	8	0	5	4	1	0	4	1	0	2	2	2	0	0",
            "0	0	0	8	8	0	5	4	1	0	4	1.5	1	2	2	2	0	0",
            "0	0	0	8	8	0	5	4	1	8	4	1	0	2	2	2	0	0",
            "0	0	0	8	8	0	1	4	1	8	4	2	0	2	2	2	0	0",
            "0	0	8	8	8	0	1	4	1	8	4.5	2	0	0	2	2	0	0",
            "0	0	8	8	8	0	1	4	1	1	4	2	0	2	2	2	0	0",
            "0	0	8	8	4	0	1	1	1	1	4	2	0	2	2	2	0	0",
            "0	0	8	8	4	0	1	1	1	1	4	2	0	2	2	2	0	0",
            "0	0	8	8	4	0	2	1	1	1.5	2	2	0	2	2	2	0	0",
            "0	0	8	8	4	0	2	4	1	1	2	2	0	2	0	2	0	0",
            "0	0	8	8	4	0	2	4	1	1	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	1	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	1	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	1	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	1	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	1	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	1	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	1	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	1.5	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	1	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	2	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	2	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	2	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	2	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	2	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	4	2	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	1	2	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	1	2	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	1	2	2	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	0	2.5	1.5	2	0	2	0	0	2	0	2	0	0",
            "0	0	8	0	4	0	2	1	2	0	2	2	0	2	0	2	0	0",
            "0	4	8	0	4	0	2	1	2	0	2	2	0	2	0	2	0	0",
            "0	4	8	0	4	5	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	4	8	0	4	5	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	4	8	0	4	5	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	4	8	0	4	5	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	4	8	0	4	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	0	8	0	4	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	0	8	0	1	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	0	8	0	1	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	0	8	0	1	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	0	8	0	1	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	4	8	0	1	1.5	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	4	8	0	1	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	4	8	0	1	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	4	8	0	1	1	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	4	8	0	1	2	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	6	8	0	1	2	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	6	8	0	1	2	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	6	8	8	4.5	2	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	6	8	8	2	2	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	6	8	8	2	2	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	6	8	8	2	2	2	2	2	0	2	2	0	2	0	2	0	0",
            "0	6	8	8	2	2	2	2	2	0	2	4	0	2	0	2	0	0",
            "5	6	8	8	2	2	2	2	2	0	2	2	0	2	0	2	0	0",
            "5	0	8	1	2	2	2	2	2	0	2	2	0	2	0	2	0	0",
            "5	0	8	1.5	2	2	2	2	2	0	2	0	0	2	0	2	0	0",
            "5	0	8.5	8	2	2	2	2	2	0	2	0	0	2	0	2	0	0",
            "5	0	8	2	2	2	2	2	2	0	2	0	0	2	0	2	0	0",
            "5	0	8	2	2	2	2	2	2	0	2	0	0	2	0	2	0	0",
            "5	0	8	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	0	8	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	0	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	0	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	0	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	0	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	5.5	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	5	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	2	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	0	2	2	0	0	2	0	0	2	0	2	0	0",
            "5.5	2	2	2	2	0	2	2	0	0	2	0	0	2	0	2	0	0",
            "5	2	2	2	2	0	2	2	0	0	2	0	0	2	0	2	0	0",
            "0	2	2	2	2	0	2	2	0	0	2	0	0	2	0	2	0	0",
            "0	2	2	2	2	0	2	2	0	0	2	0	0	2	0	2	0	0",
            "0	2	2	2	2	0	2	2	0	0	2	0	0	2	0	2	0	0",
            "0	2	2	2	2	0	2	2	0	0	2	0	0	2	0	2	0	0",
            "0	2	2	0	2	0	2	2	0	0	0	0	0	2	0	2	0	0",
            "0	2	2	0	2	0	2	2	0	0	0	0	0	2	0	2	0	0",
            "0	2	2	0	2	0	2	0	0	0	0	0	0	2	0	2	0	0",
            "0	2	2	0	2	0	2	0	0	0	0	0	0	2	0	2	0	0",
            "0	2	2	0	2	0	2	0	0	0	0	0	0	2	0	2	0	0",
            "0	2	2	0	2	0	2	0	0	0	0	0	0	2	0	0	0	0",
            "0	2	2	0	2	0	2	0	0	0	0	0	0	2	0	0	0	0",
            "0	2	2	0	2	0	2	0	0	0	0	0	0	2	0	0	0	0",
            "0	2	2	0	2	0	2	0	0	0	0	0	0	2	0	0	0	0",
            "0	2	2	0	2	0	2	0	0	0	0	0	0	2	0	0	0	0",
            "0	2	2	0	2	0	2	0	0	0	0	0	0	2	0	0	0	0",
            "0	2	2	0	2	0	2	0	0	0	0	0	0	2	0	0	0	0",
            "2	2	2	0	2	0	2	0	0	0	0	0	0	2	0	0	0	0",
            "2	2	2	0	2	0	2	0	0	0	0	0	0	2	0	0	0	0",
            "2	2	2	0	2	0	2	0	0	0	0	0	0	2	0	0	0	0",
    };


    // rows are patients, columns are loci.
    private final String[] sequencingData = new String[]
            {
    /* 1 urine */    "00000000000000000000000000000000000000000",
    /* 1 bal */      "00000000000000000000000000000000000001011",
    /* 1 urine */    "00000000000000000000000000000000000000000",
    /* 1 urine */    "00000000000000000000000000000000000000000",
    /* 1 urine */    "00000000000000000000000000000000000000000",
    /* 1 groin */    "00000000000000000000000000000000000001111",
    /* 1 throat */   "00000000000000000000000000000000011100000",
    /* 2  */         "00000000000000000000000001000000011110000",
    /* 3 groin */    "00000000000000000000000000000000011110000",
    /* 4 throat */   "01100100000001001101000010000000100001011",
    /* 5 */          "00000000000000000000000000000000011110000",
    /* 6 */          "00010110000001001101000010001000100001011",
    /* 7 */          "00010110000001001101000010001000100001011",
    /* 8 */          "10000000010000110000100100000000000001011",
    /* 9 */          "00000100000101001101000010110100100001011",
    /* 10 */         "00010110000001001101000010001000100001011",
    /* 11 */         "00011110000001001101000010001011100001011",
    /* 12 */         "00010111000001001111000010001000100001011",
    /* 13 */         "00010111000001001111000010001000100001011",
    /* 14 */         "00010110000001001101000010001000100001011",
    /* 15 */         "00010110000001001101000010001000100001011",
    /* 16 */         "00010110100011001101000010001000100001011",
    /* 17 */         "00010110000001001101000010001000100001011",
    /* 18 */         "00010111001001001111000010011000100001011",
            };

    protected final HashMap<Sequencing,Integer> sequenceToPatient;

    protected final Map<Integer,LocalDate> patientIdToFirstPositive;

    protected final Map<Integer, Map<LocalDate,Object>> patientTraceData;

    protected NIHOutbreakDataTestBase()
    {

        sequenceToPatient = new HashMap<Sequencing, Integer>();

        int sequencingDataIndex = 0;
        for(int patientId : Arrays.asList(1,1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))
        {
            sequenceToPatient.put(new Sequencing(sequencingData[sequencingDataIndex++].replace(" ", "")), patientId);
        }

        int numPatients = patientTraces[0].split("\\s+").length;

        patientIdToFirstPositive = new HashMap<Integer, LocalDate>();

        patientTraceData = new HashMap<Integer, Map<LocalDate, Object>>();
        for(int patientId = 1; patientId<=numPatients; patientId++)
        {
            patientTraceData.put(patientId, new HashMap<LocalDate, Object>());
        }


        LocalDate day = new LocalDate(2011, 6, 12);
        for(String patientLocationsOnDay : patientTraces)
        {
            String[] patientLocationOnDaySplit = patientLocationsOnDay.split("\\s+");
            for(int patientId = 1; patientId<=numPatients; patientId++)
            {
                String patientLocationOnDayWithTestPositiveFlag = patientLocationOnDaySplit[numPatients-patientId];
                boolean firstDayToTestPositive = patientLocationOnDayWithTestPositiveFlag.endsWith(".5");

                if(firstDayToTestPositive)
                {
                    patientIdToFirstPositive.put(patientId, day);
                }

                String patientLocationOnDay = patientLocationOnDayWithTestPositiveFlag.substring(0, 1);

                if(!patientLocationOnDay.equals("0"))
                {
                    patientTraceData.get(patientId).put(day, patientLocationOnDay);
                }
            }

            day = day.plusDays(1);
        }



    }

    protected Map<Integer,Map<Integer,Integer>> makeEpiDistanceMap(Integer eMax)
    {
        Map<Integer,Map<Integer,Integer>> epiDistanceMatrixAccum = new HashMap<Integer,Map<Integer,Integer>>();

        int numPatients = patientTraces[0].split("\\s+").length;

        for(int patientId = 1; patientId<=numPatients; patientId++)
        {
            epiDistanceMatrixAccum.put(patientId, new HashMap<Integer, Integer>());
        }

        for(Map<Integer,Integer> value : epiDistanceMatrixAccum.values())
        {
            for(int patientId = 1; patientId<=numPatients; patientId++)
            {
                value.put(patientId, Integer.MAX_VALUE);
            }
        }

        Map<Integer,Integer> patientIdToFirstPositiveDay = new HashMap<Integer, Integer>();

        for(int day = 0; day<patientTraces.length; day++)
        {
            String[] dayTrace = patientTraces[day].split("\\s+");

            for(int patientIndex = 0; patientIndex<numPatients; patientIndex++)
            {
                if(dayTrace[patientIndex].endsWith(".5"))
                {
                    patientIdToFirstPositiveDay.put(numPatients-patientIndex, day);
                }
            }
        }

        boolean[][] haveOverlapped = new boolean[numPatients][numPatients];
        for(int a = 0; a<numPatients; a++)
        {
            for(int b = 0; b<numPatients; b++)
            {
                haveOverlapped[a][b] = false;
            }
        }

        for(int day = 0; day<patientTraces.length; day++)
        {
            String[] patientLocationsOnDay = patientTraces[day].split("\\s+");

            for(int donorIndex = 0; donorIndex<numPatients; donorIndex++)
            {
                int donor = numPatients-donorIndex;
                String donorLocation = patientLocationsOnDay[donorIndex].replace(".5", "");

                if(donorLocation.equals("0"))
                    continue;

                for(int recipIndex = 0; recipIndex < numPatients; recipIndex++)
                {
                    int recip = numPatients-recipIndex;
                    String recipLocation = patientLocationsOnDay[recipIndex].replace(".5", "");

                    if(recipLocation.equals("0"))
                        continue;


                    if(donor == 1 && recip == 3)
                    {
                        int j = 0;
                    }

                    if(donorLocation.equals(recipLocation))
                    {


                        int donorFirstPositiveDay = patientIdToFirstPositiveDay.get(donor);
                        int recipFirstPositiveDay = patientIdToFirstPositiveDay.get(recip);

                        if(recipFirstPositiveDay < day)
                        {
                            if(!haveOverlapped[donorIndex][recipIndex])
                            {
                                epiDistanceMatrixAccum.get(donor).put(recip, eMax);
                            }
                            else
                            {
                                continue;
                            }
                        }


                        if(donorFirstPositiveDay < day && recipFirstPositiveDay < day)
                        {
                            continue;
                        }


                        int donorSilentCol;
                        if(donorFirstPositiveDay > day)
                        {
                            donorSilentCol =  donorFirstPositiveDay - day;
                        }
                        else
                        {
                            donorSilentCol = 0;
                        }

                        int recipSilentCol = recipFirstPositiveDay - day;

                        int totalSilentCol = donorSilentCol + recipSilentCol;

                        if(totalSilentCol < 0)
                        {
                            int j = 0;
                        }

                        if(totalSilentCol < epiDistanceMatrixAccum.get(donor).get(recip))
                        {

                            epiDistanceMatrixAccum.get(donor).put(recip, totalSilentCol);
                        }

                        haveOverlapped[donorIndex][recipIndex] = true;
                    }
                }
            }
        }

        for(int a = 0; a<numPatients; a++)
        {
            for(int b = 0; b<numPatients; b++)
            {
                if(!haveOverlapped[a][b])
                {
                    epiDistanceMatrixAccum.get(numPatients-a).put(numPatients-b, eMax);
                }
            }
        }

        return epiDistanceMatrixAccum;
    }


}
