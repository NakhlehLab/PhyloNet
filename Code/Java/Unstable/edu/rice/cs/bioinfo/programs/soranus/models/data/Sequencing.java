package edu.rice.cs.bioinfo.programs.soranus.models.data;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/3/13
 * Time: 5:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class Sequencing
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
