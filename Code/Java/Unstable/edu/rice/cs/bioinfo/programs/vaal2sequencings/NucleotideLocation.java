package edu.rice.cs.bioinfo.programs.vaal2sequencings;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 3/21/13
 * Time: 3:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class NucleotideLocation implements Comparable<NucleotideLocation>
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
