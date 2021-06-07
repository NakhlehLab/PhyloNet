package edu.rice.cs.bioinfo.library.language.vaal._1_1.out.reading.ir;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 3/21/13
 * Time: 2:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class LineDifferenceRecord {

    public final int Contig;

    public final int Position;

    public final String Sample;

    public final String Ref;

    public LineDifferenceRecord(int contig, int position, String sample, String ref)
    {
        Contig = contig;
        Position = position;
        Sample = sample;
        Ref = ref;
    }
}
