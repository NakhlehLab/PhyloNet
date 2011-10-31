package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/17/11
 * Time: 5:00 PM
 * To change this template use File | Settings | File Templates.
 */
public class TaxaMap implements PySONNode{

    public Iterable<TaxaMapEntry> Entries;

    public final int Line;

    public final int Col;

    public TaxaMap(int line, int col, Iterable<TaxaMapEntry> entries)
    {
        Line = line;
        Col = col;
        Entries = entries;
    }
}
