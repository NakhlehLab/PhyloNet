package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/17/11
 * Time: 4:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class TaxaMapEntry implements PySONNode
{
    public final Identifier Key;

    public final Iterable<Identifier> Values;

    public TaxaMapEntry(Identifier key, Iterable<Identifier> values)
    {
        Key = key;
        Values = values;
    }


}
