package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.Identifier;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.PySONNode;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentSet;

import java.util.AbstractMap;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 12/12/11
 * Time: 2:35 PM
 * To change this template use File | Settings | File Templates.
 */
public class DataBlockBody implements Block
{
    public final LinkedList<AbstractMap.SimpleImmutableEntry<Identifier, Identifier>> SequencePairs;

    public DataBlockBody(LinkedList<AbstractMap.SimpleImmutableEntry<Identifier, Identifier>> sequencePairs)
    {
        SequencePairs = sequencePairs;
    }

    public <R, T, E extends Exception> R execute(BlockAlgo<R, T, E> algo, T input) throws E
    {
        return algo.forDataBlock(this, input);
    }
}
