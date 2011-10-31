package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

import java.util.Iterator;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 5:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class IdentList extends PySONNodeLineAndCol
{
    private final LinkedList<Identifier> _elements = new LinkedList<Identifier>();

    public final Iterable<Identifier> Elements = _elements;

    public final int ElementsCount;

    public IdentList(int line, int col, Iterable<Identifier> elements)
    {
        super(line, col);
        for(Identifier ident : elements)
        {
            _elements.add(ident);
        }

        ElementsCount = _elements.size();
    }
}
