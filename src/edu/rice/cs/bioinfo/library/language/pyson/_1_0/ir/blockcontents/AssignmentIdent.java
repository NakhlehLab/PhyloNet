package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/14/12
 * Time: 2:01 PM
 * To change this template use File | Settings | File Templates.
 */
public class AssignmentIdent
{
    private final int _line;

    public int getLine()
    {
        return _line;
    }

    private final int _column;

    public int getColumn()
    {
        return _column;
    }

    public final String Identifier;

    public AssignmentIdent(int line, int column, String ident)
    {
        _line = line;
        _column = column;
        Identifier = ident;
    }
}
