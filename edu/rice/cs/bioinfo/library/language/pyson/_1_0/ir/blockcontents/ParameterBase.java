package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/3/11
 * Time: 6:49 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ParameterBase implements Parameter
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

    public ParameterBase(int line, int column)
    {
        _line = line;
        _column = column;
    }
}
