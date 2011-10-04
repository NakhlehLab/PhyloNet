package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 5:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickString implements PySONNode{

    public final String String;

    public final int LineNumber;

    public final int ColumnNumber;

    public RichNewickString(String string, int lineNumber, int columnNumber)
    {
        LineNumber = lineNumber;
        ColumnNumber = columnNumber;
        String = string;
    }
}
