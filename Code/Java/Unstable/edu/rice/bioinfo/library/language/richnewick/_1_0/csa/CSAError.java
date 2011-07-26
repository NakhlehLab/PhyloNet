package edu.rice.bioinfo.library.language.richnewick._1_0.csa;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/18/11
 * Time: 5:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class CSAError
{
    public final String Message;

    public final int LineNumber;

    public final int ColumnNumber;

    public CSAError(String message, int lineNumber, int columnNumber)
    {
        Message = message;
        LineNumber = lineNumber;
        ColumnNumber = columnNumber;
    }
}
