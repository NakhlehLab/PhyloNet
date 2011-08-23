package edu.rice.bioinfo.library.language.richnewick._1_0;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/28/11
 * Time: 2:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickReadError
{
    public final String Message;

    public final Integer LineNumber;

    public final Integer ColumnNumber;

    public RichNewickReadError(String message, Integer lineNumber, Integer columnNumber)
    {
        Message = message;
        LineNumber = lineNumber;
        ColumnNumber = columnNumber;
    }
}
