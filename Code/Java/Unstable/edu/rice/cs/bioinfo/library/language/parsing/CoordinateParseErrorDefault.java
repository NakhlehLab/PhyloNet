/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/1/11
 * Time: 2:13 PM
 * To change this template use File | Settings | File Templates.
 */

package edu.rice.cs.bioinfo.library.language.parsing;

public class CoordinateParseErrorDefault implements CoordinateParseError {

    private final String _message;

    private final int _lineNumber, _columnNumber;

    public CoordinateParseErrorDefault(String message, int lineNumber, int columnNumber)
    {
        _message = message;
        _lineNumber = lineNumber;
        _columnNumber = columnNumber;
    }

    public String getMessage() {
        return _message;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getLineNumber() {
        return _lineNumber;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getColumnNumber() {
        return _columnNumber;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
