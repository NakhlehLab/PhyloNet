package edu.rice.cs.bioinfo.library.language.parsing;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/1/11
 * Time: 4:45 PM
 * To change this template use File | Settings | File Templates.
 */
public class CoordinateParseErrorsException extends Exception {

    public final Iterable<CoordinateParseError> Errors;

    public CoordinateParseErrorsException(Iterable<CoordinateParseError> errors)
    {
        Errors = errors;
    }
}
