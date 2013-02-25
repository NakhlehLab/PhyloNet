package edu.rice.cs.bioinfo.library.language.richnewick.reading;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;

import java.io.IOException;
import java.io.InputStream;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/24/12
 * Time: 1:13 PM
 * To change this template use File | Settings | File Templates.
 */
public interface RichNewickParser<N>
{
    N parse(InputStream instream) throws IOException, CoordinateParseErrorsException;
}
