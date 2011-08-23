package edu.rice.bioinfo.library.language.richnewick._1_0;

import java.util.Collections;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/28/11
 * Time: 2:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickReadException extends Exception
{
    public final List<RichNewickReadError> Errors;

    public RichNewickReadException(List<RichNewickReadError> errors)
    {
        Errors = Collections.unmodifiableList(errors);
    }
}
