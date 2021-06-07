package edu.rice.cs.bioinfo.library.language.vaal.out.reading;

import edu.rice.cs.bioinfo.library.language.vaal.out.reading.csa.CSAError;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 3/21/13
 * Time: 2:01 PM
 * To change this template use File | Settings | File Templates.
 */
public interface OutReadResult<T>
{
    T getDifferences();

    Iterable<CSAError> getContextErrors();
}
