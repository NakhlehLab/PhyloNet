package edu.rice.cs.bioinfo.library.language.richnewick._1_0;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.csa.CSAError;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/8/11
 * Time: 2:03 PM
 * To change this template use File | Settings | File Templates.
 */
public interface RichNewickReadResult<N> {

    public N getNetworks();

    public Iterable<CSAError> getContextErrors();
}
