package edu.rice.bioinfo.library.language.richnewick._1_0;

import edu.rice.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilder;

import java.io.InputStream;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/28/11
 * Time: 10:15 AM
 * To change this template use File | Settings | File Templates.
 */
public interface RichNewickReader<T>
{
    public <N> T read(InputStream instream, GraphBuilder<N> graphBuilder) throws RichNewickReadException;
}
