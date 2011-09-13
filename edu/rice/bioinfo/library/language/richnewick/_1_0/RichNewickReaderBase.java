package edu.rice.bioinfo.library.language.richnewick._1_0;

import edu.rice.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilder;
import edu.rice.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilderNoAction;

import java.io.IOException;
import java.io.InputStream;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/8/11
 * Time: 2:05 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class RichNewickReaderBase<T> implements RichNewickReader<T> {

    public <N> RichNewickReadResult<T> read(InputStream instream) throws CoordinateParseErrorsException, IOException {
        return read(instream, GraphBuilderNoAction.Singleton);
    }

    public abstract <N> RichNewickReadResult<T> read(InputStream instream, GraphBuilder<N> graphBuilder) throws CoordinateParseErrorsException, IOException;
}
