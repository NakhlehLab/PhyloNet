package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorDefault;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.RichNewickReaderBase;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.csa.ASTContextAnalyser;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.csa.CSAError;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilder;
import edu.rice.cs.bioinfo.library.programming.Func1;

import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/28/11
 * Time: 10:20 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class RichNewickReaderAST extends RichNewickReaderBase<Networks>
{
    public <N> RichNewickReadResult<Networks> read(InputStream instream, GraphBuilder<N> graphBuilder) throws IOException, CoordinateParseErrorsException {

        final Networks networks = parse(instream);

        final List<CSAError> errors = new LinkedList<CSAError>();

        for(Network network : networks.Networks)
        {

            for(CSAError error : ASTContextAnalyser.analyse(network))
            {
                errors.add(error);
            }

            if(errors.size() > 0)
            {
                LinkedList<CoordinateParseError> readErrors = new LinkedList<CoordinateParseError>();
                for(CSAError csaError : errors)
                {
                    readErrors.add(new CoordinateParseErrorDefault(csaError.Message, csaError.LineNumber, csaError.ColumnNumber) {
                    });
                }
                throw new CoordinateParseErrorsException(readErrors);
            }

            DAGFactory.makeDAG(network, graphBuilder);
        }

        return new RichNewickReadResult<Networks>() {
            public Networks getNetworks() {
                return networks;
            }

            public Iterable<CSAError> getContextErrors() {
                return errors;  //To change body of implemented methods use File | Settings | File Templates.
            }
        };
    }

    protected abstract Networks parse(InputStream instream) throws CoordinateParseErrorsException, IOException;
}
