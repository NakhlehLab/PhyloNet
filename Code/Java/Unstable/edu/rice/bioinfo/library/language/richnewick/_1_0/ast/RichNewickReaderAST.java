package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReadError;
import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReadException;
import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReader;
import edu.rice.bioinfo.library.language.richnewick._1_0.csa.CSAError;
import edu.rice.bioinfo.library.language.richnewick._1_0.csa.ContextAnalyser;
import edu.rice.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilder;
import edu.rice.bioinfo.library.programming.Func1;

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
public abstract class RichNewickReaderAST implements RichNewickReader<Networks>
{
    public <N> Networks read(InputStream instream, GraphBuilder<N> graphBuilder) throws RichNewickReadException {
        Networks networks = parse(instream);


        for(Network network : networks.Networks)
        {
            final ASTNetworkInspector inspector = new ASTNetworkInspector(network);
            Func1<Object,NetworkInfo> networkNodeToPrimarySyntaxNode = new Func1<Object, NetworkInfo>() {
                public NetworkInfo execute(Object networkNode) {

                    return inspector.getPrimarySyntaxNode(networkNode);
                }
            };

            List<CSAError> errors = ContextAnalyser.Analyse(inspector.getSyntaxNodes(), inspector, inspector.getNetworkNodes(), inspector, networkNodeToPrimarySyntaxNode);

            if(errors.size() > 0)
            {
                LinkedList<RichNewickReadError> readErrors = new LinkedList<RichNewickReadError>();
                for(CSAError csaError : errors)
                {
                    readErrors.add(new RichNewickReadError(csaError.Message, csaError.LineNumber, csaError.ColumnNumber));
                }
                throw new RichNewickReadException(readErrors);
            }

            DAGFactory.makeDAG(network, graphBuilder);
        }

        return networks;
    }

    protected abstract Networks parse(InputStream instream) throws RichNewickReadException;
}
