package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea.acceptancetesting.jung;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.jung.GraphBuilderDirectedSparse;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea.acceptancetesting.ReaSearchFromRichNewick;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func5;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.io.ByteArrayInputStream;
import java.math.BigDecimal;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 3:41 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReaSearchFromRichNewickJung extends ReaSearchFromRichNewick<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>>
{
    private RichNewickReaderAST _reader;

    public ReaSearchFromRichNewickJung(RichNewickReaderAST reader)
    {
        _reader = reader;

    }

    private Func1<String,String> _makeNode = new Func1<String, String>() {
        public String execute(String input) {
            return input;
        }
    };

    private Func5<String, String, BigDecimal, BigDecimal, BigDecimal, PhyloEdge<String>> _makeEdge = new Func5<String, String, BigDecimal, BigDecimal, BigDecimal, PhyloEdge<String>>() {
        public PhyloEdge<String> execute(String source, String destination, BigDecimal branchLength, BigDecimal support, BigDecimal probability) {

           PhyloEdge<String> edge = new PhyloEdge<String>(source, destination);
           if(branchLength != null)
                edge.setBranchLength(branchLength.doubleValue());
           if(support != null)
                edge.setSupport(support.doubleValue());

           if(probability != null)
                edge.setProbability(probability.doubleValue());

            return edge;
        }
    };


    private Func1<PhyloEdge<String>, Tuple<String,String>> _edgeToTuple = new Func1<PhyloEdge<String>, Tuple<String, String>>() {
        public Tuple<String, String> execute(PhyloEdge<String> edge) {

            return new Tuple<String, String>(edge.Source, edge.Destination);
        }
    };

    @Override
    protected DirectedGraphToGraphAdapter<String, PhyloEdge<String>> makeNetwork(String richNewick) {

        try
        {
            GraphBuilderDirectedSparse<String, PhyloEdge<String>> graphBuilder = new GraphBuilderDirectedSparse(_makeNode, _makeEdge );
            RichNewickReadResult<Networks> readResult = _reader.read(new ByteArrayInputStream(richNewick.getBytes()), graphBuilder);

            if(readResult.getContextErrors().iterator().hasNext())
            {
                throw new IllegalArgumentException(readResult.getContextErrors().iterator().next().Message);
            }

            Iterator<NetworkNonEmpty> networks = readResult.getNetworks().Networks.iterator();
            if(networks.hasNext())
            {
                NetworkNonEmpty network = networks.next();

                return new DirectedGraphToGraphAdapter<String, PhyloEdge<String>>(graphBuilder.Graph, new Func1<PhyloEdge<String>, Tuple<String, String>>() {
                    public Tuple<String, String> execute(PhyloEdge<String> edge) {
                        return edge.NodesOfEdge;
                    }
                });
            }

        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }



        return null;

    }
}
