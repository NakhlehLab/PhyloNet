package edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.jung.GraphBuilderDirectedSparse;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.MDCOnNetworkYFFromRichNewick;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func5;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.io.ByteArrayInputStream;
import java.math.BigDecimal;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/14/12
 * Time: 2:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class MDCOnNetworkYFFromRichNewickJung extends MDCOnNetworkYFFromRichNewick<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>>
{
    private RichNewickReaderAST _reader;

    public MDCOnNetworkYFFromRichNewickJung(RichNewickReaderAST reader)
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



    @Override
    protected DirectedGraphToGraphAdapter<String, PhyloEdge<String>> makeNetwork(String richNewick) {

          GraphBuilderDirectedSparse<String, PhyloEdge<String>> graphBuilder = new GraphBuilderDirectedSparse(_makeNode, _makeEdge );
          _reader.readAnyErrorToRuntimeException(new ByteArrayInputStream(richNewick.getBytes()), graphBuilder);
          return new DirectedGraphToGraphAdapter<String, PhyloEdge<String>>(graphBuilder.Graph, new Func1<PhyloEdge<String>, Tuple<String, String>>() {
              public Tuple<String, String> execute(PhyloEdge<String> edge) {
                  return edge.NodesOfEdge;
              }
          });
    }
}
