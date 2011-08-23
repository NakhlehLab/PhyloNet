package edu.rice.bioinfo.programs.richnewickwebview;

import com.sun.org.apache.bcel.internal.generic.NEW;
import edu.rice.bioinfo.library.language.richnewick._1_0.HybridNodeType;
import edu.rice.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilder;
import edu.rice.bioinfo.library.programming.Func;
import edu.rice.bioinfo.library.programming.Func1;
import edu.uci.ics.jung.algorithms.layout.CircleLayout;
import edu.uci.ics.jung.algorithms.layout.DAGLayout;
import edu.uci.ics.jung.algorithms.layout.KKLayout;
import edu.uci.ics.jung.algorithms.layout.TreeLayout;
import edu.uci.ics.jung.algorithms.layout.Layout;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseMultigraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;
import edu.uci.ics.jung.visualization.BasicVisualizationServer;
import edu.uci.ics.jung.visualization.decorators.ToStringLabeller;
import edu.uci.ics.jung.visualization.renderers.Renderer;

import javax.swing.text.Position;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collection;
import java.util.concurrent.BrokenBarrierException;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/27/11
 * Time: 3:39 PM
 * To change this template use File | Settings | File Templates.
 */
public class JungGraphBuilder implements GraphBuilder<Object> {

    private DirectedGraph<Object,Object> _dag = new DirectedSparseMultigraph<Object,Object>();

    private Func1<DirectedGraph<Object,Object>,Layout<Object,Object>> _makeLayout;

    JungGraphBuilder(Func1<DirectedGraph<Object,Object>,Layout<Object,Object>> makeLayout)
    {
       _makeLayout = makeLayout;
    }

    public Object createNode(final String label) {

        Object node = new Object()
        {
            @Override
            public String toString()
            {
                return label == null ? "" : label;
            }

        };
        _dag.addVertex(node);
        return node;
    }

    public Object createHybridNode(final String label, final HybridNodeType hybridType, final BigInteger hybridNodeIndex) {
        Object node = new Object()
        {
            @Override
            public String toString()
            {
                String stringLabel = "(" + hybridNodeIndex + (makeHybridTypeLabel(hybridType)) + ")";

                if(label != null)
                    stringLabel = label + " " + stringLabel;

                return stringLabel;
            }

        };
        _dag.addVertex(node);
        return node;
    }

    private String makeHybridTypeLabel(HybridNodeType hybridType) {

        if(hybridType == null || hybridType == HybridNodeType.Unspecified)
            return "";

        if(hybridType == HybridNodeType.LateralGeneTransfer)
            return "LGT";

        if(hybridType == HybridNodeType.Hybridization)
            return "H";

        if(hybridType == HybridNodeType.Recombination)
            return "R";

         throw new IllegalArgumentException("Unexpected hybrid type.");
    }

    public void createDirectedEdge(Object tail, Object tip, BigDecimal branchLength, BigDecimal bootstrap, BigDecimal probability) {

         String label = "";

        if(branchLength != null || bootstrap != null || probability != null)
        {

            ArrayList<String> adj = new ArrayList<String>();

            if(branchLength != null)
                adj.add("bl:" + branchLength.toPlainString());

            if(bootstrap != null)
                adj.add("bs:" + bootstrap.toPlainString());

            if(probability != null)
                adj.add("p:" + probability);

            for(int i = 0; i<adj.size(); i++)
            {
                label+= adj.get(i).replace('"', '\'');

                if(i != adj.size() -1)
                    label+=" ";
            }


        }

        final String finalString = label;
        _dag.addEdge(new Object()
        {
            @Override
            public String toString()
            {
                return finalString;
            }
        }, tail, tip, EdgeType.DIRECTED);
    }

    public BasicVisualizationServer<Object,Object> getVisualization()
    {
        BasicVisualizationServer<Object, Object> vv = new BasicVisualizationServer<Object, Object>(_makeLayout.execute(_dag));
        vv.getRenderContext().setVertexLabelTransformer(new ToStringLabeller());
        vv.getRenderContext().setEdgeLabelTransformer(new ToStringLabeller());
        return vv;
    }
}
