package edu.rice.bioinfo.programs.richnewickwebview;

import com.mxgraph.view.mxGraph;
import edu.rice.bioinfo.library.language.richnewick._1_0.HybridNodeType;
import edu.rice.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilder;

import java.math.BigDecimal;
import java.math.BigInteger;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/27/11
 * Time: 2:11 PM
 * To change this template use File | Settings | File Templates.
 */
class MxGraphBuilder implements GraphBuilder<Object>{

    private mxGraph _graph = new mxGraph();

    private Object _parent = _graph.getDefaultParent();

    MxGraphBuilder()
    {
        _graph.getModel().beginUpdate();
    }

    public Object createNode(String label)
    {

            return _graph.insertVertex(_parent, null, label, 20, 20, 80, 30);

    }

    public Object createHybridNode(String label, HybridNodeType hybridNodeType, BigInteger bigInteger) {
        return _graph.insertVertex(_parent, null, label, 20, 20, 80, 30);
    }

    public void createDirectedEdge(Object tail, Object tip, BigDecimal bigDecimal, BigDecimal bigDecimal1, BigDecimal bigDecimal2) {
        _graph.insertEdge(_parent, null, "Edge", tail, tip);
    }

    public mxGraph getGraph()
    {
        _graph.getModel().endUpdate();
        return _graph;
    }
}
