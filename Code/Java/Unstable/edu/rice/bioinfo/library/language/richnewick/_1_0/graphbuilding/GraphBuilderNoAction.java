package edu.rice.bioinfo.library.language.richnewick._1_0.graphbuilding;

import edu.rice.bioinfo.library.language.richnewick._1_0.HybridNodeType;

import java.math.BigDecimal;
import java.math.BigInteger;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/1/11
 * Time: 5:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class GraphBuilderNoAction implements GraphBuilder<Object>
{
    public static final GraphBuilderNoAction Singleton = new GraphBuilderNoAction();

    private GraphBuilderNoAction()
    {
    }

    public Object createNode(String label)
    {
        return null;
    }

    public Object createHybridNode(String label, HybridNodeType hybridType, BigInteger hybridNodeIndex)
    {
        return null;
    }

    public void createDirectedEdge(Object tail, Object tip, BigDecimal branchLength, BigDecimal support, BigDecimal probability)
    {

    }
}
