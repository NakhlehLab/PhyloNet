package edu.rice.cs.bioinfo.library.language.richnewick._1_0.graphbuilding;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.security.PublicKey;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/26/11
 * Time: 10:52 AM
 * To change this template use File | Settings | File Templates.
 */
public interface GraphBuilder<N>
{
    public N createNode(String label);

    public N createHybridNode(String label, HybridNodeType hybridType, BigInteger hybridNodeIndex);

    public void createDirectedEdge(N tail, N tip, BigDecimal branchLength, BigDecimal support, BigDecimal probability);
}
