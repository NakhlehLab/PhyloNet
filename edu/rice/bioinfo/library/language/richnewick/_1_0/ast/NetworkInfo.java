package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 4:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkInfo implements AbstractSyntaxNode {

    public final NodeLabel NodeLabel;

    public final HybridNodeQualifier HybridNodeQualifier;

    public final BranchLength BranchLength;

    public final Support Support;

    public final Probability Probability;

    public NetworkInfo(NodeLabel nodeLabel, HybridNodeQualifier qualifier, BranchLength branchLength,
                       Support support, Probability probability)
    {
        NodeLabel = nodeLabel;
        HybridNodeQualifier = qualifier;
        BranchLength = branchLength;
        Support = support;
        Probability = probability;
    }
}
