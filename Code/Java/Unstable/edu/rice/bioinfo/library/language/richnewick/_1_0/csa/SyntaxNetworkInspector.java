package edu.rice.bioinfo.library.language.richnewick._1_0.csa;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 10:48 AM
 * To change this template use File | Settings | File Templates.
 */
public interface SyntaxNetworkInspector<N>
{
    String getNodeLabelText(N node);

    int getNodeLabelTextLineNumber(N node);

    int getNodeLabelTextColumnNumber(N node);

    String getHybridNodeIndexText(N node);

    int getHybridNodeIndexLineNumber(N node);

    int getHybridNodeIndexColumnNumber(N node);

    String getHybridNodeType(N node);

    int getHybridNodeTypeLineNumber(N node);

    int getHybridNodeTypeColumnNumber(N node);

    String getBranchLengthText(N node);

    int getBranchLengthLineNumber(N node);

    int getBranchLengthColumnNumber(N node);

    String getSupportText(N node);

    int getSupportLineNumber(N node);

    int getSupportColumnNumber(N node);

    String getProbabilityText(N node);

    int getProbabilityLineNumber(N node);

    int getProbabilityColumnNumber(N node);
}
