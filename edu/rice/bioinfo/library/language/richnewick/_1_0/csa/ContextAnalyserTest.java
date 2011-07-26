package edu.rice.bioinfo.library.language.richnewick._1_0.csa;
import edu.rice.bioinfo.library.programming.Func1Null;
import org.junit.Assert;
import org.junit.Test;
import org.mockito.Matchers;

import static org.mockito.Mockito.*;

import java.lang.reflect.Array;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 11:07 AM
 * To change this template use File | Settings | File Templates.
 */
public class ContextAnalyserTest
{
    @Test
    public void testAnalyse()
    {
        CSAError[] result;
        SyntaxNetworkInspector syntaxInspector;
        NetworkInspector networkInspector;

        ArrayList<Object> oneNode = new ArrayList<Object>();
        oneNode.add(new Object());

        ArrayList<Object> twoNodes = new ArrayList<Object>();
        twoNodes.add(new Object());
        twoNodes.add(new Object());

        ArrayList<Object> oneEdge = oneNode;
        ArrayList<Object> twoEdges = twoNodes;

        // no errors over no nodes;
        result = ContextAnalyser.Analyse(new ArrayList<Object>(), null, new ArrayList<Object>(), null, null);
        Assert.assertEquals(0, result.length);

        // one node, all zeros and hybrid type H index 1
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBootstrapLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = ContextAnalyser.Analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null);
        Assert.assertEquals(0, result.length);

        // branch length must be a number
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("dogs");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getBootstrapText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBootstrapLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = ContextAnalyser.Analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null);
        Assert.assertEquals(1, result.length);
        CSAError error = result[0];
        Assert.assertTrue(error.Message.contains("dogs"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);

        // bootstrap must be a number
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapText(anyObject())).thenReturn("dogs");
        when(syntaxInspector.getBootstrapLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getBootstrapColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = ContextAnalyser.Analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null);
        Assert.assertEquals(1, result.length);
        error = result[0];
        Assert.assertTrue(error.Message.contains("dogs"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);

        // bootstrap must be a number between zero and one
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapText(anyObject())).thenReturn("1.5");
        when(syntaxInspector.getBootstrapLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getBootstrapColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = ContextAnalyser.Analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null);
        Assert.assertEquals(1, result.length);
        error = result[0];
        Assert.assertTrue(error.Message.contains("1.5"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);


        // probability must be a number
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBootstrapLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("dogs");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = ContextAnalyser.Analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null);
        Assert.assertEquals(1, result.length);
        error = result[0];
        Assert.assertTrue(error.Message.contains("dogs"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);


        // probability must be a number between zero and one
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBootstrapLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("1.5");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = ContextAnalyser.Analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null);
        Assert.assertEquals(1, result.length);
        error = result[0];
        Assert.assertTrue(error.Message.contains("1.5"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);


        // hybrid node type must be H, R or LGT
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBootstrapLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("dogs");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = ContextAnalyser.Analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null);
        Assert.assertEquals(1, result.length);
        error = result[0];
        Assert.assertTrue(error.Message.contains("dogs"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);


        // hybrid node index must be a number
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBootstrapLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("dogs");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(2);

        result = ContextAnalyser.Analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null);
        Assert.assertEquals(1, result.length);
        error = result[0];
        Assert.assertTrue(error.Message.contains("dogs"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);

         // hybrid node index must be greater than or equal to one
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBootstrapLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBootstrapColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("0");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(2);

        result = ContextAnalyser.Analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null);
        Assert.assertEquals(1, result.length);
        error = result[0];
        Assert.assertTrue(error.Message.contains("0"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);

        // an edge between a node and its only parent may be probability one
        networkInspector = mock(NetworkInspector.class);
        when(networkInspector.getAllInEdges(anyObject())).thenReturn(oneEdge);
        when(networkInspector.getProbabilityText(anyObject())).thenReturn("1");

        result = ContextAnalyser.Analyse(new ArrayList<Object>(), null, oneNode, networkInspector, null);
        Assert.assertEquals(0, result.length);

        // an edge between a node and its only parent must be probability one
        networkInspector = mock(NetworkInspector.class);
        when(networkInspector.getAllInEdges(anyObject())).thenReturn(oneEdge);
        when(networkInspector.getProbabilityText(anyObject())).thenReturn(".3");

        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getNodeLabelText(anyObject())).thenReturn(null);

        result = ContextAnalyser.Analyse(new ArrayList<Object>(), syntaxInspector, oneNode, networkInspector, Func1Null.Singleton);
        Assert.assertEquals(1, result.length);
        error = result[0];
        Assert.assertTrue(error.Message.contains(".3"));

        // probabilities of node's in edges may sum to one
        networkInspector = mock(NetworkInspector.class);
        when(networkInspector.getAllInEdges(anyObject())).thenReturn(twoEdges);
        when(networkInspector.getProbabilityText(anyObject())).thenReturn(".5");

        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getNodeLabelText(anyObject())).thenReturn(null);

        result = ContextAnalyser.Analyse(new ArrayList<Object>(), syntaxInspector, twoNodes, networkInspector, Func1Null.Singleton);
        Assert.assertEquals(0, result.length);

    }
}
