/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.csa;

import edu.rice.cs.bioinfo.library.language.richnewick.reading.csa.CSAError;
import edu.rice.cs.bioinfo.library.programming.Func1Null;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.mockito.Mockito.*;

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
        List<CSAError> result;
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
        result = new ContextAnalyser().analyse(new ArrayList<Object>(), null, new ArrayList<Object>(), null, null, true);
        Assert.assertEquals(0, result.size());

        // one node, all zeros and hybrid type H index 1
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportText(anyObject())).thenReturn("0");
        when(syntaxInspector.getSupportLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result =  new ContextAnalyser().analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null, true);
        Assert.assertEquals(0, result.size());

        // branch length must be a number
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("dogs");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getSupportText(anyObject())).thenReturn("0");
        when(syntaxInspector.getSupportLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = new ContextAnalyser().analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null, true);
        Assert.assertEquals(1, result.size());
        CSAError error = result.get(0);
        Assert.assertTrue(error.Message.contains("dogs"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);

        // support must be a number
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportText(anyObject())).thenReturn("dogs");
        when(syntaxInspector.getSupportLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getSupportColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = new ContextAnalyser().analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null, true);
        Assert.assertEquals(1, result.size());
        error = result.get(0);
        Assert.assertTrue(error.Message.contains("dogs"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);

        // support must be a number between zero and one
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportText(anyObject())).thenReturn("1.5");
        when(syntaxInspector.getSupportLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getSupportColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = new ContextAnalyser().analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null, true);
        Assert.assertEquals(1, result.size());
        error = result.get(0);
        Assert.assertTrue(error.Message.contains("1.5"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);


        // probability must be a number
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportText(anyObject())).thenReturn("0");
        when(syntaxInspector.getSupportLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("dogs");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = new ContextAnalyser().analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null, true);
        Assert.assertEquals(1, result.size());
        error = result.get(0);
        Assert.assertTrue(error.Message.contains("dogs"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);


        // probability must be a number between zero and one
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportText(anyObject())).thenReturn("0");
        when(syntaxInspector.getSupportLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("1.5");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = new ContextAnalyser().analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null, true);
        Assert.assertEquals(1, result.size());
        error = result.get(0);
        Assert.assertTrue(error.Message.contains("1.5"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);


        // hybrid node type must be H, R or LGT
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportText(anyObject())).thenReturn("0");
        when(syntaxInspector.getSupportLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("dogs");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(2);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("1");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(0);

        result = new ContextAnalyser().analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null, true);
        Assert.assertEquals(1, result.size());
        error = result.get(0);
        Assert.assertTrue(error.Message.contains("dogs"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);


        // hybrid node index must be a number
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportText(anyObject())).thenReturn("0");
        when(syntaxInspector.getSupportLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("dogs");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(2);

        result = new ContextAnalyser().analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null, true);
        Assert.assertEquals(1, result.size());
        error = result.get(0);
        Assert.assertTrue(error.Message.contains("dogs"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);

         // hybrid node index must be greater than or equal to one
        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getBranchLengthText(anyObject())).thenReturn("0");
        when(syntaxInspector.getBranchLengthLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getBranchLengthColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportText(anyObject())).thenReturn("0");
        when(syntaxInspector.getSupportLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getSupportColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityText(anyObject())).thenReturn("0");
        when(syntaxInspector.getProbabilityLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getProbabilityColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeType(anyObject())).thenReturn("H");
        when(syntaxInspector.getHybridNodeTypeLineNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeTypeColumnNumber(anyObject())).thenReturn(0);
        when(syntaxInspector.getHybridNodeIndexText(anyObject())).thenReturn("0");
        when(syntaxInspector.getHybridNodeIndexLineNumber(anyObject())).thenReturn(1);
        when(syntaxInspector.getHybridNodeIndexColumnNumber(anyObject())).thenReturn(2);

        result = new ContextAnalyser().analyse(oneNode, syntaxInspector, new ArrayList<Object>(), null, null, true);
        Assert.assertEquals(1, result.size());
        error = result.get(0);
        Assert.assertTrue(error.Message.contains("0"));
        Assert.assertEquals(1, error.LineNumber);
        Assert.assertEquals(2, error.ColumnNumber);

        // an edge between a node and its only parent may be probability one
        networkInspector = mock(NetworkInspector.class);
        when(networkInspector.getAllInEdges(anyObject())).thenReturn(oneEdge);
        when(networkInspector.getEdgeProbabilityText(anyObject())).thenReturn("1");

        result = new ContextAnalyser().analyse(new ArrayList<Object>(), mock(SyntaxNetworkInspector.class), oneNode, networkInspector, Func1Null.Singleton,  true);
        Assert.assertEquals(0, result.size());

        // an edge between a node and its only parent must be probability one
        networkInspector = mock(NetworkInspector.class);
        when(networkInspector.getAllInEdges(anyObject())).thenReturn(oneEdge);
        when(networkInspector.getEdgeProbabilityText(anyObject())).thenReturn(".3");

        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getNodeLabelText(anyObject())).thenReturn(null);

        result = new ContextAnalyser().analyse(new ArrayList<Object>(), syntaxInspector, oneNode, networkInspector, Func1Null.Singleton, true);
        Assert.assertEquals(1, result.size());
        error = result.get(0);
        Assert.assertTrue(error.Message.contains(".3"));

        // probabilities of node's in edges may sum to one
        networkInspector = mock(NetworkInspector.class);
        when(networkInspector.getAllInEdges(anyObject())).thenReturn(twoEdges);
        when(networkInspector.getEdgeProbabilityText(anyObject())).thenReturn(".5");

        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getNodeLabelText(anyObject())).thenReturn(null);

        result = new ContextAnalyser().analyse(new ArrayList<Object>(), syntaxInspector, twoNodes, networkInspector, Func1Null.Singleton, true);
        Assert.assertEquals(0, result.size());

        // probabilities of node's in edges must sum to one
        networkInspector = mock(NetworkInspector.class);
        when(networkInspector.getAllInEdges(anyObject())).thenReturn(twoEdges);
        when(networkInspector.getEdgeProbabilityText(same(twoEdges.get(0)))).thenReturn(".5");
        when(networkInspector.getEdgeProbabilityText(same(twoEdges.get(1)))).thenReturn(".6");

        syntaxInspector = mock(SyntaxNetworkInspector.class);
        when(syntaxInspector.getNodeLabelText(anyObject())).thenReturn(null);

        result = new ContextAnalyser().analyse(new ArrayList<Object>(), syntaxInspector, oneNode, networkInspector, Func1Null.Singleton, true);
        Assert.assertEquals(1, result.size());
        error = result.get(0);
        Assert.assertTrue(error.Message.contains("1.1"));

    }
}
