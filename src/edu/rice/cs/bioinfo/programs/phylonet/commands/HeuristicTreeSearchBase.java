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

package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.IsRooted;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.DeepCoalescencesCounter;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseGraph;

import java.io.IOException;
import java.io.StringReader;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 4/12/12
 * Time: 3:15 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class HeuristicTreeSearchBase extends CommandBaseFileOut
{
    protected Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> inputTree;

    private final Random _rand;

     private double _bootstrap = 1.0;

    protected int numIters;

    protected Func1<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,Double> scoreTree;

    protected Func2<Double,Double,Boolean> isScoreBetter;

    private List<Tree> _geneTrees = new ArrayList<Tree>();

    private HashMap<String,String> _taxonMap = null;

    private  boolean isRooted = false;

    protected
    Predicate1<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>> isTreeRooted = new Predicate1<Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]>>() {
            public boolean execute(Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]> labelGraph) {
                return isRooted;
            }
        };

    private  Func1<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,Double> _calGTProbScoreTree = new Func1<Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]>, Double>() {
        public Double execute(Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]> tree) {

        Network<Double> network = toNetwork(tree);

	// kliu - example on how to call GeneTreeProbability

        GeneTreeProbability gtp = new GeneTreeProbability();
        Iterator<Double> probList = gtp.calculateGTDistribution(network, _geneTrees, _taxonMap, false).iterator();
     //   result.append("\n");
        double total = 0;
        for(Tree gt: _geneTrees){
            double prob = probList.next();
	    // kliu - why no log calculations in GeneTreeProbability?
	    // base 10 unnecessary
            total += Math.log10(prob);
        }

           return Math.pow(10, total);

        }
    };

    private  Func1<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,Double> _mdcScoreTree = new Func1<Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]>, Double>()
    {

        public Double execute(Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]> tree) {
            Tree network = toTree(tree);
            List<MutableTuple<Tree,Double>> _gtWithWeight = new ArrayList<MutableTuple<Tree, Double>>();
            for(Tree tr: _geneTrees){
                _gtWithWeight.add(new MutableTuple<Tree, Double>(tr, 1.0));
            }
            return (double) DeepCoalescencesCounter.countExtraCoal(_gtWithWeight, network, false, _bootstrap);
        }
    };

    private Tree toTree(Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]> tree) {

        try
        {
            String rn = JUNGToRN.toRichNewick(tree, isRooted);
            return new NewickReader(new StringReader(rn)).readTree();
        }
        catch (Exception e)
        {
            throw new RuntimeException(e);
        }

    }

    private Network<Double> toNetwork(Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]> tree) {

        try
        {
            String rn = JUNGToRN.toRichNewick(tree, isRooted);
            rn = rn.replace("[&U]", "");
            return new ExNewickReader<Double>(new StringReader(rn)).readNetwork();
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }

    }

     public HeuristicTreeSearchBase(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                                Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader, Random rand) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
        _rand = rand;
    }

    @Override
    protected int getMinNumParams() {
        return 4;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 5;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = true;

        Parameter startingInput = this.params.get(0);

        inputTree = startingInput.execute(new ParameterAlgo<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>, Object, RuntimeException>() {
            public Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> forIdentifier(ParameterIdent parameterIdent, Object o) throws RuntimeException {
                final NetworkNonEmpty nexusTree = assertAndGetNetwork(0);

                isRooted = nexusTree.RootageQualifier.execute(new IsRooted(), null);

                if(nexusTree == null)
                {
                    return null;
                }

                return new NetworkToJUNG().toTreeWithoutEdgeProperties(nexusTree);

            }

            public Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> forIdentList(ParameterIdentList parameterIdentList, Object o) throws RuntimeException {
                assertAndGetNetwork(0);
                return null;
            }

            public Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> forQuote(ParameterQuote parameterQuote, Object o) throws RuntimeException {
                assertAndGetNetwork(0);
                return null;
            }

            public Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Object o) throws RuntimeException {
                assertAndGetNetwork(0);
                return null;
            }

            public Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> forIdentSet(ParameterIdentSet parameterIdentSet, Object o) throws RuntimeException {

                inputTree = new SparseGraph<NetworkToJUNG.Label, NetworkToJUNG.Label[]>();

                ArrayList<NetworkToJUNG.Label> nodes = new ArrayList<NetworkToJUNG.Label>();


                for(String taxonName : parameterIdentSet.Elements)
                {
                    nodes.add(new NetworkToJUNG.Label(taxonName));
                }

                while(nodes.size() > 2)
                {
                    Collections.shuffle(nodes);
                    NetworkToJUNG.Label c1 = nodes.remove(0);
                    NetworkToJUNG.Label c2 = nodes.remove(0);
                    NetworkToJUNG.Label parent = new NetworkToJUNG.Label(null);
                    inputTree.addVertex(parent);
                    inputTree.addEdge(new NetworkToJUNG.Label[] { parent, c1 }, parent, c1);
                    inputTree.addEdge(new NetworkToJUNG.Label[] { parent, c2 }, parent, c2);
                    nodes.add(parent);

                }

                inputTree.addEdge(new NetworkToJUNG.Label[] { nodes.get(0), nodes.get(1) }, nodes.get(0), nodes.get(1));

                return inputTree;

            }

            public Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> forTaxaMap(ParameterTaxaMap parameterTaxaMap, Object o) throws RuntimeException {
                assertAndGetNetwork(0);
                return null;
            }
        }, null);

        noError = noError && inputTree != null;

        ParameterIdentSet geneTreeParam = this.assertParameterIdentSet(1);
        noError = noError && geneTreeParam != null;
        List<NetworkNonEmpty> geneTrees = new LinkedList<NetworkNonEmpty>();
        for(String ident : geneTreeParam.Elements)
        {
            noError = noError && this.assertNetworkExists(ident, geneTreeParam.getLine(), geneTreeParam.getColumn());
            if(noError)
            {
                geneTrees.add(this.sourceIdentToNetwork.get(ident));
            }
        }

         for(NetworkNonEmpty geneTree : geneTrees){
            String phylonetGeneTree = NetworkTransformer.toENewickTree(geneTree);
            NewickReader nr = new NewickReader(new StringReader(phylonetGeneTree));
            STITree<Double> gt = new STITree<Double>(true);
            try
            {
                nr.readTree(gt);
                _geneTrees.add(gt);
            }
            catch(Exception e)
            {
                errorDetected.execute(e.getMessage(),
                        this._motivatingCommand.getLine(), this._motivatingCommand.getColumn());
            }
        }


        ParameterIdent methodParam = this.assertParameterIdent(2);
        noError = noError && methodParam != null;

        if(methodParam != null)
        {
            if(methodParam.Content.toLowerCase().equals("gtprob"))
            {
                scoreTree = _calGTProbScoreTree;
                isScoreBetter = new Func2<Double, Double, Boolean>() {
                    public Boolean execute(Double aDouble, Double aDouble1) {
                        return aDouble > aDouble1;
                    }
                };
            }
            else if(methodParam.Content.toLowerCase().equals("mdc"))
            {
                scoreTree = _mdcScoreTree;
                isScoreBetter = new Func2<Double, Double, Boolean>() {
                    public Boolean execute(Double aDouble, Double aDouble1) {
                        return aDouble < aDouble1;
                    }
                };

            }
            else
            {
                errorDetected.execute("Unknown scoring method.", methodParam.getLine(), methodParam.getColumn());
            }
        }

        ParameterIdent numInterchangesParam = this.assertParameterIdent(3);
        noError = noError && numInterchangesParam != null;


        if(noError)
        {
            try
            {
                numIters = Integer.parseInt(numInterchangesParam.Content);
            }
            catch(NumberFormatException e)
            {
                errorDetected.execute("Interchange count must be an integer.  Found: " + numInterchangesParam.Content, numInterchangesParam.getLine(), numInterchangesParam.getColumn());
            }
        }

        ParamExtractorAllelMap aParam = new ParamExtractorAllelMap("a", this.params, this.errorDetected);
        if(aParam.ContainsSwitch){
            noError = noError && aParam.IsValidMap;
            if(aParam.IsValidMap){
                _taxonMap = aParam.ValueMap;
            }
        }

        this.checkForUnknownSwitches("a");
        this.checkAndSetOutFile(aParam);

        return noError;
    }
}
