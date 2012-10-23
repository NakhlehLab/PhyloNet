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

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;

import java.io.*;
import java.util.*;
import java.util.zip.DeflaterOutputStream;


/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 3/13/12
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class AssignBranchLengthsMaxGTProb extends CommandBaseFileOut{
    private HashMap<String,String> _taxonMap = null;
    private boolean  _printDetail = false;
    private NetworkNonEmpty _speciesNetwork;
    private List<NetworkNonEmpty> _geneTrees;
    private ParameterIdentList _geneTreeParam;
    private double _maxBranchLength;
    private int _maxAssigmentAttemptsPerBranchParam;
    private int _assigmentRounds = -1;

    public AssignBranchLengthsMaxGTProb(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                                        Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    @Override
    protected int getMinNumParams(){
        return 4;
    }

    @Override
    protected int getMaxNumParams(){
        return 8;
    }

    @Override
    protected boolean checkParamsForCommand(){
        boolean noError = true;

        _speciesNetwork = this.assertAndGetNetwork(0);
        noError = noError && _speciesNetwork != null;

        ParameterIdent maxBranchLengthParam = this.assertParameterIdent(1);
        noError = noError && maxBranchLengthParam != null;

        if(maxBranchLengthParam != null)
        {
            _maxBranchLength = Double.parseDouble(maxBranchLengthParam.Content);
        }

        ParameterIdent maxAssigmentAttemptsPerBranchParam = this.assertParameterIdent(2);
        noError = noError && maxAssigmentAttemptsPerBranchParam != null;

        if(maxAssigmentAttemptsPerBranchParam != null)
        {
            _maxAssigmentAttemptsPerBranchParam = Integer.parseInt(maxAssigmentAttemptsPerBranchParam.Content);
        }

        ParameterIdent assignmentRoundsParam = this.assertParameterIdent(3);
        noError = noError && assignmentRoundsParam  != null;

        if(assignmentRoundsParam != null)
        {
            _assigmentRounds = Integer.parseInt(assignmentRoundsParam.Content);
        }

        _geneTreeParam = this.assertParameterIdentList(4);
        noError = noError && _geneTreeParam != null;
        _geneTrees = new LinkedList<NetworkNonEmpty>();
        for(String ident : _geneTreeParam.Elements)
        {
            noError = noError && this.assertNetworkExists(ident, _geneTreeParam.getLine(), _geneTreeParam.getColumn());
            if(noError)
            {
                _geneTrees.add(this.sourceIdentToNetwork.get(ident));
            }
        }

        ParamExtractorAllelMap aParam = new ParamExtractorAllelMap("a", this.params, this.errorDetected);
        if(aParam.ContainsSwitch){
            noError = noError && aParam.IsValidMap;
            if(aParam.IsValidMap){
                _taxonMap = aParam.ValueMap;
            }
        }

        ParamExtractor pParam = new ParamExtractor("p", this.params, this.errorDetected);
        if(pParam.ContainsSwitch)
        {
            _printDetail = true;
        }

        noError = noError && checkForUnknownSwitches("p", "a");
        checkAndSetOutFile(aParam, pParam);

        return  noError;
    }

    @Override
    protected String produceResult() {
        StringBuffer result = new StringBuffer();
        
        final List<Tree> geneTrees = new ArrayList<Tree>();
        final List<Integer> counter = new ArrayList<Integer>();
        for(NetworkNonEmpty geneTree : _geneTrees){
            String phylonetGeneTree = NetworkTransformer.toENewickTree(geneTree);
            NewickReader nr = new NewickReader(new StringReader(phylonetGeneTree));
            STITree<Double> newtr = new STITree<Double>(true);
            try
            {
                nr.readTree(newtr);
            }
            catch(Exception e)
            {
                errorDetected.execute(e.getMessage(),
                        this._motivatingCommand.getLine(), this._motivatingCommand.getColumn());
            }
            boolean found = false;
            int index = 0;
            for(Tree tr: geneTrees){
                if(Trees.haveSameRootedTopology(tr, newtr)){
                    found = true;
                    break;
                }
                index++;
            }
            if(found){
                counter.set(index, counter.get(index)+1);
            }
            else{
                geneTrees.add(newtr);
                counter.add(1);
            }
        }

        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        final Network<Double> speciesNetwork = transformer.makeNetwork(_speciesNetwork);

       // int gtIndex = 0;


            for(NetNode<Double> parent : speciesNetwork.dfs())  // make the branch length of every edge initially  1
            {
                for(NetNode<Double> child : parent.getChildren())
                {
                    double initialBL = child.getParentDistance(parent);
                    if(initialBL == NetNode.NO_DISTANCE || Double.isNaN(initialBL))
                        initialBL = 1.0;
                    child.setParentDistance(parent, initialBL);
                }
            }

        for(int assigmentRound = 0; assigmentRound <_assigmentRounds; assigmentRound++)
        {
            for(final NetNode<Double> parent : speciesNetwork.bfs())  // for each edge, find best branch length
            {
                for(final NetNode<Double> child : parent.getChildren())
                {
                    final Container<Double> bestFoundBranchLength = new Container<Double>(null);
                    final Container<Double> bestFoundBranchLengthCorrespondingGTProb = new Container<Double>(null);
                    UnivariateFunction functionToOptimize = new UnivariateFunction() {
                        public double value(double suggestedBranchLength) {


                            child.setParentDistance(parent, suggestedBranchLength);

                            double prob = computeGTProb(speciesNetwork, geneTrees, counter);

                            if(bestFoundBranchLengthCorrespondingGTProb.getContents() == null || bestFoundBranchLengthCorrespondingGTProb.getContents() < prob)
                            {
                                bestFoundBranchLength.setContents(suggestedBranchLength);
                                bestFoundBranchLengthCorrespondingGTProb.setContents(prob);
                            }
                            return prob;
                        }
                    };
                    BrentOptimizer optimizer = new BrentOptimizer(.000000000001,.0000000000000001);
                    double initialBL = child.getParentDistance(parent);

                    try
                    {
                        UnivariatePointValuePair maxFoundValue = optimizer.optimize(_maxAssigmentAttemptsPerBranchParam, functionToOptimize, GoalType.MAXIMIZE, Double.MIN_VALUE, _maxBranchLength, initialBL);
                    }
                    catch(TooManyEvaluationsException e)
                    {
                    }
                    child.setParentDistance(parent, bestFoundBranchLength.getContents());

                }
            }
        }

         RnNewickPrinter<Double> rnNewickPrinter = new RnNewickPrinter<Double>();
         StringWriter sw = new StringWriter();
         rnNewickPrinter.print(speciesNetwork, sw);

        double gtProb = computeGTProb(speciesNetwork, geneTrees, counter);
        result.append("\nTotal log probability: " + gtProb + ": " + sw.toString());

        return result.toString();

    }

    private double computeGTProb(Network<Double> speciesNetwork, List<Tree> geneTrees, List<Integer> counter)
    {
        GeneTreeProbability gtp = new GeneTreeProbability();
        Iterator<Double> probList = gtp.calculateGTDistribution(speciesNetwork, geneTrees, _taxonMap, _printDetail).iterator();
        Iterator<Integer> counterIt = counter.iterator();
        double total = 0;

        for(Tree gt: geneTrees){
            double prob = probList.next();
            int count = counterIt.next();
            total += Math.log(prob)*count;
        }

        if(Double.isNaN(total))
        {
            throw new RuntimeException();
        }

        return total;
    }
}
