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
import edu.rice.cs.bioinfo.library.phylogenetics.IsLeaf;
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

import javax.management.openmbean.InvalidOpenTypeException;
import java.io.*;
import java.math.BigDecimal;
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
    private int _maxAssigmentAttemptsPerBranchParam = -1;
    private int _assigmentRounds = Integer.MAX_VALUE;
    private double _improvementThreshold;

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

        ParameterIdent improvementThresholdParam = this.assertParameterIdent(2);
        noError = noError && improvementThresholdParam != null;

        if(improvementThresholdParam != null)
        {
            _improvementThreshold = Double.parseDouble(improvementThresholdParam.Content);
        }


        ParameterIdent maxAssigmentAttemptsPerBranchParam = this.assertParameterIdent(3);
        noError = noError && maxAssigmentAttemptsPerBranchParam != null;

        if(maxAssigmentAttemptsPerBranchParam != null)
        {
            _maxAssigmentAttemptsPerBranchParam = Integer.parseInt(maxAssigmentAttemptsPerBranchParam.Content);
        }

        /*
        ParameterIdent assignmentRoundsParam = this.assertParameterIdent(3);
        noError = noError && assignmentRoundsParam  != null;

        if(assignmentRoundsParam != null)
        {
            _assigmentRounds = Integer.parseInt(assignmentRoundsParam.Content);
        } */

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

        /*
         * Make the branch length of every edge initially  1 if no initial user value is specified.
         * Make the hybrid prob of every hybrid edge initially .5 if no initial user value is specified.
         */
        for(NetNode<Double> parent : speciesNetwork.bfs())
        {
            for(NetNode<Double> child : parent.getChildren())
            {

                double initialBL = child.getParentDistance(parent);
                if(initialBL == NetNode.NO_DISTANCE || Double.isNaN(initialBL)) // no specification from user  on branch length
                    initialBL = 1.0;
                child.setParentDistance(parent, initialBL);

                if(child.getParentNumber() == 2)
                {
                    for(NetNode<Double> hybridParent : child.getParents())
                    {
                        if(child.getParentProbability(hybridParent) == 1.0) // no specification from user  on hybrid prob
                        {
                            child.setParentProbability(parent, 0.5);
                        }
                    }
                }
            }
        }

        /*
         * Try to assign branch lengths and hybrid probs to increase GTProb from the initial network.
         * Except branch lengths of leaf edges.  They don't impact GTProb.
         */

        // def: a round is an attempt to tweak each branch length and each hybrid prob.
        boolean continueRounds = true; // keep trying to improve network
        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(computeGTProb(speciesNetwork, geneTrees, counter));  // records the GTProb of the network at all times

        int assigmentRound = 0;
        for(; assigmentRound <_assigmentRounds && continueRounds; assigmentRound++)
        {
            System.out.println("\nround " + assigmentRound);
            /*
             * Prepare a random ordering of network edge examinations each of which attempts to change a branch length or hybrid prob to improve the GTProb score.
             */

            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();
            List<Proc> assigmentActions = new ArrayList<Proc>(); // store adjustment commands here.  Will execute them one by one later.

            // add branch length adjustments to the list
            for(final NetNode<Double> parent : speciesNetwork.bfs())
            {
                for(final NetNode<Double> child : parent.getChildren())
                {
                    if(child.isLeaf()) // leaf edge, skip
                        continue;

                    assigmentActions.add(new Proc()
                    {
                        public void execute()
                        {
                            UnivariateFunction functionToOptimize = new UnivariateFunction() {
                                public double value(double suggestedBranchLength) {  // brent suggests a new branch length


                                    double incumbentBranchLength = child.getParentDistance(parent);

                                    // mutate and see if it yields an improved network
                                    child.setParentDistance(parent, suggestedBranchLength);
                                    double lnProb = computeGTProb(speciesNetwork, geneTrees, counter);

                                    System.out.print("(" + parent.getName() + ", " + child.getName() + ")" + " to bl " + suggestedBranchLength + " yields " + lnProb + " vs " + lnGtProbOfSpeciesNetwork.getContents() );

                                    if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                                    {
                                        System.out.println(" (improved)");
                                        lnGtProbOfSpeciesNetwork.setContents(lnProb);
                                    }
                                    else  // didn't improve, roll back change
                                    {
                                        System.out.println("");
                                        child.setParentDistance(parent, incumbentBranchLength);
                                    }
                                    return lnProb;
                                }
                            };
                            BrentOptimizer optimizer = new BrentOptimizer(.000000000001,.0000000000000001); // very small numbers so we control when brent stops, not brent.

                            try
                            {
                                optimizer.optimize(_maxAssigmentAttemptsPerBranchParam, functionToOptimize, GoalType.MAXIMIZE, Double.MIN_VALUE, _maxBranchLength);
                            }
                            catch(TooManyEvaluationsException e) // _maxAssigmentAttemptsPerBranchParam exceeded
                            {
                            }
                        }
                    });
                }
            }

            // add hybrid probs to hybrid edges
            for(final NetNode<Double> child : speciesNetwork.bfs()) // find every hybrid node
            {
                if(child.isRoot()) // calling getParentNumber on root causes NPE. Bug workaround.
                    continue;

                if(child.getParentNumber() == 2)  // hybrid node
                {
                    Iterator<NetNode<Double>> hybridParents = child.getParents().iterator();
                    final NetNode<Double> hybridParent1 = hybridParents.next();
                    final NetNode<Double> hybridParent2 = hybridParents.next();

                    assigmentActions.add(new Proc()
                    {
                        public void execute()
                        {
                            UnivariateFunction functionToOptimize = new UnivariateFunction() {
                                public double value(double suggestedProb) {

                                    double incumbentHybridProbParent1 = child.getParentProbability(hybridParent1);

                                    // try new pair of hybrid probs
                                    child.setParentProbability(hybridParent1, suggestedProb);
                                    child.setParentProbability(hybridParent2, 1.0 - suggestedProb);


                                    double lnProb = computeGTProb(speciesNetwork, geneTrees, counter);


                                    System.out.print("(" + hybridParent1.getName() + ", " + child.getName() + ")" + " to hp " + suggestedProb + " yields " + lnProb + " vs " + lnGtProbOfSpeciesNetwork.getContents() );


                                    if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // change improved GTProb, keep it
                                    {
                                        System.out.println(" (improved)");
                                        lnGtProbOfSpeciesNetwork.setContents(lnProb);
                                    }
                                    else // change did not improve, roll back
                                    {
                                        System.out.println("");
                                        child.setParentProbability(hybridParent1, incumbentHybridProbParent1);
                                        child.setParentProbability(hybridParent2, 1.0 - incumbentHybridProbParent1);
                                    }
                                    return lnProb;
                                }
                            };
                            BrentOptimizer optimizer = new BrentOptimizer(.000000000001,.0000000000000001); // very small numbers so we control when brent stops, not brent.

                            try
                            {
                                optimizer.optimize(_maxAssigmentAttemptsPerBranchParam, functionToOptimize, GoalType.MAXIMIZE, 0, 1.0);
                            }
                            catch(TooManyEvaluationsException e)  // _maxAssigmentAttemptsPerBranchParam exceeded
                            {
                            }
                        }
                    });

                }
            }

            Collections.shuffle(assigmentActions); // randomize the order we will try to ajust network edge properties

            for(Proc assigment : assigmentActions)   // for each change attempt, perform attempt
            {
                assigment.execute();
            }


            if( ((double)lnGtProbOfSpeciesNetwork.getContents()) == lnGtProbLastRound)  // if no improvement was made wrt to last around, stop trying to find a better assignment
            {
                continueRounds = false;
            }
            else if (lnGtProbOfSpeciesNetwork.getContents() > lnGtProbLastRound) // improvement was made, ensure it is large enough wrt to improvement threshold to continue searching
            {
                double improvementPercentage = Math.pow(Math.E, (lnGtProbOfSpeciesNetwork.getContents() - lnGtProbLastRound)) - 1.0;  // how much did we improve over last round
                if(improvementPercentage < _improvementThreshold  )  // improved, but not enough to keep searching
                {
                    continueRounds = false;
                }
            }
            else
            {
                throw new IllegalStateException("Should never have decreased prob.");
            }
        }



        RnNewickPrinter<Double> rnNewickPrinter = new RnNewickPrinter<Double>();
        StringWriter sw = new StringWriter();
        rnNewickPrinter.print(speciesNetwork, sw);

        result.append("\nTotal log probability: " + lnGtProbOfSpeciesNetwork.getContents() + ": " + sw.toString());

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
