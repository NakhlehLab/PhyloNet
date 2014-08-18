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

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Container;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;

import java.io.StringReader;
import java.io.StringWriter;
import java.util.*;


/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 3/13/12
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("SearchBranchLengthsMaxGTProb")
public class SearchBranchLengthsMaxGTProb extends CommandBaseFileOut{
    private HashMap<String,String> _taxonMap = null;
    private boolean  _printDetail = false;
    private NetworkNonEmpty _speciesNetwork;
    private List<NetworkNonEmpty> _geneTrees;
    private ParameterIdentList _geneTreeParam;
    private double _maxBranchLength;
    private int _maxAssigmentAttemptsPerBranchParam = -1;
    private int _assigmentRounds = Integer.MAX_VALUE;
    private double _improvementThreshold;

    public SearchBranchLengthsMaxGTProb(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                                        Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                                        Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
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

        ParamExtractor amParam = new ParamExtractor("am", this.params, this.errorDetected);
        if(amParam.ContainsSwitch)
        {
            _computeGTProbStrategy = this._computeGTProbStrategyBox;
        }
        else
        {
            _computeGTProbStrategy = this._computeGTProbStrategyCalGTProb;
        }


        ParamExtractor pParam = new ParamExtractor("p", this.params, this.errorDetected);
        if(pParam.ContainsSwitch)
        {
            _printDetail = true;
        }

        noError = noError && checkForUnknownSwitches("p", "a", "am");
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

                if(child.getParentCount() == 2)
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
        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(_computeGTProbStrategy.execute(speciesNetwork, geneTrees, counter));  // records the GTProb of the network at all times

        int assigmentRound = 0;
        for(; assigmentRound <_assigmentRounds && continueRounds; assigmentRound++)
        {

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
                    if(!parent.isRoot())  // jd tmp
                        continue;

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
                                    double lnProb = _computeGTProbStrategy.execute(speciesNetwork, geneTrees, counter);


                                    RnNewickPrinter<Double> rnNewickPrinter = new RnNewickPrinter<Double>();
                                    StringWriter sw = new StringWriter();
                                    rnNewickPrinter.print(speciesNetwork, sw);
                                    //   String inferredNetwork = sw.toString();
                                    //    System.out.println(inferredNetwork + "\t" + lnProb);

                                    if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                                    {
                                        lnGtProbOfSpeciesNetwork.setContents(lnProb); // System.out.println("(improved)");
                                    }
                                    else  // didn't improve, roll back change
                                    {
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

                            //   System.out.println("-----------------------------------------------------------------------");


                        }
                    });
                }
            }

            // add hybrid probs to hybrid edges
            for(final NetNode<Double> child : speciesNetwork.bfs()) // find every hybrid node
            {
                if(child.isRoot()) // calling getParentNumber on root causes NPE. Bug workaround.
                    continue;

                if(child.getParentCount() == 2)  // hybrid node
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


                                    double lnProb = _computeGTProbStrategy.execute(speciesNetwork, geneTrees, counter);


                                    if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // change improved GTProb, keep it
                                    {

                                        lnGtProbOfSpeciesNetwork.setContents(lnProb);
                                    }
                                    else // change did not improve, roll back
                                    {

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

            //     Collections.shuffle(assigmentActions); // randomize the order we will try to adjust network edge properties

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
        String inferredNetwork = sw.toString();

        this.richNewickGenerated(inferredNetwork);

        result.append("\nTotal log probability: " + lnGtProbOfSpeciesNetwork.getContents() + ": " + inferredNetwork);

        return result.toString();

    }

    private Func3<Network<Double>, List<Tree>, List<Integer>, Double> _computeGTProbStrategyCalGTProb = new Func3<Network<Double>, List<Tree>, List<Integer>, Double>() {
        public Double execute(Network<Double> speciesNetwork, List<Tree> geneTrees, List<Integer> counter) {

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
    };

    private Func3<Network<Double>, List<Tree>, List<Integer>, Double> _computeGTProbStrategyBox = new Func3<Network<Double>, List<Tree>, List<Integer>, Double>() {
        public Double execute(Network<Double> speciesNetwork, List<Tree> geneTrees, List<Integer> counter) {

            double t0True = 1.1;
            double t1True = 1.1;
            double gammaTrue = .5;
            double t0Star = speciesNetwork.findNode("D").getParentDistance(speciesNetwork.getRoot());
            double t1Star = speciesNetwork.findNode("F").getParentDistance(speciesNetwork.getRoot());
            double gammaStar = speciesNetwork.findNode("E").getParentProbability(speciesNetwork.findNode("D"));

            return new Box(t0True, t1True, gammaTrue, t0Star, t1Star, gammaStar, _geneTrees.size()).callnLikelihood();
        }
    };

    public class Box {
        public double t0,t1,gamma;
        public double P1,P2,P3;
        public double t0star, t1star, gammastar;
        public double P1star,P2star,P3star;
        public double lnLikelihood;
        public double MaxlnLikelihood, finalt0star, finalt1star, finalgammastar;
        public int n;

        public Box(double t0, double t1, double gamma, double t0star, double t1star, double gammastar, int n){
            this.t0 = t0;
            this.t1 = t1;
            this.gamma = gamma;
            this.t0star = t0star;
            this.t1star = t1star;
            this.gammastar = gammastar;
            this.n = n;
            P1 = (1-gamma)*(1-2.0/3.0*Math.exp(-t1))+gamma * Math.exp(-t0)/3.0;
            P2 = gamma*(1-2.0/3.0*Math.exp(-t0))+(1-gamma)*Math.exp(-t1)/3.0;
            P3 = (1-gamma)*Math.exp(-t1)/3.0+gamma*Math.exp(-t0)/3.0;
            callnLikelihood(); // initialize
        }

        // compute log likelihood function
        public double callnLikelihood(){
            P1star = (1-gammastar)*(1-2.0/3*Math.exp(-t1star))+gammastar*Math.exp(-t0star)/3.0;;
            P2star = gammastar*(1-2.0/3.0*Math.exp(-t0star))+(1-gammastar)*Math.exp(-t1star)/3.0;
            P3star = (1-gammastar)*Math.exp(-t1star)/3.0+gammastar*Math.exp(-t0star)/3.0;
            lnLikelihood = n*(P1*Math.log(P1star)+P2*Math.log(P2star)+P3*Math.log(P3star));
            return lnLikelihood;
        }
        /*
          // A simple method to find the max log likelihood function for one branch
      public void findRestrictedMax(int whichone){
          double first = Double.MIN_VALUE;
          double second = Double.MIN_VALUE;
          double third = Double.MIN_VALUE;
          if (whichone == 0) { // t0star
              for (int i=1; i<=1000; i++) {
                  t0star = 0.01*i;
                  callnLikelihood();
                  third = lnLikelihood;

                  if (second>first && second> third) { // we found the max
                      lnLikelihood = second;
                      t0star = 0.01*(i-1);
                      break;
                  }
                  else {
                      first = second;	second = third;
                  }
              }
          } // t0 case

          else if (whichone == 1) { //t1star
              for (int i=1;i<=1000;i++) {
                  t1star = 0.01*i;
                  callnLikelihood();
                  third = lnLikelihood;
                  if (second>first && second> third) { // we found the max
                      lnLikelihood = second;
                      t1star = 0.01*(i-1);
                      break;
                  }
                  else {
                      first = second;	second = third;
                  }
              }
          } // t1star case

          else if (whichone == 2) {  // gammastar
              for (int i=1;i<=1000;i++) {
                  gammastar = 0.001*i;
                  callnLikelihood();
                  third = lnLikelihood;
                  if (second>first && second> third) { // we found the max
                      lnLikelihood = second;
                      gammastar = 0.001*(i-1);
                      break;
                  }
                  else {
                      first = second;	second = third;
                  }
              }
          }  // gammastar case
      }

      public void OptEdgeLen() {
          double lnLikelihoodCurBest = lnLikelihood; // get the initial lnLikelihood value
          double improveratio;
          while (true) {
              findRestrictedMax(0);
              assert lnLikelihood >= lnLikelihoodCurBest;
              improveratio = Math.exp(lnLikelihood - lnLikelihoodCurBest);
              if (improveratio < 1.001) {
                  break;
              }
              else
                  lnLikelihoodCurBest = lnLikelihood;

              findRestrictedMax(1);
              assert lnLikelihood >= lnLikelihoodCurBest;
              improveratio = Math.exp(lnLikelihood - lnLikelihoodCurBest);
              if (improveratio < 1.001) {
                  break;
              }
              else
                  lnLikelihoodCurBest = lnLikelihood;

          }
          MaxlnLikelihood = lnLikelihood;
          finalt0star = t0star;
          finalt1star = t1star;
          finalgammastar = gammastar;
      }        */
    }

    private Func3<Network<Double>, List<Tree>, List<Integer>, Double> _computeGTProbStrategy = _computeGTProbStrategyCalGTProb;


}
