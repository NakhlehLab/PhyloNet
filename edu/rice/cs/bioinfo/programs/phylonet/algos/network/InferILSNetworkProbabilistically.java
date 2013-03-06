/*
 * Copyright (c) 2013 Rice University.
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

package edu.rice.cs.bioinfo.programs.phylonet.algos.network;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.RichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.phylogenetics.FindRoot;
import edu.rice.cs.bioinfo.library.phylogenetics.GetDirectSuccessors;
import edu.rice.cs.bioinfo.library.phylogenetics.GetInDegree;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkWholeNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours.AllNeighboursHillClimberSteepestAscent;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;

import java.io.ByteArrayInputStream;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferILSNetworkProbabilistically extends MDCOnNetworkYFFromRichNewickJung {
    private Network[] _optimalNetworks;
    private double[] _optimalScores;
    int _maxRounds;
    int _maxTryPerBranch;
    double _improvementThreshold;
    double _maxBranchLength;
    double _Brent1;
    double _Brent2;


    public InferILSNetworkProbabilistically(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }

    public void setBrentParameter(int maxRounds, int maxTryPerBranch, double improvementThreshold, double maxBranchLength, double Brent1, double Brent2){
        _maxRounds = maxRounds;
        _maxTryPerBranch = maxTryPerBranch;
        _improvementThreshold = improvementThreshold;
        _maxBranchLength = maxBranchLength;
        _Brent1 = Brent1;
        _Brent2 = Brent2;
    }

    public List<Tuple<Network,Double>> inferNetwork(List<Tree> gts, Map<String,List<String>> species2alleles, Long maxExaminations, Long maxReticulations, int diameterLimit, Network startNetwork, int numSol){
        _optimalNetworks = new Network[numSol];
        _optimalScores = new double[numSol];
        Arrays.fill(_optimalScores, Double.NEGATIVE_INFINITY);

        List<Tree> distinctTrees = new ArrayList<Tree>();
        List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList = new ArrayList<Tuple3<Tree, Double, List<Integer>>>();
        summarizeGeneTrees(gts, distinctTrees, nbTreeAndCountAndBinaryIDList);

        DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = getStartNetwork(gts, species2alleles,startNetwork);
        NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
        AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);

        Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> scorer = getScoreFunction(distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);
        Comparator<Double> comparator = getDoubleScoreComparator();
        HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, maxExaminations, maxReticulations, diameterLimit); // search starts here
        List<Tuple<Network, Double>> resultList = new ArrayList<Tuple<Network, Double>>();
        for(int i=0; i<numSol; i++){
            resultList.add(new Tuple<Network, Double>(_optimalNetworks[i], _optimalScores[i]));
        }
        //System.out.println("\n #Networks " + result.ExaminationsCount);
        return resultList;
    }


    private DirectedGraphToGraphAdapter<String,PhyloEdge<String>> getStartNetwork(List<Tree> gts, Map<String,List<String>> species2alleles, Network<Object> startingNetwork){
        if(startingNetwork == null){
            Map<String,String> allele2species = null;
            if(species2alleles!=null){
                allele2species = new HashMap<String, String>();
                for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                    String species = entry.getKey();
                    for(String allele: entry.getValue()){
                        allele2species.put(allele,species);
                    }
                }
            }
            MDCInference_DP mdc = new MDCInference_DP();
            Solution sol;
            if(allele2species==null){
                sol = mdc.inferSpeciesTree(gts, false, 1, false, 100, true, -1).get(0);
            }
            else{
                sol = mdc.inferSpeciesTree(gts, allele2species, false, 1, false, 100, true, -1).get(0);
            }

            Tree startingTree= Trees.generateRandomBinaryResolution(sol._st);
            startingNetwork = string2Network(startingTree.toString());
        }



        int index = 1;
        for(NetNode<Object> node: startingNetwork.dfs()){
            if(node.getName()==null || node.getName().equals("")){
                String name;
                do{
                    name = "i" + (index++);
                }while(startingNetwork.findNode(name)!=null);
                node.setName(name);
            }
        }

        String newNetwork = network2String(startingNetwork);

        return makeNetwork(newNetwork);

    }

    private void summarizeGeneTrees(List<Tree> originalGTs, List<Tree> distinctGTs, List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        for(Tree tr: originalGTs){
            double weight = ((STINode<Double>)tr.getRoot()).getData();
            int index = 0;
            Tuple3<Tree, Double, List<Integer>> newTuple = null;
            for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
                if(Trees.haveSameRootedTopology(tr, triple.Item1)){
                    newTuple = new Tuple3<Tree, Double, List<Integer>>(tr, triple.Item2 + weight, triple.Item3);
                    break;
                }
                index++;
            }
            if(newTuple!=null){
                nbTreeAndCountAndBinaryIDList.set(index, newTuple);
            }
            else{
                List<Integer> binaryIDs = new ArrayList<Integer>();
                for(Tree btr: Trees.getAllBinaryResolution(tr)){
                    index = 0;
                    boolean exist = false;
                    for(Tree exTr: distinctGTs){
                        if(Trees.haveSameRootedTopology(btr, exTr)){
                            binaryIDs.add(index);
                            exist = true;
                            break;
                        }
                        index++;
                    }
                    if(!exist){
                        distinctGTs.add(btr);
                        binaryIDs.add(index);
                    }
                }
                newTuple = new Tuple3<Tree, Double, List<Integer>>(tr, weight, binaryIDs);
                nbTreeAndCountAndBinaryIDList.add(newTuple);
            }
        }
    }

    private Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }


    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                Network<Object> speciesNetwork = networkNew2Old(network);
                double score = findOptimalBranchLength(speciesNetwork, distinctTrees, species2alleles,nbTreeAndCountAndBinaryIDList);

                if(score > _optimalScores[_optimalNetworks.length-1]){
                    boolean exist = false;
                    for(int i=0; i<_optimalNetworks.length; i++){
                        if(_optimalNetworks[i]==null)break;

                        if(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeClusterDistance(speciesNetwork, _optimalNetworks[i])[2]<0.000001 &&
                                edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeTripartitionDistance(speciesNetwork, _optimalNetworks[i])[2]<0.000001){
                            exist = true;
                            break;
                        }
                    }
                    if(!exist){
                        int index = -1;
                        for(int i=0; i<_optimalScores.length; i++){
                            if(score > _optimalScores[i]){
                                index = i;
                                break;
                            }
                        }
                        for(int i=_optimalScores.length-1; i>index; i--){
                            _optimalNetworks[i] = _optimalNetworks[i-1];
                            _optimalScores[i] = _optimalScores[i-1];
                        }
                        _optimalScores[index] = score;
                        _optimalNetworks[index] = string2Network(network2String(speciesNetwork));
                        //System.out.println(network2String(speciesNetwork) + ": "+score);
                    }
                }
                return score;
            }
        };
    }


    private double findOptimalBranchLength(final Network<Object> speciesNetwork, final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        boolean continueRounds = true; // keep trying to improve network

        //boolean isNetwork = false;
        for(NetNode<Object> node: speciesNetwork.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent,1.0);
                if(node.isNetworkNode()){
                    node.setParentProbability(parent, 0.5);
                    //isNetwork = true;
                }
            }
        }

        //final boolean hasNetworkNode = isNetwork;
        //long start = System.currentTimeMillis();
        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList));  // records the GTProb of the network at all times
        //System.out.println("first one : " + (System.currentTimeMillis()-start)/1000.0);
        //final Container<Double> timeTotal = new Container<Double>(System.currentTimeMillis()/1000.0);  // records the GTProb of the network at all times
        //final Container<Integer> timeCount = new Container<Integer>(0);  // records the GTProb of the network at all times


        int roundIndex = 0;
        for(; roundIndex <_maxRounds && continueRounds; roundIndex++)
        {
            /*
            * Prepare a random ordering of network edge examinations each of which attempts to change a branch length or hybrid prob to improve the GTProb score.
            */

            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();
            List<Proc> assigmentActions = new ArrayList<Proc>(); // store adjustment commands here.  Will execute them one by one later.

            // add branch length adjustments to the list

            for(final NetNode<Object> parent : speciesNetwork.bfs())
            {

                for(final NetNode<Object> child : parent.getChildren())
                {
                    if(child.isLeaf()) // leaf edge, skip
                        continue;
                    //System.out.println(parent.getName() + " -> " + child.getName());
                    // records the GTProb of the network at all times


                    assigmentActions.add(new Proc()
                    {
                        public void execute()
                        {
                            //timeTotal.setContents(System.currentTimeMillis()/1000.0);
                            //final Container<Integer> callCount = new Container<Integer>(0);
                            UnivariateFunction functionToOptimize = new UnivariateFunction() {
                                public double value(double suggestedBranchLength) {  // brent suggests a new branch length
                                    //callCount.setContents(callCount.getContents()+1);
                                    double incumbentBranchLength = child.getParentDistance(parent);


                                    // mutate and see if it yields an improved network


                                    child.setParentDistance(parent, suggestedBranchLength);

                                    double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, parent);
/*
                                            timeCount.setContents(timeCount.getContents()+1);
                                            if(timeCount.getContents()%1000==0){
                                                //timeTotal.setContents(timeTotal.getContents()+(System.currentTimeMillis()/1000.0-timeTotal.getContents()));
                                                System.out.println(timeCount.getContents() + ": " +(System.currentTimeMillis()/1000.0-timeTotal.getContents())/timeCount.getContents());
                                            }
*/
                                    /*
                                    RnNewickPrinter<Object> rnNewickPrinter = new RnNewickPrinter<Object>();
                                    StringWriter sw = new StringWriter();
                                    rnNewickPrinter.print(speciesNetwork, sw);
                                       String inferredNetwork = sw.toString();
                                        System.out.println(inferredNetwork + "\t" + lnProb);
                                    */
                                    if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                                    {
                                        lnGtProbOfSpeciesNetwork.setContents(lnProb); // System.out.println("(improved)");
                                        //System.out.println(network2String(speciesNetwork) + ":" + lnProb);
                                    }
                                    else  // didn't improve, roll back change
                                    {
                                        child.setParentDistance(parent, incumbentBranchLength);
                                    }
                                    return lnProb;
                                }
                            };
                            BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                            try
                            {
                                optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, Double.MIN_VALUE, _maxBranchLength);
                            }
                            catch(TooManyEvaluationsException e) // _maxAssigmentAttemptsPerBranchParam exceeded
                            {
                            }

                            //System.out.println(computeProbability(speciesNetwork, geneTrees, counter, child, parent) + " vs. " + computeProbability(speciesNetwork, geneTrees, counter));

                            //   System.out.println("-----------------------------------------------------------------------");
                            //System.out.println(callCount.getContents());
                            computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, parent);
                        }
                    });
                }
            }


            // add hybrid probs to hybrid edges
            for(final NetNode<Object> child : speciesNetwork.bfs()) // find every hybrid node
            {
                if(child.isRoot()) // calling getParentNumber on root causes NPE. Bug workaround.
                    continue;

                if(child.getParentNumber() == 2)  // hybrid node
                {
                    Iterator<NetNode<Object>> hybridParents = child.getParents().iterator();
                    final NetNode hybridParent1 = hybridParents.next();
                    final NetNode hybridParent2 = hybridParents.next();
                    //child.setParentProbability(hybridParent1,0.5);
                    //child.setParentProbability(hybridParent2,0.5);
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

                                    //long start = System.currentTimeMillis();
                                    double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, null);
                                    //System.out.print((System.currentTimeMillis()-start)/1000.0 + " ");
                                    //double lnProb = computeProbability(speciesNetwork, geneTrees, counter);
                                    /*
                                    timeCount.setContents(timeCount.getContents()+1);
                                    if(timeCount.getContents()%1000==0){
                                        //timeTotal.setContents(timeTotal.getContents()+(System.currentTimeMillis()/1000.0-timeTotal.getContents()));
                                        System.out.println(timeCount.getContents() + ": " +(System.currentTimeMillis()/1000.0-timeTotal.getContents())/timeCount.getContents());
                                    }
                                    */
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
                            BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                            try
                            {
                                optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, 0, 1.0);
                            }
                            catch(TooManyEvaluationsException e)  // _maxAssigmentAttemptsPerBranchParam exceeded
                            {
                            }
                            computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, null);
                            //System.out.println(computeProbability(speciesNetwork, geneTrees, counter, child, null) + " vs. " + computeProbability(speciesNetwork, geneTrees, counter));

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

        //System.out.println(standardComputeProbability(speciesNetwork, geneTrees, counter) + " vs. " + lnGtProbOfSpeciesNetwork.getContents());
        return lnGtProbOfSpeciesNetwork.getContents();
    }

    private double computeProbability(Network speciesNetwork, List<Tree> geneTrees, Map<String, List<String>> species2alleles, List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList) {
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        List<Double> probList = gtp.calculateGTDistribution(speciesNetwork, geneTrees, species2alleles);
        double total = 0;
        for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double maxProb = 0;
            for(int id: triple.Item3){
                maxProb = Math.max(maxProb, probList.get(id));
            }
            total += Math.log(maxProb) * triple.Item2;
        }
        return total;
    }




    public double computeProbability(Network speciesNetwork, List<Tree> geneTrees, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList, NetNode child, NetNode parent) {
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        List<Double> probList = gtp.calculateGTDistribution(speciesNetwork, geneTrees, child, parent);
        double total = 0;
        for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double maxProb = 0;
            for(int id: triple.Item3){
                maxProb = Math.max(maxProb, probList.get(id));
            }
            total += Math.log(maxProb) * triple.Item2;
        }
        return total;
    }

    private String network2String(Network speciesNetwork){
        RnNewickPrinter<Double> rnNewickPrinter = new RnNewickPrinter<Double>();
        StringWriter sw = new StringWriter();
        rnNewickPrinter.print(speciesNetwork, sw);
        return sw.toString();
    }

    private String network2String(final DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork){
        Func1<String,String> _getNetworkNodeLabel  = new Func1<String,String>()
        {
            public String execute(String node) {
                return node;
            }
        };

        Func1<String, Iterable<String>> _getDestinationNodes = new Func1<String, Iterable<String>>() {
            public Iterable<String> execute(String node) {
                return new GetDirectSuccessors<String,PhyloEdge<String>>().execute(speciesNetwork, node);
            }
        };


        Func2<String, String, String> _getNetworkDistanceForPrint = new Func2<String, String, String>() {
            public String execute(String parent, String child) {
                PhyloEdge<String> edge = speciesNetwork.getEdge(parent, child);
                if(edge.getBranchLength()==null){
                    return null;
                }
                return edge.getBranchLength()+"";
            }
        };

        Func2<String, String, String> _getProbabilityForPrint = new Func2<String, String, String>() {
            public String execute(String parent, String child) {
                PhyloEdge<String> edge = speciesNetwork.getEdge(parent, child);
                if(edge.getProbability()==null){
                    return null;
                }
                return edge.getProbability()+"";
            }
        };

        Func2<String, String, String> _getSupportForPrint = new Func2<String, String, String>() {
            public String execute(String parent, String child) {
                PhyloEdge<String> edge = speciesNetwork.getEdge(parent, child);
                if(edge.getSupport()==null){
                    return null;
                }
                return edge.getSupport()+"";
            }
        };

        Func1<String, HybridNodeType> _getHybridTypeForPrint = new Func1<String, HybridNodeType>()
        {
            public HybridNodeType execute(String node)
            {
                int inDegree = new GetInDegree<String,PhyloEdge<String>>().execute(speciesNetwork, node);
                return inDegree == 2 ? HybridNodeType.Hybridization : null;
            }
        };

        Func1<String,String> _getHybridNodeIndexForPrint = new Func1<String, String>() {
            List<String> hybridNodes = new ArrayList<String>();

            public String execute(String node) {
                int inDegree = new GetInDegree<String,PhyloEdge<String>>().execute(speciesNetwork, node);
                if(inDegree == 2){
                    int index = hybridNodes.indexOf(node) + 1;
                    if(index == 0){
                        hybridNodes.add(node);
                        return hybridNodes.size()+"";
                    }
                    else{
                        return index + "";
                    }
                }
                else{
                    return null;
                }
            }
        };

        try{
            StringWriter sw = new StringWriter();
            //   new RichNewickPrinterCompact<String>().print(true, "R", _getNetworkNodeLabel, _getDestinationNodes, _getNetworkDistanceForPrint, _getSupportForPrint, _getProbabilityForPrint, _getHybridNodeIndexForPrint, _getHybridTypeForPrint, sw);
            RichNewickPrinterCompact<String> printer = new RichNewickPrinterCompact<String>();
            printer.setGetBranchLength(_getNetworkDistanceForPrint);
            printer.setGetProbability(_getProbabilityForPrint);
            printer.setGetSupport(_getSupportForPrint);

            printer.print(true, new FindRoot<String>().execute(speciesNetwork), _getNetworkNodeLabel, _getDestinationNodes, _getHybridNodeIndexForPrint, _getHybridTypeForPrint, sw);
            sw.flush();
            sw.close();
            return sw.toString();
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return null;
    }




    private Network networkNew2Old(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork){
        Network<Object> bniNetwork = null;
        try{

            String networkString = network2String(speciesNetwork);
            bniNetwork = string2Network(networkString);
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return bniNetwork;
    }

    private Network string2Network(String networkString){
        try{
            RichNewickReaderAST reader = new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
            reader.setHybridSumTolerance(BigDecimal.valueOf(0.00001));
            NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
            RichNewickReadResult<Networks> readResult = reader.read(new ByteArrayInputStream(networkString.getBytes()));
            return transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return null;
    }

}
