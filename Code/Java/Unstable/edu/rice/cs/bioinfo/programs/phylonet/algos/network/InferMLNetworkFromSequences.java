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
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkNeighbourhoodRandomWalkGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkWholeNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours.AllNeighboursHillClimberFirstBetter;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours.AllNeighboursHillClimberSteepestAscent;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import edu.rice.cs.bioinfo.programs.phmm.src.phylogeny.Felsenstein;
import edu.rice.cs.bioinfo.programs.phmm.src.substitutionModel.GTRSubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.integration.FullLikelihoodFromSequence;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BfsSearch;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;

import java.io.ByteArrayInputStream;
import java.io.StringReader;
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
public class InferMLNetworkFromSequences extends MDCOnNetworkYFFromRichNewickJung {
    private Network[] _optimalNetworks;
    private double[] _optimalScores;
    private int _maxRounds = 100;
    private int _maxTryPerBranch = 100;
    private double _improvementThreshold = 0.001;
    private double _maxBranchLength = 6;
    private double _maxHeight;
    private double _Brent1 = 0.01;
    private double _Brent2 = 0.001;
    private Long _maxExaminations = null;
    private long _maxFailure = 100;
    private int _diameterLimit;
    private Network _startNetwork;
    protected double[] _operationWeight = {0.15,0.15,0.2,0.5};
    protected int _numRuns = 1;
    protected int _numThreads = 1;
    private int _bins;
    private int _numSamples;
    private Long _seed = null;

    private int _currentSite;

    //private Map<SpeciesPair, Double> _pair2time;


    public InferMLNetworkFromSequences(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }


    public void setSearchParameter(int maxRounds, int maxTryPerBranch, double improvementThreshold, double maxBranchLength, double Brent1, double Brent2, Long maxExaminations, Long maxFailure, int diameterLimit, int parallel, Network startNetwork, double[] operationWeight, int numRuns, Long seed){
        _maxRounds = maxRounds;
        _maxTryPerBranch = maxTryPerBranch;
        _improvementThreshold = improvementThreshold;
        _maxBranchLength = maxBranchLength;
        _Brent1 = Brent1;
        _Brent2 = Brent2;
        _maxExaminations = maxExaminations;
        _diameterLimit = diameterLimit;
        _startNetwork = startNetwork;
        _maxFailure = maxFailure;
        _numThreads = parallel;
        _operationWeight = operationWeight;
        _numRuns = numRuns;
        _seed = seed;
    }

    public void setStartingNetwork(String network){
        _startNetwork = string2Network(network);
    }

    public void setIntegralParameter(int bins, int samples, double maxHeight){
        _bins = bins;
        _numSamples = samples;
        _maxHeight = maxHeight;
    }

    private List<Tuple<char[], Integer>> summarizeSquences(List<char[][]> sequences){
        List<Tuple<char[], Integer>> distinctSequences = new ArrayList<Tuple<char[], Integer>>();
        Map<String, Integer> alignment2count = new HashMap<String, Integer>();
        for (char[][] sa: sequences) {
            for(int i=0; i<sa[0].length; i++){
                String site = "";
                for(int j=0; j<sa.length; j++){
                    site += sa[j][i];
                }
                Integer count = alignment2count.get(site);
                if(count == null){
                    alignment2count.put(site, 1);
                }
                else{
                    alignment2count.put(site, count+1);
                }
            }
        }
        for(Map.Entry<String, Integer> entry: alignment2count.entrySet()){
            distinctSequences.add(new Tuple<char[], Integer>(entry.getKey().toCharArray(), entry.getValue()));
        }
        return distinctSequences;
    }


    public List<Tuple<Network,Double>> inferNetwork(String[] gtTaxa, List<char[][]> sequences, Map<String,List<String>> species2alleles, GTRSubstitutionModel gtrsm, double theta, int maxReticulations, int numSol){
        _optimalNetworks = new Network[numSol];
        _optimalScores = new double[numSol];
        Arrays.fill(_optimalScores, Double.NEGATIVE_INFINITY);

        List<Tuple<char[], Integer>> distinctSequences = summarizeSquences(sequences);
        List<Tree> allTrees = Trees.generateAllBinaryTrees(gtTaxa);
        String startingNetwork = network2String(_startNetwork);

        for(int i=0; i<_numRuns; i++) {
            DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = makeNetwork(startingNetwork);
            NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>> allNeighboursStrategy = new NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>>(_operationWeight, makeNode, makeEdge, _seed);
            AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>, Double> searcher = new AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);
            Func1<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, Double> scorer = getScoreFunction(gtTaxa, distinctSequences, allTrees, species2alleles, gtrsm, theta);
            Comparator<Double> comparator = getDoubleScoreComparator();
            HillClimbResult<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, Double> result = searcher.search(speciesNetwork, scorer, comparator, _maxExaminations, maxReticulations, _maxFailure, _diameterLimit); // search starts here
         }

        List<Tuple<Network, Double>> resultList = new ArrayList<Tuple<Network, Double>>();
        for(int i=0; i<numSol; i++){
            if(_optimalNetworks[i]!=null)
                resultList.add(new Tuple<Network, Double>(_optimalNetworks[i], _optimalScores[i]));
        }
        //System.out.println("\n #Networks " + result.ExaminationsCount);
        return resultList;
    }

    private Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }


    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final String[] gtTaxa, final List<Tuple<char[], Integer>> sequences, final List<Tree> allTrees, final Map<String, List<String>> species2alleles, final GTRSubstitutionModel gtrsm, final double theta){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {

                Network<Object> speciesNetwork = networkNew2Old(network);

                double score = findUltrametricOptimalBranchLength(speciesNetwork, gtTaxa, sequences, allTrees, species2alleles, gtrsm, theta);
                //System.out.println(network2String(speciesNetwork) + ": "+score);

                if(score > _optimalScores[_optimalNetworks.length-1]){
                    boolean exist = false;
                    for(int i=0; i<_optimalNetworks.length; i++){
                        if(_optimalNetworks[i]==null)break;

                        NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
                        if(metric.computeDistanceBetweenTwoNetworks(speciesNetwork, _optimalNetworks[i]) == 0){
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

                        System.out.println("Optimal Networks");
                        for(int i=0; i<_optimalNetworks.length; i++){
                            if(_optimalNetworks[i]!=null)
                                System.out.println(_optimalScores[i] + ": " + network2String(_optimalNetworks[i]));
                        }
                    }
                }
                System.out.println();
                //System.out.println(network2String(speciesNetwork) + ": "+score);
                //System.out.println();
                //System.out.println("End scoring ..." + (System.currentTimeMillis()-start)/1000.0);
                //System.exit(0);
                return score;
            }
        };
    }

    private void initializeNetwork(Network<Object> speciesNetwork, Map<NetNode, Double> node2height){
        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            double height = 0;
            if(!node.isLeaf()){
                for(NetNode child: node.getChildren()){
                    height = Math.max(height, node2height.get(child)+1);
                }
            }
            node2height.put(node, height);
            for(NetNode child: node.getChildren()){
                child.setParentDistance(node, height - node2height.get(child));
            }
        }
    }


    private double findUltrametricOptimalBranchLength(final Network<Object> speciesNetwork, final String[] gtTaxa, final List<Tuple<char[], Integer>> sequences, final List<Tree> allTrees, final Map<String, List<String>> species2alleles, final GTRSubstitutionModel gtrsm, final double theta){
        boolean continueRounds = true; // keep trying to improve network

        Map<NetNode, Double> node2height = new HashMap<NetNode, Double>();
        initializeNetwork(speciesNetwork, node2height);

        double initialProb;
        if(_numThreads == 1){
            initialProb = computeProbability(speciesNetwork, gtTaxa, sequences, allTrees, species2alleles, gtrsm, theta);
        }else{
            initialProb = computeProbabilityParallel(speciesNetwork, gtTaxa, sequences, allTrees,  species2alleles, gtrsm, theta);
        }

        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(initialProb);  // records the GTProb of the network at all times
        final Container<Map<NetNode, Double>> node2heightContainer = new Container<Map<NetNode, Double>>(node2height);

        int roundIndex = 0;
        for(; roundIndex <_maxRounds && continueRounds; roundIndex++)
        {
            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();

            for(final NetNode<Object> child : speciesNetwork.getNetworkNodes()) // find every hybrid node
            {


                Iterator<NetNode<Object>> hybridParents = child.getParents().iterator();
                final NetNode hybridParent1 = hybridParents.next();
                final NetNode hybridParent2 = hybridParents.next();

                UnivariateFunction functionToOptimize = new UnivariateFunction() {
                    public double value(double suggestedProb) {
                        double incumbentHybridProbParent1 = child.getParentProbability(hybridParent1);
                        child.setParentProbability(hybridParent1, suggestedProb);
                        child.setParentProbability(hybridParent2, 1.0 - suggestedProb);

                        double lnProb;
                        if(_numThreads == 1){
                            lnProb = computeProbability(speciesNetwork, gtTaxa, sequences, allTrees, species2alleles, gtrsm, theta);
                        }else{
                            lnProb = computeProbabilityParallel(speciesNetwork, gtTaxa, sequences, allTrees, species2alleles, gtrsm, theta);
                        }

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

            }

            for(final NetNode<Object> node : edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork))
            {
                if(node.isLeaf()){
                    continue;
                }

                final Container<Double> minHeight = new Container<Double>(0.0);
                final Container<Double> maxHeight = new Container<Double>(Double.MAX_VALUE);


                for(NetNode<Object> child: node.getChildren()){
                    double childHeight = node2heightContainer.getContents().get(child);
                    minHeight.setContents(Math.max(minHeight.getContents(), childHeight));
                }

                if(node.isRoot()){
                    maxHeight.setContents(minHeight.getContents() + _maxBranchLength);

                }
                else{
                    for(NetNode<Object> parent: node.getParents()){
                        double parentHeight = node2heightContainer.getContents().get(parent);
                        maxHeight.setContents(Math.min(maxHeight.getContents(), parentHeight));
                    }
                }


                //System.out.println("\nChanging node " + node.getName() + " from " + minHeight.getContents() + " to " + maxHeight.getContents());
                UnivariateFunction functionToOptimize = new UnivariateFunction() {

                    public double value(double suggestedHeight) {  // brent suggests a new branch length
                        double incumbentHeight = node2heightContainer.getContents().get(node);


                        for(NetNode<Object> child: node.getChildren()){
                            child.setParentDistance(node, suggestedHeight - node2heightContainer.getContents().get(child));
                        }

                        if(!node.isRoot()){
                            for(NetNode<Object> parent: node.getParents()){
                                node.setParentDistance(parent, node2heightContainer.getContents().get(parent) - suggestedHeight);
                            }
                        }

                        //System.out.println("trying " + network2String(speciesNetwork));
                        double lnProb;
                        if(_numThreads == 1){
                            lnProb = computeProbability(speciesNetwork, gtTaxa, sequences, allTrees, species2alleles, gtrsm, theta);
                        }else{
                            lnProb = computeProbabilityParallel(speciesNetwork, gtTaxa, sequences, allTrees, species2alleles, gtrsm, theta);
                        }
                        //System.out.print("suggest: "+ suggestedHeight + " " + lnProb + " ");
                        if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                        {
                            lnGtProbOfSpeciesNetwork.setContents(lnProb);
                            node2heightContainer.getContents().put(node, suggestedHeight);
                            //System.out.println( " better ");

                        }
                        else  // didn't improve, roll back change
                        {
                            for(NetNode<Object> child: node.getChildren()){
                                child.setParentDistance(node, incumbentHeight - node2heightContainer.getContents().get(child));
                            }
                            if(!node.isRoot()){
                                for(NetNode<Object> parent: node.getParents()){
                                    node.setParentDistance(parent, node2heightContainer.getContents().get(parent) - incumbentHeight);
                                }
                            }
                            //System.out.println( " worse ");
                        }
                        return lnProb;
                    }
                };
                BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                try
                {
                    optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, minHeight.getContents(), maxHeight.getContents());
                }
                catch(TooManyEvaluationsException e) // _maxAssigmentAttemptsPerBranchParam exceeded
                {
                }


                //System.out.println(network2String(speciesNetwork) + " " + lnGtProbOfSpeciesNetwork.getContents());
            }




            if( ((double)lnGtProbOfSpeciesNetwork.getContents()) == lnGtProbLastRound)  // if no improvement was made wrt to last around, stop trying to find a better assignment
            {
                continueRounds = false;
            }
            else if (lnGtProbOfSpeciesNetwork.getContents() > lnGtProbLastRound) // improvement was made, ensure it is large enough wrt to improvement threshold to continue searching
            {

                double improvementPercentage = Math.pow(Math.E, (lnGtProbOfSpeciesNetwork.getContents() - lnGtProbLastRound)) - 1.0;  // how much did we improve over last round
                //System.out.println(improvementPercentage + " vs. " + _improvementThreshold);
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

         return lnGtProbOfSpeciesNetwork.getContents();
    }


    private class MyThread extends Thread{
        Network _speciesNetwork;
        String[] _gtTaxa;
        List<Tuple<char[], Integer>> _sequences;
        List<Tree> _allTrees;
        Map<String, List<String>> _species2alleles;
        GTRSubstitutionModel _gtrsm;
        double _theta;
        double[] _probs;


        public MyThread(String speciesNetwork, String[] gtTaxa, List<Tuple<char[], Integer>> sequences, List<Tree> allTrees,  Map<String, List<String>> species2alleles, GTRSubstitutionModel gtrsm, double theta, double[] probs){
            _speciesNetwork = string2Network(speciesNetwork);
            _gtTaxa = gtTaxa;
            _sequences = sequences;
            _allTrees = allTrees;
            _species2alleles = species2alleles;
            _gtrsm = gtrsm;
            _theta = theta;
            _probs = probs;
        }


        public void run() {
            int locusID;
            while((locusID = getNextLocus())!=-1){
                _probs[locusID] = computeProbabilityForOneSite(_speciesNetwork, _gtTaxa, _sequences.get(locusID), _allTrees, _species2alleles, _gtrsm, _theta);
            }

        }

        private synchronized int getNextLocus(){
            _currentSite++;
            if(_currentSite >= _sequences.size()){
                return -1;
            }
            else{
                return _currentSite;
            }
        }
    }

    protected double computeProbabilityForOneSite(Network speciesNetwork, String[] gtTaxa, Tuple<char[],Integer> seq, List<Tree> allTrees,  Map<String, List<String>> species2alleles, GTRSubstitutionModel gtrsm, double theta){
        double locusProb = 0;
        Map<String, Character> sequenceMap = new HashMap<String, Character>();
        for(int i=0; i< gtTaxa.length; i++){
            sequenceMap.put(gtTaxa[i], seq.Item1[i]);
        }
        for(Tree gt: allTrees) {
            FullLikelihoodFromSequence integrator = new FullLikelihoodFromSequence(speciesNetwork, gt, species2alleles, sequenceMap, gtrsm, theta, 0, _maxHeight);
            //TODO
            //integrator.setPrintDetails(printDetails);
            //integrator.setRandomSeed(100);
            locusProb += integrator.computeLikelihoodWithIntegral(_numSamples, _bins);
        }
        return Math.log(locusProb) * seq.Item2;
    }


    protected double computeProbability(Network speciesNetwork, String[] gtTaxa, List<Tuple<char[], Integer>> sequences, List<Tree> allTrees,  Map<String, List<String>> species2alleles, GTRSubstitutionModel gtrsm, double theta) {
        double totalProb = 0;
        for(Tuple <char[], Integer> seq: sequences){
            double siteProb = computeProbabilityForOneSite(speciesNetwork, gtTaxa, seq, allTrees, species2alleles, gtrsm, theta);
            totalProb += siteProb;
        }
        return totalProb;
    }


    protected double computeProbabilityParallel(Network speciesNetwork, String[] gtTaxa, List<Tuple<char[], Integer>> sequences, List<Tree> allTrees,  Map<String, List<String>> species2alleles, GTRSubstitutionModel gtrsm, double theta) {
        double[] probArray = new double[sequences.size()];
        Thread[] myThreads = new Thread[_numThreads];

        String networkExp = network2String(speciesNetwork);
        _currentSite = -1;
        for(int i=0; i<_numThreads; i++){
            List<Tree> trCopies = new ArrayList<Tree>();
            for(Tree tr: allTrees){
                Tree copy = new STITree(tr);
                trCopies.add(copy);
            }
            myThreads[i] = new MyThread(networkExp, gtTaxa, sequences, trCopies, species2alleles, gtrsm, theta, probArray);
            myThreads[i].start();
        }

        for(int i=0; i<_numThreads; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }

        double totalProb = 0;
        for(double prob: probArray){
            totalProb += prob;
        }

        return totalProb;
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

    public static void main(String[] args){
        GTRSubstitutionModel gtrsm = new GTRSubstitutionModel();
        double[] rates = {1.05, 4.7346, 1.2441, 0.5014, 4.8320};
        double[] freqs = {0.2112, 0.2888, 0.2896, 0.2104};
        gtrsm.setSubstitutionRates(rates, freqs);

        String[] taxa = {"A","B","C"};
        char[][] sequenceAlignment1 = new char[3][5];
        String l1A = "AATCG";
        sequenceAlignment1[0] = l1A.toCharArray();
        String l1B = "ACGAC";
        sequenceAlignment1[1] = l1B.toCharArray();
        String l1C = "CATCG";
        sequenceAlignment1[2] = l1C.toCharArray();
        char[][] sequenceAlignment2 = new char[3][5];
        String l2A = "AATCC";
        sequenceAlignment2[0] = l2A.toCharArray();
        String l2B = "AATAG";
        sequenceAlignment2[1] = l2B.toCharArray();
        String l2C = "CCGCG";
        sequenceAlignment2[2] = l2C.toCharArray();
        List<char[][]> sequences = new ArrayList<char[][]>();
        sequences.add(sequenceAlignment2);
        sequences.add(sequenceAlignment1);
        InferMLNetworkFromSequences inference = new InferMLNetworkFromSequences();
        inference.setIntegralParameter(1, 5, 10);
        inference.setStartingNetwork("((B,C)i1,A)i2;");
        inference.inferNetwork(taxa, sequences, null, gtrsm, 1, 0, 1);
    }

  }
