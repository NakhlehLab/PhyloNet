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
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours.AllNeighboursHillClimberFirstBetter;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.integration.FullLikelihoodFromSequence;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
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

import java.io.*;
import java.math.BigDecimal;
import java.util.*;
import java.util.concurrent.*;

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
    private double _improvementThreshold = 0.01;
    private double _maxBranchLength = 1;
    private double _Brent1 = 0.05;
    private double _Brent2 = 0.005;
    private Long _maxExaminations = null;
    private long _maxFailure = 100;
    private int _diameterLimit;
    private Network _startNetwork;
    protected double[] _operationWeight = {0.15,0.15,0.2,0.5};
    protected int _numRuns = 1;
    protected int _numThreads = 1;
    protected Long _seed = null;
    private File resultFile = null;
    private List<Tuple<Network,Double>> _networksTried = new ArrayList<Tuple<Network, Double>>();


    public void setResultFile(File file){
        resultFile = file;
    }


    public void setParallel(int numThreads){
        _numThreads = numThreads;
    }


    public void setPreviousResult(Network[] optimalNetworks, double[] optimalScores, int numRuns){
        _optimalNetworks = optimalNetworks;
        _optimalScores = optimalScores;
        _numRuns = numRuns;
    }


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

    public void setStartingNetwork(Network network){
        _startNetwork = network;
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




    public List<Tuple<Network,Double>> inferNetwork(String[] gtTaxa, List<char[][]> sequences, Map<String,String> allele2species, RateModel gtrsm, double theta, int maxReticulations, int numSol){
        if(_optimalNetworks==null) {
            _optimalNetworks = new Network[numSol];
            _optimalScores = new double[numSol];
            Arrays.fill(_optimalScores, Double.NEGATIVE_INFINITY);
        }

        List<Tuple<char[], Integer>> distinctSequences = summarizeSquences(sequences);
        String startingNetwork = _startNetwork.toString();

        for(int i=0; i<_numRuns; i++) {
            System.out.println("Run #" + i);
            DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = makeNetwork(startingNetwork);
            NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>> allNeighboursStrategy = new NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>>(_operationWeight, makeNode, makeEdge, _seed);
            AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>, Double> searcher = new AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);
            Func1<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, Double> scorer = getScoreFunction(gtTaxa, allele2species, distinctSequences, gtrsm, theta);
            Comparator<Double> comparator = getDoubleScoreComparator();
            HillClimbResult<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, Double> result = searcher.search(speciesNetwork, scorer, comparator, _maxExaminations, maxReticulations, _maxFailure, _diameterLimit); // search starts here
            if(resultFile!=null) {
                try {
                    BufferedWriter bw = new BufferedWriter(new FileWriter(resultFile, true));
                    bw.append("Run #" + i + "\n");
                    for (int j = 0; j < _optimalNetworks.length; j++) {
                        bw.append(_optimalScores[j] + ": " + _optimalNetworks[j].toString() + "\n");
                    }
                    bw.newLine();
                    bw.close();
                } catch (Exception e) {
                    System.err.println(e.getMessage());
                    e.getStackTrace();
                }
            }
        }

        List<Tuple<Network, Double>> resultList = new ArrayList<Tuple<Network, Double>>();
        for(int i=0; i<numSol; i++){
            if(_optimalNetworks[i]!=null)
                resultList.add(new Tuple<Network, Double>(_optimalNetworks[i], _optimalScores[i]));
        }
        //System.out.println("\n #Networks " + result.ExaminationsCount);
        return resultList;
    }

    private Double hasAlreadyComputed(Network network){
        for(Tuple<Network,Double> networkTried: _networksTried){
            if(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.hasTheSameTopology(networkTried.Item1, network)){
                return networkTried.Item2;
            }
        }
        return null;
    }

    private Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }


    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final String[] gtTaxa, final Map<String, String> allele2species, final List<Tuple<char[], Integer>> sequences, final RateModel rModel, final double theta){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {

                Network<Object> speciesNetwork = networkNew2Old(network);
                Double score = hasAlreadyComputed(speciesNetwork);
                if(score == null) {
                    long start = System.currentTimeMillis();
                    score = findNonUltrametricOptimalBranchLength(gtTaxa, speciesNetwork, allele2species, sequences, rModel, theta);
                    System.out.println(System.currentTimeMillis() - start);
                    _networksTried.add(new Tuple<Network, Double>(cloneNetwork(speciesNetwork), score));
                }
                System.out.println(speciesNetwork.toString() + ": "+score);

                if(score > _optimalScores[_optimalNetworks.length-1]){
                    boolean exist = false;
                    for(int i=0; i<_optimalNetworks.length; i++){
                        if(_optimalNetworks[i]==null)break;

                        NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
                        if(metric.computeDistanceBetweenTwoNetworks(speciesNetwork, _optimalNetworks[i]) == 0){
                            exist = true;
                            if(i != 0 && score > _optimalScores[i-1]) {
                                int index = i;
                                while(index >= 1 && score > _optimalScores[index-1]){
                                    _optimalNetworks[index] = _optimalNetworks[index-1];
                                    _optimalScores[index] = _optimalScores[index-1];
                                    index--;
                                }
                                _optimalNetworks[index] = cloneNetwork(speciesNetwork);
                                _optimalScores[index] = score;
                            }
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
                        _optimalNetworks[index] = cloneNetwork(speciesNetwork);


                         System.out.println("Optimal Networks");
                         for(int i=0; i<_optimalNetworks.length; i++){
                             if(_optimalNetworks[i]!=null)
                                 System.out.println(_optimalScores[i] + ": " + _optimalNetworks[i].toString());
                         }

                    }
                }
                //System.out.println();
                //System.out.println(network2String(speciesNetwork) + ": "+score);
                //System.out.println();
                //System.out.println("End scoring ..." + (System.currentTimeMillis()-start)/1000.0);
                //System.exit(0);
                return score;
            }
        };
    }

    private Network cloneNetwork(Network network){
        return edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork(network.toString());

    }



    protected double findNonUltrametricOptimalBranchLength(final String[] gtTaxa, final Network<Object> speciesNetwork, final Map<String, String> allele2species, final List<Tuple<char[], Integer>> sequences, final RateModel rModel, final double theta){
        boolean continueRounds = true; // keep trying to improve network


        for(NetNode<Object> node: speciesNetwork.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent,theta);
                if(node.isNetworkNode()){
                    //node.setParentDistance(parent,0.0);
                    node.setParentProbability(parent, 0.5);
                }
            }
        }


        //long start = System.currentTimeMillis();
        double initalProb = computeProbability(gtTaxa, speciesNetwork, allele2species, sequences, rModel, theta);
        //System.out.println(initalProb);
        //System.out.print("\n"+(System.currentTimeMillis()-start));
        //computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);

        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(initalProb);  // records the GTProb of the network at all times
        //final Container<Integer> callCount = new Container<Integer>(0);


        int roundIndex = 0;
        for(; roundIndex <_maxRounds && continueRounds; roundIndex++)
        {
             /*
             * Prepare a random ordering of network edge examinations each of which attempts to change a branch length or hybrid prob to improve the GTProb score.
             */
            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();
            List<Proc> assigmentActions = new ArrayList<Proc>(); // store adjustment commands here.  Will execute them one by one later.


            for(final NetNode<Object> parent : edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork))
            {

                for(final NetNode<Object> child : parent.getChildren())
                {

                    assigmentActions.add(new Proc()
                    {
                        public void execute()
                        {

                            UnivariateFunction functionToOptimize = new UnivariateFunction() {
                                public double value(double suggestedBranchLength) {
                                    //System.out.print(" l/"+child.getName()+" ");
                                    //callCount.setContents(callCount.getContents()+1);
                                    double incumbentBranchLength = child.getParentDistance(parent);

                                    child.setParentDistance(parent, suggestedBranchLength);

                                    double lnProb = computeProbability(gtTaxa, speciesNetwork, allele2species, sequences, rModel, theta);
                                    //System.out.println("Changing branch ("+parent.getName()+","+child.getName()+") to " + suggestedBranchLength);
                                    //System.out.println(network2String(speciesNetwork)+": " + lnProb);

                                    if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                                    {
                                        lnGtProbOfSpeciesNetwork.setContents(lnProb);

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

                        }
                    });
                }
            }


            for(final NetNode<Object> child : speciesNetwork.getNetworkNodes()) // find every hybrid node
            {

                Iterator<NetNode<Object>> hybridParents = child.getParents().iterator();
                final NetNode hybridParent1 = hybridParents.next();
                final NetNode hybridParent2 = hybridParents.next();

                assigmentActions.add(new Proc()
                {
                    public void execute()
                    {
                        UnivariateFunction functionToOptimize = new UnivariateFunction() {
                            public double value(double suggestedProb) {
                                //callCount.setContents(callCount.getContents()+1);
                                //System.out.print(" p/"+child.getName()+" ");
                                double incumbentHybridProbParent1 = child.getParentProbability(hybridParent1);

                                child.setParentProbability(hybridParent1, suggestedProb);
                                child.setParentProbability(hybridParent2, 1.0 - suggestedProb);

                                double lnProb = computeProbability(gtTaxa, speciesNetwork, allele2species, sequences, rModel, theta);
                                //System.out.println("Changing node probability to "+ suggestedProb);
                                //System.out.println(network2String(speciesNetwork)+": " + lnProb);
                                //System.out.println(Math.abs(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList) - lnProb));
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
                        //System.out.println(network2String(speciesNetwork) + " : " + lnGtProbOfSpeciesNetwork.getContents());
                    }
                });


            }


            // add hybrid probs to hybrid edges
            if(_seed==null) {
                //Collections.shuffle(assigmentActions);
            }
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
        //System.out.println(callCount.getContents());
        //System.out.println(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList) + " vs. " + lnGtProbOfSpeciesNetwork.getContents());
        return lnGtProbOfSpeciesNetwork.getContents();
    }



    protected double computeProbability(String[] gtTaxa, Network speciesNetwork, Map<String, String> allele2species, List<Tuple<char[], Integer>> sequences, RateModel rModel, double theta) {
        //long startTime = System.currentTimeMillis();
        double totalProb = 0;
        ExecutorService executor = Executors.newFixedThreadPool(_numThreads);
        List<Future<Double>> futureList = new ArrayList<Future<Double>>();
        int chunkSize = sequences.size()/_numThreads;
        int extraNumber = sequences.size()%_numThreads;
        int extra = 0;
        int startingPoint = 0;
        for(int i=0; i<_numThreads; i++) {
            int size = extra++<extraNumber? chunkSize+1: chunkSize;
            MyThread myThread = new MyThread(speciesNetwork.toString(),gtTaxa,allele2species,sequences.subList(startingPoint, startingPoint+size),rModel, theta);
            startingPoint += size;
            Future<Double> future = executor.submit(myThread);
            futureList.add(future);
        }
        for(Future<Double> future : futureList){
            try {
                totalProb += future.get();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        executor.shutdown();
        //System.out.println((System.currentTimeMillis()-startTime)/1000.0);
        System.out.println(speciesNetwork.toString() + ": " + totalProb);
        return totalProb;
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
            bniNetwork = edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork(networkString);
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return bniNetwork;
    }


    private class MyThread implements Callable<Double>{
        String _network;
        String[] _gtTaxa;
        Map<String, String> _allele2species;
        List<Tuple<char[], Integer>> _sequences;
        RateModel _rModel;
        double _theta;


        public MyThread(String network, String[] gtTaxa, Map<String,String> allele2species, List<Tuple<char[], Integer>> sequences, RateModel rModel, double theta){
            _network = network;
            _gtTaxa = gtTaxa;
            _allele2species = allele2species;
            _sequences = sequences;
            _rModel = rModel;
            _theta = theta;
        }


        public Double call(){
            double totalProb = 0;
            for(Tuple<char[], Integer> seqTuple: _sequences){
                Map<String,Character> snapMap = new HashMap<String, Character>();
                for(int i=0; i<_gtTaxa.length; i++){
                    snapMap.put(_gtTaxa[i], seqTuple.Item1[i]);
                }
                OneNucleotideObservation converter = new OneNucleotideObservation(snapMap);
                /*
                Tree st = null;
                try{
                    NewickReader nr = new NewickReader(new StringReader(_network));
                    st = nr.readTree();
                }catch (Exception e){
                    System.err.println(e.getMessage());
                    e.getStackTrace();
                }
                SNAPPAlgorithm run = new SNAPPAlgorithm(st, _rModel, _theta);
                */

                SNAPPAlgorithm run = new SNAPPAlgorithm(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork(_network), _rModel, _theta);
                totalProb += Math.log(run.getProbability(converter)) * seqTuple.Item2;

            }

            return totalProb;
        }
    }


    public static void main(String[] args){
        long start = System.currentTimeMillis();
        GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});

        char[] nucleotides = {'A','C','T','G'};
        String[] gtTaxa = {"A","B","C","D","E"};
        char[][] sequences = new char[5][(int)(Math.pow(4, gtTaxa.length))];
        int index = 0;
        for (char a : nucleotides) {
            for (char b : nucleotides) {
                for (char c : nucleotides) {
                    for (char d : nucleotides) {
                        for (char e : nucleotides) {
                            sequences[0][index] = a;
                            sequences[1][index] = b;
                            sequences[2][index] = c;
                            sequences[3][index] = d;
                            sequences[4][index++] = e;
                        }
                    }
                }
            }
        }
        List<char[][]> allLoci = new ArrayList<char[][]>();
        allLoci.add(sequences);

        InferMLNetworkFromSequences inference = new InferMLNetworkFromSequences();
        inference.setParallel(3);
        inference.setStartingNetwork(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork("((((B:1)X#H1:0::0.2,A:1)i1:1,(X#H1:1::0.8,Y#H2:1::0.2)i4:1)i5:1,((D:1,E:1)i2:1,(C:1)Y#H2:1::0.8)i3:1)i6;"));
        System.out.println(inference.inferNetwork(gtTaxa, allLoci, null, gtrModel, 1, 0, 1));
        System.out.println((System.currentTimeMillis()-start)/1000.0);

    }

}
