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
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

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
public class InferILSNetworkParsimoniouslyParallel extends MDCOnNetworkYFFromRichNewickJung {
    private Network[] _optimalNetworks;
    private double[] _optimalScores;
    protected Long _maxFailure;
    protected Long _maxExaminations;
    protected int _diameterLimit;
    private Set<String> _fixedHybrid;
    protected Network<Object> _startNetwork;
    private int _numProcessors;
    protected double[] _operationWeight;
    private int _numRuns;
    private Long _seed;

    public InferILSNetworkParsimoniouslyParallel(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }


    public void setSearchParameter(Long maxExaminations, Long maxFailure, int diameterLimit, Network startNetwork, Set<String> fixedHybrid, int numProcessors, double[] operationWeight, int numRuns, Long seed){
        _maxExaminations = maxExaminations;
        _diameterLimit = diameterLimit;
        _startNetwork = startNetwork;
        _maxFailure = maxFailure;
        _fixedHybrid = fixedHybrid;
        _numProcessors = numProcessors;
        _operationWeight = operationWeight;
        _numRuns = numRuns;
        _seed = seed;
    }


    /*
    private void temp(List<Tree> gts){
        String[] leaves = {"Nyala2","LKudu1","GEland1","Eland3","Kudu4","MNyala1","Sitat2","Bongo1","Script1","Script2","Script3","Sylvat1","Sylvat2","Sylvat3","Sylvat4"};
        List<String> leafList = Arrays.asList(leaves);
        int index = 0;
        for(Tree tr: gts){
            ((STITree)tr).constrainByLeaves(leafList);
            System.out.println("\nTREE tre"+ index++ + "=" + tr.toString());
        }

    }
    */

    private void checkNetworkWithHybrids(Network<Object> startNetwork){
        if(startNetwork == null || _fixedHybrid.size() == 0){
            return;
        }

        if(startNetwork != null){
            for(NetNode<Object> node: startNetwork.getNetworkNodes()){
                NetNode hybridSpecies = node.getChildren().iterator().next();
                if(!(hybridSpecies.isLeaf() && _fixedHybrid.contains(hybridSpecies.getName()))){
                    throw new IllegalArgumentException("The starting network contains hybrid that is not in the specified hybrid set.");
                }
            }
        }
    }


    public List<Tuple<Network,Double>> inferNetwork(List<MutableTuple<Tree,Double>> gts, Map<String,List<String>> species2alleles, int maxReticulations, int numSol){
        _optimalNetworks = new Network[numSol];
        _optimalScores = new double[numSol];
        Arrays.fill(_optimalScores, Integer.MAX_VALUE);

        String startingNetwork = getStartNetwork(gts, species2alleles,_fixedHybrid, _startNetwork);

        for(int i=1; i<=_numRuns; i++){
            DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = makeNetwork(startingNetwork);
            NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(_operationWeight, makeNode, makeEdge, _seed);
            AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);
            //NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
            //AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);

            Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> scorer = getScoreFunction(gts, species2alleles);
            Comparator<Double> comparator = getDoubleScoreComparator();
            HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, _maxExaminations, maxReticulations, _maxFailure, _diameterLimit); // search starts here
            //HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, _maxExaminations, maxReticulations, _diameterLimit); // search starts here
/*
            System.out.println("Results after run #" + i);
            for(int j=0; j<_optimalNetworks.length; j++){
                System.out.println(_optimalScores[j] + ": " + network2String(_optimalNetworks[j]));
            }
            System.out.println("===============================\n\n");
*/
        }

        List<Tuple<Network, Double>> resultTuples= postProcessResult(gts, species2alleles);
        return resultTuples;
    }

    private List<Tuple<Network, Double>> postProcessResult(List<MutableTuple<Tree,Double>> gts, Map<String,List<String>> species2alleles){
        List<Tuple<Network, Double>> resultList = new ArrayList<Tuple<Network, Double>>();
        for(int i=0; i<_optimalNetworks.length; i++){
            MDCOnNetworkYF scorer = new MDCOnNetworkYF();
            int[] xls = new int[gts.size()];
            scorer.countExtraCoal(_optimalNetworks[i], gts, species2alleles, xls);
            resultList.add(new Tuple<Network, Double>(_optimalNetworks[i], _optimalScores[i]));
        }
        return resultList;
    }


    protected void createHybrid(Network<Object> network, String hybrid){
        List<Tuple<NetNode,NetNode>> edgeList = new ArrayList<Tuple<NetNode,NetNode>>();
        Tuple<NetNode,NetNode> destinationEdge = null;
        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(network)){
            for(NetNode child: node.getChildren()){
                if(child.isLeaf() && child.getName().equals(hybrid)){
                    if(node.isNetworkNode()){
                        return;
                    }
                    destinationEdge = new Tuple<NetNode, NetNode>(node, child);
                }
                else{
                    edgeList.add(new Tuple<NetNode, NetNode>(node, child));
                }
            }

        }

        int numEdges = edgeList.size();
        Tuple<NetNode,NetNode> sourceEdge = edgeList.get((int)(Math.random() * numEdges));
        NetNode insertedSourceNode = new BniNetNode();
        insertedSourceNode.adoptChild(sourceEdge.Item2, NetNode.NO_DISTANCE);
        sourceEdge.Item1.removeChild(sourceEdge.Item2);
        sourceEdge.Item1.adoptChild(insertedSourceNode, NetNode.NO_DISTANCE);
        NetNode insertedDestinationNode = new BniNetNode();
        insertedDestinationNode.adoptChild(destinationEdge.Item2, NetNode.NO_DISTANCE);
        destinationEdge.Item1.removeChild(destinationEdge.Item2);
        destinationEdge.Item1.adoptChild(insertedDestinationNode, NetNode.NO_DISTANCE);
        insertedSourceNode.adoptChild(insertedDestinationNode, NetNode.NO_DISTANCE);
    }


    private String getStartNetwork(List<MutableTuple<Tree,Double>> gts, Map<String,List<String>> species2alleles, Set<String> hybridSpecies, Network<Object> startingNetwork){
        checkNetworkWithHybrids(startingNetwork);

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
                sol = mdc.inferSpeciesTree(gts, false, 1, false, true, -1).get(0);
            }
            else{
                sol = mdc.inferSpeciesTree(gts, allele2species, false, 1, false, true, -1).get(0);
            }
            Tree startingTree= Trees.generateRandomBinaryResolution(sol._st);
            startingNetwork = string2Network(startingTree.toString());
        }

        for(String hybrid: hybridSpecies){
            createHybrid(startingNetwork, hybrid);
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

        return newNetwork;

    }


    private Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o2, o1);
            }
        };
    }



    private int[] computeXLParallel(Network network, List<MutableTuple<Tree,Double>> gts, Map<String, List<String>> species2alleles, int[] xls){
        Thread[] myThreads = new Thread[_numProcessors];

        MDCOnNetworkYF mdc = new MDCOnNetworkYF();
        mdc.setParallel(true);
        mdc.processNetwork(network);

        for(int i=0; i<_numProcessors; i++){
            myThreads[i] = new MyThread(network, gts, species2alleles, xls, mdc);
            myThreads[i].start();
        }

        for(int i=0; i<_numProcessors; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }

        return xls;

    }



    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final List<MutableTuple<Tree,Double>> gts, final Map<String, List<String>> species2alleles){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                Network bniNetwork = networkNew2Old(network);

                int[] scores = new int[gts.size()];
                if(_numProcessors == 1){
                    MDCOnNetworkYF scorer = new MDCOnNetworkYF();
                    scorer.countExtraCoal(bniNetwork, gts, species2alleles, scores);
                }
                else{
                    computeXLParallel(bniNetwork, gts, species2alleles, scores);
                }

                double total = 0;
                Iterator<MutableTuple<Tree,Double>> weightIt = gts.iterator();
                for(int score: scores){
                    total += score * weightIt.next().Item2;
                }


                if(total < _optimalScores[_optimalNetworks.length-1]){
                    boolean exist = false;
                    for(int i=0; i<_optimalNetworks.length; i++){
                        if(_optimalNetworks[i]==null)break;
                        if(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeClusterDistance(bniNetwork, _optimalNetworks[i])[2]<0.00001 &&
                                edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeTripartitionDistance(bniNetwork, _optimalNetworks[i])[2]<0.00001){
                            exist = true;
                            break;
                        }
                    }
                    if(!exist){
                        int index = -1;
                        for(int i=0; i<_optimalScores.length; i++){
                            if(total < _optimalScores[i]){
                                index = i;
                                break;
                            }
                        }
                        for(int i=_optimalScores.length-1; i>index; i--){
                            _optimalNetworks[i] = _optimalNetworks[i-1];
                            _optimalScores[i] = _optimalScores[i-1];
                        }
                        _optimalScores[index] = total;
                        //System.out.println(network2String(bniNetwork));
                        _optimalNetworks[index] = string2Network(network2String(bniNetwork));
                    }
                }

                /*
                System.out.println(total + ": "+network2String(bniNetwork));
                System.out.println(_optimalScores[0] + ": " + network2String(_optimalNetworks[0]));
                System.out.println();
                */
                return total;

            }
        };
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

    private class MyThread extends Thread{
        Network _speciesNetwork;
        List<MutableTuple<Tree,Double>> _geneTrees;
        Map<String, List<String>> _species2alleles;
        int[] _xls;
        MDCOnNetworkYF _mdc;


        public MyThread(Network speciesNetwork, List<MutableTuple<Tree,Double>> geneTrees, Map<String, List<String>> species2alleles, int[] xls, MDCOnNetworkYF mdc){
            _speciesNetwork = speciesNetwork;
            _geneTrees = geneTrees;
            _species2alleles = species2alleles;
            _xls = xls;
            _mdc = mdc;
        }


        public void run() {
            _mdc.countExtraCoal(_speciesNetwork, _geneTrees, _species2alleles, _xls);
        }
    }


}
