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
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
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
public class InferILSNetworkParsimoniously extends MDCOnNetworkYFFromRichNewickJung {
    private Network[] _optimalNetworks;
    private double[] _optimalScores;

    public InferILSNetworkParsimoniously(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }

    public List<Tuple<Network,Double>> inferNetwork(List<Tree> gts, Map<String,List<String>> species2alleles, Long maxExaminations, int maxReticulations, int diameterLimit, Network startNetwork, int numSol){
        _optimalNetworks = new Network[numSol];
        _optimalScores = new double[numSol];
        Arrays.fill(_optimalScores, Integer.MAX_VALUE);

        List<Tree> distinctGTs = new ArrayList<Tree>();
        List<Double> gtCounter = new ArrayList<Double>();
        summarizeGeneTrees(gts, distinctGTs, gtCounter);

        DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = getStartNetwork(gts, species2alleles,startNetwork);
        NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
        AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);

        Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> scorer = getScoreFunction(distinctGTs, gtCounter, species2alleles);
        Comparator<Double> comparator = getDoubleScoreComparator();
        HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, maxExaminations, maxReticulations, diameterLimit); // search starts here
        //DirectedGraphToGraphAdapter<String,PhyloEdge<String>> resultNetwork = result.BestExaminedNetwork;
        //System.out.println("\n #Networks " + result.ExaminationsCount);
        List<Tuple<Network, Double>> resultTuples= postProcessResult(gts, species2alleles);
        return resultTuples;
    }

    private List<Tuple<Network, Double>> postProcessResult(List<Tree> gts, Map<String,List<String>> species2alleles){
        List<Tuple<Network, Double>> resultList = new ArrayList<Tuple<Network, Double>>();
        for(int i=0; i<_optimalNetworks.length; i++){
            MDCOnNetworkYF scorer = new MDCOnNetworkYF();
            scorer.computeInheritanceProb(_optimalNetworks[i], gts, species2alleles);
            resultList.add(new Tuple<Network, Double>(_optimalNetworks[i], _optimalScores[i]));
        }
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

    private void summarizeGeneTrees(List<Tree> originalGTs, List<Tree> distinctGTs, List<Double> counter){
        for(Tree tr: originalGTs){
            Double weight = ((STINode<Double>)tr.getRoot()).getData();
            if(weight == null){
                weight = 1.0;
            }
            int index = 0;
            boolean exist = false;
            for(Tree exTr: distinctGTs){
                if(Trees.haveSameRootedTopology(tr, exTr)){
                    exist = true;
                    break;
                }
                index++;
            }
            if(exist){
                counter.set(index, counter.get(index)+weight);
            }
            else{
                distinctGTs.add(tr);
                counter.add(weight);
            }
        }
    }

    private Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o2, o1);
            }
        };
    }


    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final List<Tree> gts, final List<Double> counters, final Map<String, List<String>> species2alleles){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                Network bniNetwork = networkNew2Old(network);

                MDCOnNetworkYF scorer = new MDCOnNetworkYF();
                List<Integer> scores = scorer.countExtraCoal(bniNetwork, gts, species2alleles);

                double total = 0;
                Iterator<Double> counter = counters.iterator();
                for(int score: scores){
                    total += score*counter.next();
                    //total += score;
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


}