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
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;

import java.io.*;
import java.math.BigDecimal;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class InferNetworkML extends MDCOnNetworkYFFromRichNewickJung {
    protected Network[] _optimalNetworks;
    protected double[] _optimalScores;
    protected int _maxRounds;
    protected int _maxTryPerBranch;
    protected double _improvementThreshold;
    protected double _maxBranchLength;
    protected double _Brent1;
    protected double _Brent2;
    protected Long _maxFailure;
    protected Long _maxExaminations;
    protected int _diameterLimit;
    protected Network<Object> _startNetwork;
    protected Set<String> _fixedHybrid;
    protected double[] _operationWeight;
    protected int _numThreads;
    protected int _numRuns;
    protected Long _seed;

    protected File resultFile = null;
    protected File logFile = null;
    protected List<Tuple<Network,Double>> _networksTried = new ArrayList<Tuple<Network, Double>>();

    private boolean _printDetails = false;


/*
    protected Network[] _optimalNetworks;
    protected double[] _optimalScores;
    protected int _maxRounds = 100;
    protected int _maxTryPerBranch = 100;
    protected double _improvementThreshold = 0.001;
    protected double _maxBranchLength = 6;
    protected double _Brent1 = 0.01;
    protected double _Brent2 = 0.001;
    protected Long _maxExaminations = null;
    protected long _maxFailure = 100;
    protected int _diameterLimit;
    protected Network _startNetwork;
    protected double[] _operationWeight = {0.15,0.15,0.2,0.5};
    protected int _numRuns = 10;
    protected int _numThreads = 1;
    protected Long _seed = null;
    protected Set<String> _fixedHybrid = new HashSet<String>();
*/

    public void setResultFile(File file){
        resultFile = file;
    }

    public void setLogFile(File file){
        logFile = file;
    }


    public void setStartNetwork(Network startNetwork){
        _startNetwork = startNetwork;
    }


    public void setParallel(int numThreads){
        _numThreads = numThreads;
    }

    public InferNetworkML(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }

    public void setSearchParameter(int maxRounds, int maxTryPerBranch, double improvementThreshold, double maxBranchLength, double Brent1, double Brent2, Long maxExaminations, Long maxFailure, int diameterLimit, int parallel, Network startNetwork, Set<String> fixedHybrid, double[] operationWeight, int numRuns, Long seed){
        _maxRounds = maxRounds;
        _maxTryPerBranch = maxTryPerBranch;
        _improvementThreshold = improvementThreshold;
        _maxBranchLength = maxBranchLength;
        //_maxBranchLength = 12;
        _Brent1 = Brent1;
        _Brent2 = Brent2;
        _maxExaminations = maxExaminations;
        _diameterLimit = diameterLimit;
        _startNetwork = startNetwork;
        _maxFailure = maxFailure;
        _numThreads = parallel;
        _fixedHybrid = fixedHybrid;
        _operationWeight = operationWeight;
        _numRuns = numRuns;
        _seed = seed;
    }



    public void setPreviousResult(Network[] optimalNetworks, double[] optimalScores, int numRuns){
        _optimalNetworks = optimalNetworks;
        _optimalScores = optimalScores;
        _numRuns = numRuns;
    }


    public void readPreviousTriedNetwork(){
        try{
            if(logFile.exists()) {
                BufferedReader br = new BufferedReader(new FileReader(logFile));
                String line;
                while ((line=br.readLine())!=null){
                    String[] fields = line.split(" ");
                    _networksTried.add(new Tuple<Network, Double>(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork(fields[1].trim()),Double.parseDouble(fields[0])));
                }
                br.close();
            }
        }catch(Exception e) {
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
    }


    public List<Tuple<Network,Double>> inferNetwork(List originalData, Map<String,List<String>> species2alleles, int maxReticulations, int numSol){
        if(_optimalNetworks==null) {
            _optimalNetworks = new Network[numSol];
            _optimalScores = new double[numSol];
            Arrays.fill(_optimalScores, Double.NEGATIVE_INFINITY);
        }

        Map<String,String> allele2species = null;
        if(species2alleles!=null){
            allele2species = new HashMap<String, String>();
            for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                for(String allele: entry.getValue()){
                    allele2species.put(allele, entry.getKey());
                }
            }
        }

        List dataForNetworkInference = new ArrayList();
        List dataForStartingNetwork = new ArrayList();
        List dataCorrespondence = new ArrayList();
        summarizeData(originalData, allele2species, dataForStartingNetwork, dataForNetworkInference, dataCorrespondence);

        String startingNetwork = getStartNetwork(dataForStartingNetwork, species2alleles, _fixedHybrid, _startNetwork);
        dataForStartingNetwork.clear();



        for(int i=1; i<=_numRuns; i++){
            if(_printDetails) {
                System.out.println("\n\nRun #" + i + ":");
            }
            DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = makeNetwork(startingNetwork);
            NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(_operationWeight, makeNode, makeEdge, _seed);
            AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);

            Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> scorer = getScoreFunction(dataForNetworkInference, species2alleles, dataCorrespondence);
            Comparator<Double> comparator = getDoubleScoreComparator();
            HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, _maxExaminations, maxReticulations, _maxFailure, _diameterLimit); // search starts here

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
                }
            }

            if(_printDetails) {
                System.out.println("Results after run #" + i);
                for (int j = 0; j < _optimalNetworks.length; j++) {
                    System.out.println(_optimalScores[j] + ": " + _optimalNetworks[j].toString());
                }
                System.out.println("===============================\n\n");
            }

        }
        List<Tuple<Network, Double>> resultList = new ArrayList<Tuple<Network, Double>>();
        for(int i=0; i<numSol; i++){
            resultList.add(new Tuple(_optimalNetworks[i], _optimalScores[i]));
        }
        return resultList;
    }


    abstract protected String getStartNetwork(List dataForStartingNetwork, Map<String,List<String>> species2alleles, Set<String> hybridSpecies, Network<Object> startingNetwork);


    abstract protected void summarizeData(List originalData, Map<String,String> allele2species, List dataForStartingNetwork, List dataForInferNetwork, List dataCorrespondences);

    protected Tuple<Network,Double> hasAlreadyComputed(Network network){
        for(Tuple<Network,Double> networkTried: _networksTried){
            if(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.hasTheSameTopology(networkTried.Item1, network)){
                return networkTried;
            }
        }
        return null;
    }



    protected Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }


    protected Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final List summarizedData, final Map<String, List<String>> species2alleles, final List dataCorrespondences){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                //_counter++;
                Network<Object> speciesNetwork = networkNew2Old(network);
                Tuple<Network,Double> previousResult = hasAlreadyComputed(speciesNetwork);
                double score;
                if(previousResult == null) {
                    score = findOptimalBranchLength(speciesNetwork, species2alleles, summarizedData, dataCorrespondences);
                    String networkString = speciesNetwork.toString();
                    _networksTried.add(new Tuple<Network, Double>(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork(networkString), score));
                    if(logFile!=null) {
                        try {
                            BufferedWriter bw = new BufferedWriter(new FileWriter(logFile, true));
                            bw.append(score + " " + networkString + "\n");
                            bw.close();
                        } catch (Exception e) {
                            System.err.println(e.getMessage());
                            e.getStackTrace();
                        }
                    }
                }
                else{
                    speciesNetwork = previousResult.Item1;
                    score = previousResult.Item2;
                }
                //System.out.println(speciesNetwork + ": " + score);
                boolean optimalChanged = false;
                int originalIndex  = -1;
                if(score > _optimalScores[_optimalNetworks.length-1]){

                    for(int i=0; i<_optimalNetworks.length; i++){
                        if(_optimalNetworks[i]==null)break;
                        if(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.hasTheSameTopology(speciesNetwork, _optimalNetworks[i])){
                            originalIndex = i;
                            break;
                        }
                    }
                    if(originalIndex == -1 || score>_optimalScores[originalIndex]){
                        optimalChanged = true;
                        int index = -1;
                        for(int i=0; i<_optimalScores.length; i++){
                            if(score > _optimalScores[i]){
                                index = i;
                                break;
                            }
                        }
                        int startingMovingPoint = originalIndex == -1? _optimalScores.length-1:originalIndex;
                        for(int i=startingMovingPoint; i>index; i--){
                            _optimalNetworks[i] = _optimalNetworks[i-1];
                            _optimalScores[i] = _optimalScores[i-1];
                        }
                        _optimalScores[index] = score;
                        _optimalNetworks[index] = cloneNetwork(speciesNetwork);
                    }
                }

                System.gc();

                if(_printDetails) {
                    System.out.println(score + ": " + speciesNetwork.toString());
                    if (optimalChanged) {
                        System.out.println("Optimal ones:");
                        for (int i = 0; i < _optimalNetworks.length; i++) {
                            if (_optimalNetworks[i] != null) {
                                System.out.println(_optimalScores[i] + ": " + _optimalNetworks[i].toString());
                            }
                        }
                    }

                    System.out.println();
                }


                return score;
            }
        };
    }


    abstract protected double findOptimalBranchLength(final Network<Object> speciesNetwork, final Map<String, List<String>> species2alleles, final List summarizedData, final List dataCorrespondence);




    protected String network2String(final DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork){
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




    protected Network networkNew2Old(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork){
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



    protected Network cloneNetwork(Network network){
        return edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork(network.toString());
    }


}
