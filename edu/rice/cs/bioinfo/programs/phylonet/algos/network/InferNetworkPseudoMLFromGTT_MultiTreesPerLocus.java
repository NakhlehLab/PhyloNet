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

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing.SimulatedAnnealingSalterPearL;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomTopologyNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferNetworkPseudoMLFromGTT_MultiTreesPerLocus extends InferNetworkPseudoMLFromGTT {
    public InferNetworkPseudoMLFromGTT_MultiTreesPerLocus(){
        _fullLikelihoodCalculator = new NetworkLikelihoodFromGTT_MultiTreesPerLocus();
        _likelihoodCalculator = new NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus();
        if(_batchSize!=0){
            ((NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus)_likelihoodCalculator).setBatchSize(_batchSize);
        }
    }

    protected void summarizeData(List originalData, Map<String,String> allele2species, List dataForStartingNetwork, List dataForInferNetwork, List dataCorrespondences){
        NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus likelihoodComputer = new NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus();
        likelihoodComputer.summarizeData(originalData, allele2species, dataForStartingNetwork, dataForInferNetwork, dataCorrespondences);

    }



    /*
    public void inferNetworkForBird(String gtFile, Map<String,String> allele2species, int maxReticulations, int numSol, LinkedList<Tuple<Network,Double>> resultList){

        List dataForNetworkInference = new ArrayList();
        List dataCorrespondence = new ArrayList();
        summarizeDataForBird(gtFile, dataForNetworkInference, dataCorrespondence);

        //System.out.println("Finish reading");
        Set<String> singleAlleleSpecies = new HashSet<>();
        if(_optimizeBL){
            _topologyVsParameterOperation[0] = 1;
            _topologyVsParameterOperation[1] = 0;
        }
        else{
            String[] taxa = {"CAPCA","PYGAD","BALRE","FALPE","PTEGU","MERNU","HALLE","PODCR","EGRGA","CATAU","APTFO","CALAN","PHACA","PHORU","TAUER","MANVI","PELCR","ACACH","COLST","LEPDI","CHLUN","MELUN","HALAL","TAEGU","COLLI","OPHHO","GALGA","MESUN","TYTAL","EURHE","CORBR","CARCR","BUCRH","APAVI","PHALE","CHAPE","NIPNI","MELGA","PICPU","GAVST","NESNO","FULGL","GEOFO","CHAVO","CUCCA","ANAPL"};
            for(String taxon: taxa){
                singleAlleleSpecies.add(taxon);
            }

        }

        for(Object nodeO: _startNetwork.dfs()){
            NetNode node = (NetNode)nodeO;
            for(Object parent: node.getParents()){
                if(Double.isNaN(node.getParentDistance((NetNode)parent)))
                    node.setParentDistance((NetNode)parent,1.0);
                if(node.isNetworkNode()){
                    if(Double.isNaN(node.getParentProbability((NetNode)parent)))
                        node.setParentProbability((NetNode)parent, 0.5);
                }
            }
        }

        String startingNetwork = _startNetwork.toString();

        NetworkRandomNeighbourGenerator allNeighboursStrategy = new NetworkRandomNeighbourGenerator(new NetworkRandomTopologyNeighbourGenerator(_topologyOperationWeight, maxReticulations, _moveDiameter, _reticulationDiameter, _fixedHybrid, _seed), _topologyVsParameterOperation[0], new NetworkRandomParameterNeighbourGenerator(singleAlleleSpecies), _topologyVsParameterOperation[1], _seed);
        Comparator<Double> comparator = getDoubleScoreComparator();
        //SimpleHillClimbing searcher = new SimpleHillClimbing(comparator, allNeighboursStrategy);

        SimulatedAnnealingSalterPearL searcher = new SimulatedAnnealingSalterPearL(comparator, allNeighboursStrategy, _seed);
        searcher.setLogFile(_logFile);
        searcher.setIntermediateFile(_intermediateResultFile.getAbsolutePath());
        Func1<Network, Double> scorer = getScoreFunctionForBird(dataForNetworkInference, dataCorrespondence);
        Network speciesNetwork = Networks.readNetwork(startingNetwork);
        searcher.search(speciesNetwork, scorer, numSol, _numRuns, _maxExaminations, _maxFailure, _optimizeBL, resultList); // search starts here

    }



    protected double computeLikelihoodForBird(final Network<Object> speciesNetwork, final List summarizedGTs, final List treeCorrespondence){
        //System.out.println(speciesNetwork.toString());
        NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus likelihoodComputer = new NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus();
        likelihoodComputer.setBatchSize(_batchSize);
        likelihoodComputer.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength, _Brent1, _Brent2, _numThreads);
        double prob = likelihoodComputer.computeProbabilityForBird(speciesNetwork, summarizedGTs, treeCorrespondence);
        return prob;
    }


    protected void summarizeDataForBird(String gtFile, List dataForInferNetwork, List dataCorrespondences){
        NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus likelihoodComputer = new NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus();
        likelihoodComputer.summarizeDataForBird4(gtFile, dataForInferNetwork, dataCorrespondences);
    }


    protected Func1<Network, Double> getScoreFunctionForBird(final List summarizedData, final List dataCorrespondences){
        return new Func1<Network, Double>() {
            public Double execute(Network speciesNetwork) {
                return computeLikelihoodForBird(speciesNetwork, summarizedData, dataCorrespondences);
            }
        };
    }
    */

}
