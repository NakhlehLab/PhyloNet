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

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.*;

/**
 * Created by yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is to compute the likelihood of a species network
 */
public abstract class NetworkLikelihood extends MDCOnNetworkYFFromRichNewickJung {
    protected int _maxRounds = 100;
    protected int _maxTryPerBranch = 100;
    protected double _improvementThreshold = 0.001;
    protected double _maxBranchLength = 6;
    protected double _Brent1 = 0.01;
    protected double _Brent2 = 0.001;
    protected int _numThreads = 1;
    protected boolean _printDetails = false;
    protected int _numRuns = 5;


    /**
     * This function is to set the number of threads for parallel computing
     */
    public void setParallel(int numThreads){
        _numThreads = numThreads;
    }


    /**
     * This function is to set the number of threads for parallel computing
     */
    public NetworkLikelihood(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }


    /**
     * This function is to set all parameters used for optimizing branch lengths and inheritance probabilities of a species network when computing its likelihood
     *
     * @param maxRounds             the maximal rounds when using Brent's method to optimize branch lengths and inheritance probabilities of a network topology from gene trees
     * @param maxTryPerBranch       the maximal number trials of updating one branch length per round when using Brent's method to optimize branch lengths and inheritance probabilities of a network topology from gene trees
     * @param improvementThreshold  the threshold of likelihood improvement between rounds when using Brent's method to optimize branch lengths and inheritance probabilities of a network topology from gene trees
     * @param maxBranchLength       the upper bound of branch lengths when using Brent's method to optimize branch lengths and inheritance probabilities of a network topology from gene trees
     * @param Brent1                rel, which is original stopping criterion of Brent’s algorithm for optimizing branch lengths and inheritance probabilities of a network topology from gene trees
     * @param Brent2                abs, which is original stopping criterion of Brent’s algorithm for optimizing branch lengths and inheritance probabilities of a network topology from gene trees
     * @param numThreads            the number of threads for parallel computing
     */
    public void setSearchParameter(int maxRounds, int maxTryPerBranch, double improvementThreshold, double maxBranchLength, double Brent1, double Brent2, int numThreads){
        _maxRounds = maxRounds;
        _maxTryPerBranch = maxTryPerBranch;
        _improvementThreshold = improvementThreshold;
        _maxBranchLength = maxBranchLength;
        _Brent1 = Brent1;
        _Brent2 = Brent2;
        _numThreads = numThreads;
    }



    /**
     * This is the main function for computing the likelihood from the original data
     *
     * @param speciesNetwork        the species network
     * @param originalData          a collection of input data, gene trees or sequence alignments
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param needOptimize          whether optimizing the branch lengths and inheritance probabilities or using as is
     */
    public double computeLikelihood(Network speciesNetwork, List originalData, Map<String,List<String>> species2alleles, boolean needOptimize){
        Map<String,String> allele2species = null;
        if(species2alleles!=null){
            allele2species = new HashMap<String, String>();
            for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                for(String allele: entry.getValue()){
                    allele2species.put(allele, entry.getKey());
                }
            }
        }

        List dataCorrespondences = new ArrayList();
        List summarizedData = new ArrayList();
        summarizeData(originalData, allele2species, summarizedData, dataCorrespondences);

        Set<String> singleAlleleSpecies = new HashSet<>();
        findSingleAlleleSpeciesSet(speciesNetwork, species2alleles, singleAlleleSpecies);

        double maxProb = Double.NEGATIVE_INFINITY;
        for(int i=0; i<_numRuns; i++){
            double prob = computeLikelihood(speciesNetwork, species2alleles, summarizedData, dataCorrespondences, singleAlleleSpecies, needOptimize);
            maxProb = Math.max(maxProb, prob);
        }
        return maxProb;
    }


    /**
     * This function is to compute the likelihood from summarized data for one round
     *
     * @param speciesNetwork        the species network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param summarizedData        summarized data
     * @param dataCorrespondences   the correspondences between the summarized data and the original data
     * @param singleAlleleSpecies   to help identify which branch lengths in the species network can be ignored
     * @param needOptimize          whether optimizing the branch lengths and inheritance probabilities or using as is
     */
    protected double computeLikelihood(Network speciesNetwork, Map<String,List<String>> species2alleles, List summarizedData, List dataCorrespondences, Set<String> singleAlleleSpecies, boolean needOptimize){
        double MLScore;
        if(!needOptimize){
            MLScore = computeProbability(speciesNetwork, summarizedData, dataCorrespondences, species2alleles);
        }
        else{
            MLScore = findOptimalBranchLength(speciesNetwork, species2alleles, summarizedData, dataCorrespondences, singleAlleleSpecies);
        }
        return MLScore;
    }


    /**
     * This function is to summarize the input data
     *
     * @param originalData              original input data
     * @param allele2species            mapping from allele to species which it is sampled from
     * @param summarizedData            summarized data
     * @param dataCorrespondences       relationships between the original data and the data in dataForInferNetwork
     */
    public abstract void summarizeData(List originalData, Map<String, String> allele2species, List summarizedData, List dataCorrespondences);



    /**
     * This function is to optimize the branch lengths and inheritance probabilities of a given species network
     *
     * @param speciesNetwork            the species network
     * @param species2alleles           mapping from species to alleles which they is sampled from
     * @param summarizedData            summarized data
     * @param singleAlleleSpecies       to help identify which branch lengths in the species network can be ignored
     * @param dataCorrespondences       relationships between the original data and the data in dataForInferNetwork
     *
     * @return likelihood of the species network after its branch lengths and inheritance probabilities are optimized
     */
    abstract protected double findOptimalBranchLength(Network<Object> speciesNetwork, Map<String, List<String>> species2alleles, List summarizedData, List dataCorrespondences, Set<String> singleAlleleSpecies);



    /**
     * This function is to compute the likelihood from summarized data
     *
     * @param speciesNetwork        the species network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param summarizedData        summarized data
     * @param dataCorrespondences   the correspondences between the summarized data and the original data
     */
    abstract public double computeProbability(Network<Object> speciesNetwork, List summarizedData, List dataCorrespondences, Map<String, List<String>> species2alleles);


    /**
     * This function is to find the set of branches whose lengths cannot be estimated so that they can be ignored during the inference
     *
     * @param speciesNetwork        the species network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param singleAlleleSpecies   species that have only one allele sampled from it
     */
    abstract protected void findSingleAlleleSpeciesSet(Network speciesNetwork, Map<String,List<String>> species2alleles, Set<String> singleAlleleSpecies);

}
