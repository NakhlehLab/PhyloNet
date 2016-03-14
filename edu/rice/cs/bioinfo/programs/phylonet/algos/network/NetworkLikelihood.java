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
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
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

    //TODO
    protected int _maxNumACs;
    protected int _dataSize;


    public int getMaxNumACs(){
        return _maxNumACs;
    }

    public int getDataSize(){
        return _dataSize;
    }

    public void setParallel(int numThreads){
        _numThreads = numThreads;
    }

    public NetworkLikelihood(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }

    public void setSearchParameter(int maxRounds, int maxTryPerBranch, double improvementThreshold, double maxBranchLength, double Brent1, double Brent2, int numThreads){
        _maxRounds = maxRounds;
        _maxTryPerBranch = maxTryPerBranch;
        _improvementThreshold = improvementThreshold;
        _maxBranchLength = maxBranchLength;
        _Brent1 = Brent1;
        _Brent2 = Brent2;
        _numThreads = numThreads;
    }


    protected Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }



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


    protected double computeLikelihood(Network speciesNetwork, Map<String,List<String>> species2alleles, List summarizedData, List dataCorrespondences, Set<String> singleAlleleSpecies, boolean needOptimize){
        //long startTime = System.currentTimeMillis();
        double MLScore;
        if(!needOptimize){
            MLScore = computeProbability(speciesNetwork, summarizedData, species2alleles, dataCorrespondences);
        }
        else{
            MLScore = findOptimalBranchLength(speciesNetwork, species2alleles, summarizedData, dataCorrespondences, singleAlleleSpecies);
        }
        //System.out.println((System.currentTimeMillis()-startTime)/1000.0);
        return MLScore;
    }


    abstract protected void summarizeData(List originalData, Map<String,String> allele2species, List summarizedData, List dataCorrespondences);

    abstract protected double findOptimalBranchLength(Network<Object> speciesNetwork, Map<String, List<String>> species2alleles, List summarizedData, List dataCorrespondences, Set<String> singleAlleleSpecies);

    abstract protected double computeProbability(Network<Object> speciesNetwork, List summarizedData, Map<String, List<String>> species2alleles, List dataCorrespondences);

    abstract protected void findSingleAlleleSpeciesSet(Network speciesNetwork, Map<String,List<String>> species2alleles, Set<String> singleAlleleSpecies);

}
