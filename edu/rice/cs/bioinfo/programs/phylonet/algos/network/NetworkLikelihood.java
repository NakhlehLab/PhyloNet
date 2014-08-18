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
    protected int _maxRounds;
    protected int _maxTryPerBranch;
    protected double _improvementThreshold;
    protected double _maxBranchLength;
    protected double _Brent1;
    protected double _Brent2;
    protected int _numThreads;


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

    abstract protected void summarizeData(List originalData, Map<String,String> allele2species, List summarizedData, List dataCorrespondences);

    abstract protected double findOptimalBranchLength(Network<Object> speciesNetwork, Map<String, List<String>> species2alleles, List summarizedData, List dataCorrespondences);

    abstract protected double computeProbability(Network<Object> speciesNetwork, List summarizedData, Map<String, List<String>> species2alleles, List dataCorrespondences);


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
        return computeLikelihood(speciesNetwork, summarizedData, dataCorrespondences, species2alleles, needOptimize);
    }


    protected double computeLikelihood(Network speciesNetwork,  List summarizedData, List dataCorrespondences, Map<String,List<String>> species2alleles, boolean needOptimize){
        double MLScore;
        if(!needOptimize){
            MLScore = computeProbability(speciesNetwork, summarizedData, species2alleles, dataCorrespondences);
        }
        else{
            MLScore = findOptimalBranchLength(speciesNetwork, species2alleles, summarizedData, dataCorrespondences);
        }
        return MLScore;
    }


}
