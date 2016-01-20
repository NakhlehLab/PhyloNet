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

import edu.rice.cs.bioinfo.library.programming.Container;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus extends NetworkPseudoLikelihoodFromGTT {

    public void summarizeData(List originalGTs, Map<String,String> allele2species, List allTriplets, List tripletFrequencies){
        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        for(Object list: originalGTs) {
            for (MutableTuple<Tree, Double> gtTuple : (List<MutableTuple<Tree, Double>>)list) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
                if (existingTuple == null) {
                    existingTuple = gtTuple;
                    exp2tree.put(exp, existingTuple);

                } else {
                    existingTuple.Item2 += gtTuple.Item2;
                }

            }
        }
        List<MutableTuple<Tree,Double>> gts = new ArrayList<>();
        gts.addAll(exp2tree.values());
        computeTripleFrequenciesInGTs(gts, allele2species, allTriplets, tripletFrequencies);
    }


    public void summarizeData(List originalGTs, Map<String,String> allele2species, List distinctGTs, List allTriplets, List tripletFrequencies){
        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        for(Object list: originalGTs) {
            for (MutableTuple<Tree, Double> gtTuple : (List<MutableTuple<Tree, Double>>)list) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
                if (existingTuple == null) {
                    existingTuple = gtTuple;
                    exp2tree.put(exp, existingTuple);

                } else {
                    existingTuple.Item2 += gtTuple.Item2;
                }

            }
        }
        distinctGTs.addAll(exp2tree.values());
        computeTripleFrequenciesInGTs(distinctGTs, allele2species, allTriplets, tripletFrequencies);
    }


    public static void computeTripleFrequenciesInGTs(List<MutableTuple<Tree,Double>> gts, Map<String,String> allele2species, List<String> allTriplets, List<double[]> tripletFrequencies){

        Map<String, double[]> triple2counts = new HashMap<>();

        for(MutableTuple<Tree,Double> gt: gts){
            Map<String, double[]> result;
            if(allele2species==null) {
                result = computeTripleFrequenciesFromSingleGT(gt.Item1);

            }else{
                result = computeTripleFrequenciesFromSingleGT(gt.Item1, allele2species);
            }
            for(Map.Entry<String,double[]> entry: result.entrySet()) {
                double[] freq = triple2counts.get(entry.getKey());
                if (freq == null) {
                    freq = new double[3];
                    triple2counts.put(entry.getKey(), freq);
                }
                for (int i = 0; i < freq.length; i++) {
                    freq[i] += entry.getValue()[i] * gt.Item2;
                }

            }
        }

        for(Map.Entry<String, double[]> entry: triple2counts.entrySet()){
            allTriplets.add(entry.getKey());
            tripletFrequencies.add(entry.getValue());
            /*
            for(double value: entry.getValue()){
                System.out.print(value/5000.0 + ",");
            }
            */
        }

    }





    protected double calculateFinalLikelihood(double[][] probs, List tripletFrequencies){
        int index = 0;
        double totalProb = 0;
        for(Object o: tripletFrequencies){
            double[] freqs = (double[])o;
            for(int i=0; i<3; i++){
                if(probs[index][i] == 0)return Double.NEGATIVE_INFINITY;
                totalProb += freqs[i] * Math.log(probs[index][i]);
                //totalProb *= Math.pow(probs[index][i],freqs[i]);
            }
            index++;
        }
        return totalProb;
    }

}
