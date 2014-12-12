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
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
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
public class NetworkPseudoLikelihoodFromGTT extends NetworkLikelihood {

    protected double findOptimalBranchLength(final Network<Object> speciesNetwork, final Map<String, List<String>> species2alleles, final List tripleFrequencies, final List gtCorrespondence){
        boolean continueRounds = true; // keep trying to improve network

        for(NetNode<Object> node: speciesNetwork.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent,1.0);
                if(node.isNetworkNode()){
                    node.setParentProbability(parent, 0.5);
                }
            }
        }

        double initalProb = computeProbability(speciesNetwork, tripleFrequencies, species2alleles, gtCorrespondence);
        if(_printDetails)
            System.out.println(speciesNetwork.toString() + " : " + initalProb);

        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(initalProb);  // records the GTProb of the network at all times

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
                    if(child.isLeaf()){
                        if(species2alleles==null || species2alleles.get(child.getName()).size()<2){
                            continue;
                        }
                    }

                    assigmentActions.add(new Proc()
                    {
                        public void execute()
                        {

                            UnivariateFunction functionToOptimize = new UnivariateFunction() {
                                public double value(double suggestedBranchLength) {
                                    double incumbentBranchLength = child.getParentDistance(parent);

                                    child.setParentDistance(parent, suggestedBranchLength);

                                    double lnProb = computeProbability(speciesNetwork, tripleFrequencies, species2alleles, gtCorrespondence);
                                    //System.out.println(speciesNetwork + ": " + lnProb);
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

                            if(_printDetails)
                                System.out.println(speciesNetwork.toString() + " : " + lnGtProbOfSpeciesNetwork.getContents());

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
                                double incumbentHybridProbParent1 = child.getParentProbability(hybridParent1);

                                child.setParentProbability(hybridParent1, suggestedProb);
                                child.setParentProbability(hybridParent2, 1.0 - suggestedProb);

                                double lnProb = computeProbability(speciesNetwork, tripleFrequencies, species2alleles, gtCorrespondence);
                                //System.out.println(speciesNetwork + ": " + lnProb);
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
                        if(_printDetails)
                            System.out.println(speciesNetwork.toString() + " : " + lnGtProbOfSpeciesNetwork.getContents());
                    }
                });


            }


            // add hybrid probs to hybrid edges
            Collections.shuffle(assigmentActions);

            for(Proc assigment : assigmentActions)   // for each change attempt, perform attempt
            {
                assigment.execute();
            }
            if(_printDetails) {
                System.out.println("Round end ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                System.out.println(speciesNetwork.toString() + "\n" + lnGtProbOfSpeciesNetwork.getContents() + "\n");
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
        return lnGtProbOfSpeciesNetwork.getContents();
    }


    private class MyThread implements Callable<Double> {
        Network _network;
        GeneTreeProbabilityPseudo _calculator;
        List<MutableTuple<String, double[]>> _tripleFrequencies;


        public MyThread(Network network,GeneTreeProbabilityPseudo calculator,List<MutableTuple<String, double[]>> tripleFrequencies){
            _network = network;
            _calculator = calculator;
            _tripleFrequencies = tripleFrequencies;
        }


        public Double call(){
            //System.out.println("Start computing");
            return _calculator.computePseudoLikelihood(_network, _tripleFrequencies);
        }
    }




    protected double computeProbability(Network<Object> speciesNetwork, List tripleFrequencies, Map<String, List<String>> species2alleles, List gtCorrespondences) {
        /*
        GeneTreeProbabilityPseudo likelihood = new GeneTreeProbabilityPseudo();
        double prob = likelihood.computePseudoLikelihood(speciesNetwork, tripleFrequencies);
        */
        GeneTreeProbabilityPseudo calculator = new GeneTreeProbabilityPseudo();
        calculator.initialize(speciesNetwork);
        double totalProb = 0;
        //System.out.println(speciesNetwork);
        if(_numThreads>1) {
            calculator.setParallel(true);
            ExecutorService executor = Executors.newFixedThreadPool(_numThreads);
            List<Future<Double>> futureList = new ArrayList<Future<Double>>();
            for (int i = 0; i < _numThreads; i++) {
                MyThread myThread = new MyThread(speciesNetwork, calculator, tripleFrequencies);
                Future<Double> future = executor.submit(myThread);
                futureList.add(future);
            }
            for (Future<Double> future : futureList) {
                try {
                    totalProb += future.get();
                } catch (Exception e) {
                    System.err.println(e.getMessage());
                    e.getStackTrace();
                    System.exit(-1);
                }
            }
            executor.shutdown();
        }else{
            try {
                totalProb = calculator.computePseudoLikelihood(speciesNetwork, tripleFrequencies);
            }catch (Exception e){
                System.out.println(speciesNetwork);
                System.err.println(e.getMessage());
                e.getStackTrace();
                System.exit(-1);

            }
        }
        return totalProb;
    }



    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List tripleFrequencies, List treeCorrespondences){
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
        tripleFrequencies.addAll(computeTripleFrequenciesInGTs(gts, allele2species));
    }



    private List<MutableTuple<String, double[]>> computeTripleFrequenciesInGTs(List<MutableTuple<Tree,Double>> gts, Map<String,String> allele2species){
        String[] speciesArray;
        if(allele2species==null){
            speciesArray = gts.get(0).Item1.getLeaves();
        }
        else{
            Set<String> taxonSet = new HashSet<>();
            taxonSet.addAll(allele2species.values());
            speciesArray = new String[taxonSet.size()];
            int index = 0;
            for(String taxon: taxonSet){
                speciesArray[index++] = taxon;
            }
        }
        Map<String, double[]> triple2counts = new HashMap<>();
        for(int i=0; i<speciesArray.length; i++){
            for(int j=i+1; j<speciesArray.length; j++){
                for(int k=j+1; k<speciesArray.length; k++){
                    triple2counts.put(speciesArray[i]+"&"+speciesArray[j]+"&"+speciesArray[k], new double[3]);
                }
            }
        }
        Map<String, Integer> species2ID = new HashMap<>();
        for(int i=0; i<speciesArray.length; i++){
            species2ID.put(speciesArray[i], i);
        }

        for(MutableTuple<Tree,Double> gt: gts){
            if(allele2species==null) {

                computeTripleFrequenciesFromSingleGT(gt, speciesArray, species2ID, triple2counts);
            }else{
                Map<String, Integer> allele2speciesID = new HashMap<>();
                for(Map.Entry<String,String> entry: allele2species.entrySet()){
                    allele2speciesID.put(entry.getKey(), species2ID.get(entry.getValue()));
                }
                computeTripleFrequenciesFromSingleGT(gt, allele2speciesID, speciesArray, triple2counts);
            }
        }
        List<MutableTuple<String, double[]>> tripleFrequencies = new ArrayList<>();
        for(Map.Entry<String, double[]> entry: triple2counts.entrySet()){
            tripleFrequencies.add(new MutableTuple<String, double[]>(entry.getKey(), entry.getValue()));
        }
        return tripleFrequencies;
    }

    private void computeTripleFrequenciesFromSingleGT(MutableTuple<Tree,Double> gt, String[] taxa, Map<String, Integer> taxon2ID, Map<String, double[]> triple2counts){
        int numTaxa = taxon2ID.size();
        int[][] pairwiseDepths = new int[numTaxa][numTaxa];
        Map<TNode, Integer> node2depth = new HashMap<>();
        Map<TNode, List<String>> node2leaves = new HashMap<>();
        for(TNode node: gt.Item1.postTraverse()){
            int depth = 0;
            List<String> leaves = new ArrayList<>();
            if(node.isLeaf()){
                //leaves.add(taxon2ID.get(node.getName()));
                leaves.add(node.getName());
            }
            else{
                List<List<String>> childLeavesList = new ArrayList<>();
                for(TNode child: node.getChildren()){
                    depth = Math.max(depth, node2depth.get(child));
                    List<String> childLeaves = new ArrayList<>();
                    childLeaves.addAll(node2leaves.get(child));
                    childLeavesList.add(childLeaves);
                    leaves.addAll(childLeaves);
                }
                depth++;
                for(int i=0; i<childLeavesList.size(); i++){
                    List<String> childLeaves1 = childLeavesList.get(i);
                    for(int j=i+1; j<childLeavesList.size(); j++){
                        List<String> childLeaves2 = childLeavesList.get(j);
                        for(String leaf1: childLeaves1){
                            int id1 = taxon2ID.get(leaf1);
                            for(String leaf2: childLeaves2){
                                int id2 = taxon2ID.get(leaf2);
                                pairwiseDepths[id1][id2] = depth;
                                pairwiseDepths[id2][id1] = depth;
                            }
                        }
                    }

                }

            }
            node2depth.put(node, depth);
            node2leaves.put(node, leaves);
        }

        for(int i=0; i<numTaxa; i++){
            for(int j=i+1; j<numTaxa; j++){
                int ij = pairwiseDepths[i][j];
                String pair = taxa[i]+"&"+taxa[j];
                for(int k=j+1; k<numTaxa; k++){
                    int minIndex = -1;
                    int ik = pairwiseDepths[i][k];
                    int jk = pairwiseDepths[j][k];
                    if(ij<ik && ij<jk){
                        minIndex = 0;
                    }else if(ik<ij && ik<jk)
                    {
                        minIndex = 1;
                    }
                    else if(jk<ij && jk<ik){
                        minIndex = 2;
                    }
                    if(minIndex!=-1){
                        triple2counts.get(pair+"&"+taxa[k])[minIndex]+=gt.Item2;
                    }
                    else{
                        double[] frequencies = triple2counts.get(pair+"&"+taxa[k]);
                        for(int m=0; m<3; m++){
                            frequencies[m] += gt.Item2/3;
                        }
                    }
                }
            }
        }
    }


    private void computeTripleFrequenciesFromSingleGT(MutableTuple<Tree,Double> gt, Map<String,Integer> allele2speciesID, String[] speciesArray, Map<String, double[]> triple2counts){
        Set<String> allAlleles = new HashSet<>();
        int[] alleleNums = new int[speciesArray.length];
        for(String allele: gt.Item1.getLeaves()){
            allAlleles.add(allele);
            alleleNums[allele2speciesID.get(allele)]++;
        }
        Map<TNode,Set<String>> node2leaves = new HashMap<>();
        for(TNode node: gt.Item1.postTraverse()){
            Set<String> leavesUnder = new HashSet<>();
            node2leaves.put(node, leavesUnder);
            if(node.isLeaf()){
                leavesUnder.add(node.getName());
            }
            else{
                List<Set<String>> childLeavesList = new ArrayList<>();
                for(TNode child: node.getChildren()){
                    Set<String> childLeaves = node2leaves.get(child);
                    leavesUnder.addAll(childLeaves);
                    childLeavesList.add(childLeaves);
                }

                allAlleles.removeAll(leavesUnder);

                for(int i=0; i<childLeavesList.size(); i++){
                    Set<String> childLeaves1 = childLeavesList.get(i);
                    for(int j=i+1; j<childLeavesList.size(); j++){
                        Set<String> childLeaves2 = childLeavesList.get(j);
                        for(String allele1: childLeaves1){
                            int species1 = allele2speciesID.get(allele1);
                            for(String allele2: childLeaves2){
                                int species2 = allele2speciesID.get(allele2);
                                if(species1!=species2){
                                    for(String allele3: allAlleles){
                                        int species3 = allele2speciesID.get(allele3);
                                        if(species1!=species3 && species2!=species3){
                                            addHighestFrequency(species1, species2, species3, speciesArray, alleleNums, triple2counts);
                                        }
                                    }
                                }
                            }
                        }
                        //non-binary node
                        for(int k=j+1; k<childLeavesList.size(); k++) {
                            Set<String> childLeaves3 = childLeavesList.get(k);
                            for(String allele1: childLeaves1) {
                                int species1 = allele2speciesID.get(allele1);
                                for (String allele2 : childLeaves2) {
                                    int species2 = allele2speciesID.get(allele2);
                                    for (String allele3 : childLeaves3) {
                                        int species3 = allele2speciesID.get(allele3);
                                        addEqualFrequency(species1, species2, species3, speciesArray, alleleNums, triple2counts);
                                    }
                                }
                            }
                        }
                    }

                }

                allAlleles.addAll(leavesUnder);
            }

        }

    }

    private void addEqualFrequency(int species1, int species2, int species3, String[] speciesArray, int[] alleleNums, Map<String, double[]> triple2counts) {
        double weight = 1.0/(alleleNums[species1]*alleleNums[species2]*alleleNums[species3]);
        int[] forSort = new int[3];
        forSort[0] = species1;
        forSort[1] = species2;
        forSort[2] = species3;
        Arrays.sort(forSort);
        String exp = speciesArray[forSort[0]]+"&"+speciesArray[forSort[1]]+"&"+speciesArray[forSort[2]];
        double[] frequency = triple2counts.get(exp);
        for(int i=0; i<3; i++){
            frequency[i] += weight/3;
        }
    }



    private void addHighestFrequency(int species1, int species2, int species3, String[] speciesArray, int[] alleleNums, Map<String, double[]> triple2counts) {
        double weight = 1.0/(alleleNums[species1]*alleleNums[species2]*alleleNums[species3]);
        if(species1<species3 && species2<species3){
            String exp = speciesArray[Math.min(species1,species2)]+"&"+speciesArray[Math.max(species1, species2)]+"&"+speciesArray[species3];
            triple2counts.get(exp)[0] += weight;
        }
        else if((species1<species3&&species3<species2) || (species2<species3&&species3<species1)){
            String exp = speciesArray[Math.min(species1,species2)]+"&"+speciesArray[species3]+"&"+speciesArray[Math.max(species1,species2)];
            triple2counts.get(exp)[1] += weight;
        }
        else if(species1>species3 && species2>species3){
            String exp = speciesArray[species3]+"&"+speciesArray[Math.min(species1,species2)]+"&"+speciesArray[Math.max(species1,species2)];
            triple2counts.get(exp)[2] += weight;
        }
    }

    protected double calculateFinalLikelihood(double[] probs, List gtCorrespondences){return 0;}
}
