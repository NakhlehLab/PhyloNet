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
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NetworkPseudoLikelihoodFromGTT extends NetworkLikelihood {
    private int _batchSize;

    public void setBatchSize(int size){
        _batchSize = size;
    }


    public static Map<String, double[]> computeTripleFrequenciesFromSingleGT(Tree gt){
        List<String> taxaList = new ArrayList<>();
        Map<String,Integer> pairwiseDepths = new HashMap<>();
        Map<TNode, Integer> node2depth = new HashMap<>();
        Map<TNode, List<String>> node2leaves = new HashMap<>();
        for(TNode node: gt.postTraverse()){
            int depth = 0;
            List<String> leaves = new ArrayList<>();
            if(node.isLeaf()){
                //leaves.add(taxon2ID.get(node.getName()));
                leaves.add(node.getName());
                taxaList.add(node.getName());
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
                            for(String leaf2: childLeaves2){
                                pairwiseDepths.put(leaf1+"&"+leaf2,depth);
                                pairwiseDepths.put(leaf2+"&"+leaf1,depth);
                            }
                        }
                    }

                }

            }
            node2depth.put(node, depth);
            node2leaves.put(node, leaves);
        }

        String[] taxa = taxaList.toArray(new String[0]);
        Arrays.sort(taxa);
        int numTaxa = taxa.length;
        Map<String, double[]> triple2counts = new HashMap<>();

        for(int i=0; i<numTaxa; i++){
            String taxonI = taxa[i];
            for(int j=i+1; j<numTaxa; j++){
                String taxonJ = taxa[j];
                String pair = taxonI+"&"+taxonJ;
                int ij = pairwiseDepths.get(pair);
                for(int k=j+1; k<numTaxa; k++){
                    String taxonK = taxa[k];
                    String triplet = pair+"&"+taxonK;
                    double[] freq = new double[3];
                    triple2counts.put(triplet, freq);
                    int minIndex = -1;
                    int ik = pairwiseDepths.get(taxonI+"&"+taxonK);
                    int jk = pairwiseDepths.get(taxonJ+"&"+taxonK);
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
                        freq[minIndex] = 1;
                    }
                    else{

                        for(int m=0; m<3; m++){
                            freq[m] = 1/3.0;
                        }
                    }
                }
            }
        }

        return triple2counts;
    }


    public static Map<String, double[]>  computeTripleFrequenciesFromSingleGT(Tree gt, Map<String,String> allele2species){
        Set<String> allAlleles = new HashSet<>();
        for(String allele: gt.getLeaves()){
            allAlleles.add(allele);
        }


        Map<TNode,Set<String>> node2leaves = new HashMap<>();
        Map<String, double[]> triple2counts = new HashMap<>();
        for(TNode node: gt.postTraverse()){
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
                            String species1 = allele2species.get(allele1);
                            for(String allele2: childLeaves2){
                                String species2 = allele2species.get(allele2);
                                if(!species1.equals(species2)){
                                    for(String allele3: allAlleles){
                                        String species3 = allele2species.get(allele3);
                                        if(!species1.equals(species3) && !species2.equals(species3)){
                                            addHighestFrequency(species1, species2, species3, triple2counts);
                                        }
                                    }
                                }
                            }
                        }
                        //non-binary node
                        for(int k=j+1; k<childLeavesList.size(); k++) {
                            Set<String> childLeaves3 = childLeavesList.get(k);
                            for(String allele1: childLeaves1) {
                                String species1 = allele2species.get(allele1);
                                for (String allele2 : childLeaves2) {
                                    String species2 = allele2species.get(allele2);
                                    if(!species1.equals(species2)) {
                                        for (String allele3 : childLeaves3) {
                                            String species3 = allele2species.get(allele3);
                                            if(!species1.equals(species3) && !species2.equals(species3)) {
                                                addEqualFrequency(species1, species2, species3, triple2counts);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                }

                allAlleles.addAll(leavesUnder);
            }

        }
        return triple2counts;
    }

    private static void addEqualFrequency(String species1, String species2, String species3, Map<String, double[]> triple2counts) {
        String[] forSort = new String[3];
        forSort[0] = species1;
        forSort[1] = species2;
        forSort[2] = species3;
        Arrays.sort(forSort);
        String exp = forSort[0]+"&"+forSort[1]+"&"+forSort[2];
        double[] frequency = triple2counts.get(exp);
        if(frequency == null) {
            frequency = new double[3];
            triple2counts.put(exp, frequency);
        }
        for (int i = 0; i < 3; i++) {
            frequency[i] += 1.0 / 3;
        }
    }



    private static void addHighestFrequency(String species1, String species2, String species3, Map<String, double[]> triple2counts) {
        int index = 0;
        String[] forSort = new String[3];
        forSort[0] = species1;
        forSort[1] = species2;
        forSort[2] = species3;
        Arrays.sort(forSort);
        String exp = forSort[0]+"&"+forSort[1]+"&"+forSort[2];
        if(forSort[2].equals(species3)){
            index = 0;
        }
        else if(forSort[1].equals(species3)){
            index = 1;
        }
        else if(forSort[0].equals(species3)){
            index = 2;
        }
        double[] frequency = triple2counts.get(exp);
        if(frequency == null) {
            frequency = new double[3];
            triple2counts.put(exp, frequency);
        }
        frequency[index] ++;
    }

    protected double findOptimalBranchLength(final Network<Object> speciesNetwork, final Map<String, List<String>> species2alleles, final List tripleFrequencies, final List gtCorrespondence, final Set<String> singleAlleleSpecies){
        boolean continueRounds = true; // keep trying to improve network

        for(NetNode<Object> node: speciesNetwork.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent,1.0);
                if(node.isNetworkNode()){
                    node.setParentProbability(parent, 0.5);
                }
            }
        }

        Set<NetNode> node2ignoreForBL = findEdgeHavingNoBL(speciesNetwork);

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
                    if(node2ignoreForBL.contains(child)){
                        continue;
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
        //System.out.println(speciesNetwork + " " + lnGtProbOfSpeciesNetwork.getContents());
        return lnGtProbOfSpeciesNetwork.getContents();
    }


    protected class MyThread extends Thread{
        Network _network;
        GeneTreeProbabilityPseudo _calculator;
        List<String> _allTriplets;
        double[][] _probs;


        public MyThread(Network network,GeneTreeProbabilityPseudo calculator,List<String> allTriplets, double[][] probs){
            _network = network;
            _calculator = calculator;
            _allTriplets = allTriplets;
            _probs = probs;
        }


        public void run(){
            _calculator.computePseudoLikelihood(_network, _allTriplets, _probs);
        }
    }




    protected double computeProbability(Network<Object> speciesNetwork, List allTriplets, Map<String, List<String>> species2alleles, List tripleFrequencies) {
        /*
        GeneTreeProbabilityPseudo likelihood = new GeneTreeProbabilityPseudo();
        double prob = likelihood.computePseudoLikelihood(speciesNetwork, tripleFrequencies);
        */
        GeneTreeProbabilityPseudo calculator = new GeneTreeProbabilityPseudo();
        if(_numThreads!=0){
            int batchSize = allTriplets.size()/_numThreads;
            if(allTriplets.size()%_numThreads != 0){
                batchSize++;
            }
            calculator.setBatchSize(batchSize);
        }
        calculator.initialize(speciesNetwork);
        double[][] probs = new double[allTriplets.size()][3];
        Thread[] myThreads = new Thread[_numThreads];
        //System.out.println(speciesNetwork);
        if(_numThreads>1) {
            calculator.setParallel(true);
            for (int i = 0; i < _numThreads; i++) {
                myThreads[i] = new MyThread(speciesNetwork, calculator, allTriplets, probs);
                myThreads[i].start();
            }
            for (int i = 0; i < _numThreads; i++) {
                try {
                    myThreads[i].join();
                } catch (InterruptedException ignore) {
                }
            }
        }else{
            try {
                calculator.computePseudoLikelihood(speciesNetwork, allTriplets, probs);
            }catch (Exception e){
                System.out.println(speciesNetwork);
                System.err.println(e.getMessage());
                e.getStackTrace();
                System.exit(-1);

            }
        }
        double totalProb = calculateFinalLikelihood(probs, tripleFrequencies);
        //System.out.println(speciesNetwork);
        //System.out.println(totalProb);
        return totalProb;
    }

    abstract protected double calculateFinalLikelihood(double[][] probs, List tripleFrequencies);


    public List gettingTripleFrequencies(Network speciesNetwork, List originalData, Map<String,List<String>> species2alleles){
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
        return summarizedData;
    }


    protected void findSingleAlleleSpeciesSet(Network speciesNetwork, Map<String,List<String>> species2alleles, Set<String> singleAlleleSpecies){
        for(Object node: speciesNetwork.getLeaves()){
            singleAlleleSpecies.add(((NetNode)node).getName());
        }
    }

    private Set<NetNode> findEdgeHavingNoBL(Network network) {
        Set<NetNode> node2ignore = new HashSet<>();
        Map<NetNode, Set<String>> node2leaves = new HashMap<>();
        for (Object nodeO : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeO;
            Set<String> leaves = new HashSet<>();
            if (node.isLeaf()) {
                leaves.add(node.getName());
            }
            else {
                for (Object childO : node.getChildren()) {
                    NetNode childNode = (NetNode) childO;
                    Set<String> childLeaves = node2leaves.get(childNode);
                    leaves.addAll(childLeaves);
                }
            }
            if(leaves.size()<=1){
                node2ignore.add(node);
            }
            node2leaves.put(node, leaves);

        }
        return node2ignore;
    }

}
