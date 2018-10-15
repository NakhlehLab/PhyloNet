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


import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/**
 * Created by Hunter Tidwell
 * Based on GeneTreeProbability by Yun Yu
 * Date: 10/5/17
 *
 * This class is to compute the probability of observing a collection of gene trees given a species network based on Mul-tree.
 * This method uses only topologies of gene trees.
 * See "The probability of a gene tree topology within a phylogenetic network with applications to hybridization detection”, PLoS Genetics, 2012
 * The method is an extension of "Gene tree distributions under the coalescent process", Evolution, 2005. Also see this paper for computation details. Variable names are the same as those appeared in the paper.
 */

public class GeneTreeProbabilityIntegrated {
    private boolean _printDetails = false;
    private List<String> _netTaxa;
    private List<String> _stTaxa;
    private boolean [][] _R, _M, _S;
    private Map<String,Integer> _nname2tamount;  //map the node name in the network to the number of corresponding nodes in the species tree
    private Map<String,String> _tname2nname;	 //map the node name in the species tree to the name of the corresponding node in the network
    private Map<String,List<TNode>> _hname2tnodes;  //map the name of (nodes below) hybrid node to the corresponding nodes in the species tree
    private Tree _mulTree;

    // Tau used as hyperparameter in coalescence likelihood integration
    private double _branch_length_truncation;
    // Beta used as hyperparameter in coalescence likelihood integration
    private double _branch_length_mean;

    // Count lineages in to a reticulation node from both parents
    // private Map<String, String> _hybrid_nodes;
    private Map<String, Integer> _hybrid_node_u_count; // keys are network nodes which have 2 parents
    private Map<String, Boolean> _hybrid_node_check; // True for names of nodes in the multree corresponding to hybrid nodes
    private Map<String, String> _hybrid_node_to_netnode; // map from multree to network for hybrid nodes
    private boolean _parallel;

    /**
     * Constructor that initializes the variables.
     */
    public GeneTreeProbabilityIntegrated(){
        _netTaxa = new ArrayList<String>();
        _stTaxa = new ArrayList<String>();
        _nname2tamount = new TreeMap<String,Integer>();
        _hname2tnodes = new TreeMap<String,List<TNode>>();
        _tname2nname = new TreeMap<String,String>();
        _printDetails = false;

        _hybrid_node_to_netnode = new HashMap<String, String>();
        _hybrid_node_u_count = new HashMap<String, Integer>();
        _hybrid_node_check = new HashMap<String, Boolean>();

        _branch_length_mean = 1.0;
        _branch_length_truncation = Double.POSITIVE_INFINITY;
    }

    /**
     * Cleans variables for potential subsequent calls
     */
    private void emptyState () {
        _netTaxa.clear();
        _stTaxa.clear();
        _nname2tamount.clear();
        _hname2tnodes.clear();
        _tname2nname.clear();
        _printDetails = false;

        _hybrid_node_to_netnode.clear();
        _hybrid_node_u_count.clear();
        _hybrid_node_check.clear();

        _branch_length_mean = 1.0;
        _branch_length_truncation = Double.POSITIVE_INFINITY;
    }

    /**
     * Sets the parameters of the branch length prior distribution
     */
    public void setBranchLengthExponentialPrior(double mean, double truncation) {
        // Default = 1
        _branch_length_mean = mean;
        // Default = Double.POSITIVE_INFINITY
        _branch_length_truncation = truncation;
    }

    /**
     * Get total log likelihood (assuming equal weights for each tree)
     * @param net
     * @param trees
     * @return
     */
    public double getTotalLogLikelihood(Network net, List<Tree> trees) {
        // todo pass in the allele map
        List<Double> probabilities = this.calculateGTDistribution(net, trees, null, false);
        double logLikelihood = 0.;
        for (Double d: probabilities) {
            logLikelihood += Math.log(d);
        }
        return logLikelihood;
    }

    /**
     * Computes the probability of observing a collection of gene trees given a species network
     *
     * @param net 	the given species network
     * @param gts	the given collection of gene trees
     * @param allele2species	the mapping from the names of alleles to the names of the species, which is used when multiple alleles are sampled per species
     * @param toPrint	whether to print out the details of computation
     *
     * @return	a list of probabilities corresponding to the list of gene trees.
     */
    public List<Double> calculateGTDistribution(Network<Double> net, List<Tree> gts, Map<String,String> allele2species, boolean toPrint) {

        _printDetails = toPrint;
        networkToTree(net);
        if(_printDetails){
            System.out.println("MUL tree: " + _mulTree.toNewick());
            System.out.println();
        }


        for(NetNode leaf: net.getLeaves()){
            _netTaxa.add(leaf.getName());
        }


        for(Map.Entry<String, Integer> entry: _nname2tamount.entrySet()){
            if(entry.getValue() > 1)
                for(int i=1; i<=entry.getValue(); i++){
                    String name = entry.getKey();
                    _tname2nname.put(name + "_" + i, name);
                }
        }
        _S = calculateSorR(_mulTree);
        computeNodesUnderHybrid(_mulTree);

        List<Double> problist = new ArrayList<Double>();
        double totalprob = 0;
        for(Tree oneGT: gts){
            if(_printDetails){
                System.out.println("Gene tree " + oneGT+" :");
            }

            List<Tree> allRootingGT = new ArrayList<>();
            if(oneGT.isRooted())
                allRootingGT.add(oneGT);
            else
                allRootingGT = oneGT.getAllRootingTrees();

            double gtprob = 0;

            for(Tree gt : allRootingGT) {

                _R = calculateSorR(gt);
                List<String> gtTaxa = Arrays.asList(gt.getLeaves());
                List<List<String>> allelesList = new ArrayList<List<String>>();
                for (int i = 0; i < _netTaxa.size(); i++) {
                    allelesList.add(new ArrayList<String>());
                }

                int[] upper = new int[_netTaxa.size()];
                for (String gtleaf : gtTaxa) {
                    String nleaf = gtleaf;
                    if (allele2species != null) {
                        nleaf = allele2species.get(gtleaf);
                    }
                    int index = _netTaxa.indexOf(nleaf);
                    List<String> alleles = allelesList.get(index);
                    if (alleles.size() == 0) {
                        upper[index] = _nname2tamount.get(nleaf);
                    }
                    alleles.add(gtleaf);
                }

                List<int[]> mergeNumber = new ArrayList<int[]>();
                for (List<String> alleles : allelesList) {
                    int[] first = new int[alleles.size()];
                    Arrays.fill(first, 1);
                    mergeNumber.add(first);
                }

                do {
                    int[] mapping = new int[gtTaxa.size()];
                    for (int i = 0; i < _netTaxa.size(); i++) {
                        String baseName = _netTaxa.get(i);
                        List<String> alleles = allelesList.get(i);
                        int[] subscribes = mergeNumber.get(i);
                        for (int j = 0; j < alleles.size(); j++) {
                            mapping[gtTaxa.indexOf(alleles.get(j))] = _stTaxa.indexOf(baseName + "_" + subscribes[j]);
                        }
                    }

                    if (_printDetails) {
                        System.out.print("Mapping: ");
                        for (int i = 0; i < mapping.length; i++) {
                            System.out.print(gtTaxa.get(i) + "->" + _stTaxa.get(mapping[i]) + "\t");
                        }
                        System.out.println();
                    }

                    List<int[]> histories = computeHistories(gt, gtTaxa, mapping);
                    double gtmapprob = 0;
                    for (int[] history : histories) {
                        gtmapprob += Double.parseDouble(computeProbability(mapping, history, false));
                    }
                    gtprob += gtmapprob;

                    if (_printDetails) {
                        System.out.println("Probability of this mapping: " + gtmapprob);
                        System.out.println();
                    }


                } while (mergeNumberAddOne(mergeNumber, upper));
            }
            if(_printDetails){
                System.out.println("Total probability of gene tree " + oneGT + " : " + gtprob);
            }
            problist.add(gtprob);
            totalprob += gtprob;
        }
        if(_printDetails){
            System.out.println();
            System.out.println("Total probability of all gene trees: "+totalprob);
        }
        emptyState();
        return problist;
    }


    /**
     * Finds the optimal coalescent history of a gene tree given a species network
     * @param net 	the species network
     * @param gt	the gene tree
     * @param gtTaxa	the set of taxa of the gene tree
     * @param allele2species	the mapping from the names of alleles to the names of the species, which is used when multiple alleles are sampled per species
     *
     * @return	the optimal coalescent history. The three-tuple contains the allele mapping, coalescent history and the probability, respectively
     */
    public Tuple3<int[],int[],Double> findOptimalMapping(Network<Double> net, Tree gt, List<String> gtTaxa, Map<String,String> allele2species){
        networkToTree(net);
        for(NetNode leaf: net.getLeaves()){
            _netTaxa.add(leaf.getName());
        }

        for(Map.Entry<String, Integer> entry: _nname2tamount.entrySet()){
            if(entry.getValue() > 1)
                for(int i=1; i<=entry.getValue(); i++){
                    String name = entry.getKey();
                    _tname2nname.put(name+"_"+i, name);
                }
        }
        _S = calculateSorR(_mulTree);
        computeNodesUnderHybrid(_mulTree);

        if(_printDetails){
            System.out.println("Gene tree " + gt+" :");
        }
        _R = calculateSorR(gt);
        List<List<String>> allelesList = new ArrayList<List<String>>();
        for(int i=0; i<_netTaxa.size(); i++){
            allelesList.add(new ArrayList<String>());
        }
        int[] upper = new int[_netTaxa.size()];
        for(String gtleaf: gtTaxa){
            String nleaf = gtleaf;
            if(allele2species!=null){
                nleaf = allele2species.get(gtleaf);
            }
            int index = _netTaxa.indexOf(nleaf);
            List<String> alleles = allelesList.get(index);
            if(alleles.size()==0){
                upper[index] = _nname2tamount.get(nleaf);
            }
            alleles.add(gtleaf);
        }

        List<int[]> mergeNumber = new ArrayList<int[]>();
        for(List<String> alleles: allelesList){
            int[] first = new int[alleles.size()];
            Arrays.fill(first, 1);
            mergeNumber.add(first);
        }

        Tuple3<int[],int[],Double> optimalCoalescentHistory = new Tuple3(null, null, -1);
        do{
            int[] mapping = new int[gtTaxa.size()];
            for(int i=0; i<_netTaxa.size(); i++){
                String baseName = _netTaxa.get(i);
                List<String> alleles = allelesList.get(i);
                int[] subscribes = mergeNumber.get(i);
                for(int j=0; j<alleles.size(); j++){
                    mapping[gtTaxa.indexOf(alleles.get(j))] = _stTaxa.indexOf(baseName+"_"+subscribes[j]);
                }
            }
            if(_printDetails){
                for(int i=0; i<mapping.length; i++){
                    System.out.print(gtTaxa.get(i)+"->"+_stTaxa.get(mapping[i])+"\t");
                }
                System.out.println();
            }

            List<int[]> histories = computeHistories(gt, gtTaxa, mapping);

            for(int[] history: histories){
                double prob = Double.parseDouble(computeProbability(mapping, history, false));

                if(prob > optimalCoalescentHistory.Item3){
                    optimalCoalescentHistory = new Tuple3(mapping, history, -1);
                }

            }
            if(_printDetails)
                System.out.println("");

        }while(mergeNumberAddOne(mergeNumber,upper));

        if(_printDetails) {
            System.out.println(_mulTree);
            System.out.println("Probability:" + optimalCoalescentHistory.Item3);
            int[] mapping = optimalCoalescentHistory.Item1;
            int[] history = optimalCoalescentHistory.Item2;
            System.out.println("Mapping:");
            for (int j = 0; j < mapping.length; j++) {
                System.out.print(gtTaxa.get(j) + "->" + _stTaxa.get(mapping[j]) + "\t");
            }
            System.out.println();
            System.out.println("History:");
            for (int j = 0; j < history.length; j++) {
                if (history[j] != -1) {
                    System.out.println(gt.getNode(j).toString() + ":" + _mulTree.getNode(history[j]).getName());
                }
            }
            System.out.println();
        }

        return optimalCoalescentHistory;
    }


    /**
     * Computes the expected number of extra lineages of a collection of gene trees given a species network
     * @param net 	the given species network
     * @param gts	the given collection of gene trees
     * @param allele2species	the mapping from the names of alleles to the names of the species, which is used when multiple alleles are sampled per species
     *
     * @return	a list of expected number of extra lineages corresponding to the list of gene trees.
     */
    public List<Double> calculateExpectXL(Network<Double> net, List<Tree> gts, Map<String,String> allele2species){
        networkToTree(net);
        for(NetNode leaf: net.getLeaves()){
            _netTaxa.add(leaf.getName());
        }

        for(Map.Entry<String, Integer> entry: _nname2tamount.entrySet()){
            if(entry.getValue() > 1)
                for(int i=1; i<=entry.getValue(); i++){
                    String name = entry.getKey();
                    _tname2nname.put(name+"_"+i, name);
                }
        }
        _S = calculateSorR(_mulTree);
        computeNodesUnderHybrid(_mulTree);

        List<Double> xllist = new ArrayList<Double>();
        for(Tree gt: gts){
            if(_printDetails){
                System.out.println("Gene tree " + gt+" :");
            }
            _R = calculateSorR(gt);
            List<String> gtTaxa = Arrays.asList(gt.getLeaves());
            List<List<String>> allelesList = new ArrayList<List<String>>();
            for(int i=0; i<_netTaxa.size(); i++){
                allelesList.add(new ArrayList<String>());
            }
            int[] upper = new int[_netTaxa.size()];
            for(String gtleaf: gtTaxa){
                String nleaf = gtleaf;
                if(allele2species!=null){
                    nleaf = allele2species.get(gtleaf);
                }
                int index = _netTaxa.indexOf(nleaf);
                List<String> alleles = allelesList.get(index);
                if(alleles.size()==0){
                    upper[index] = _nname2tamount.get(nleaf);
                }
                alleles.add(gtleaf);
            }

            List<int[]> mergeNumber = new ArrayList<int[]>();
            for(List<String> alleles: allelesList){
                int[] first = new int[alleles.size()];
                Arrays.fill(first, 1);
                mergeNumber.add(first);
            }

            double gtprob = 0;
            double expectedXL = 0;
            do{
                int[] mapping = new int[gtTaxa.size()];
                for(int i=0; i<_netTaxa.size(); i++){
                    String baseName = _netTaxa.get(i);
                    List<String> alleles = allelesList.get(i);
                    int[] subscribes = mergeNumber.get(i);
                    for(int j=0; j<alleles.size(); j++){
                        mapping[gtTaxa.indexOf(alleles.get(j))] = _stTaxa.indexOf(baseName+"_"+subscribes[j]);
                    }
                }

                if(_printDetails){
                    for(int i=0; i<mapping.length; i++){
                        System.out.print(gtTaxa.get(i)+"->"+_stTaxa.get(mapping[i])+"\t");
                    }
                    System.out.println();
                }

                List<int[]> histories = computeHistories(gt, gtTaxa, mapping);

                if(_printDetails){
                    System.out.println("Mapping:");
                    for(int j=0; j<mapping.length; j++){
                        System.out.print(gtTaxa.get(j)+"->"+_stTaxa.get(mapping[j])+"\t");
                    }
                    System.out.println();
                }

                for(int[] history: histories){
                    String result = computeProbability(mapping, history, true);
                    int index = result.indexOf("|");
                    double prob = Double.parseDouble(result.substring(0, index));
                    gtprob += prob;
                    int xl = Integer.parseInt(result.substring(index+1));
                    expectedXL += xl * prob;

                    if(_printDetails){
                        System.out.println("History:");
                        for(int j=0; j<history.length; j++){
                            if(history[j] != -1){
                                System.out.println(gt.getNode(j).toString() + ":" + _mulTree.getNode(history[j]).toString());
                            }
                        }
                        System.out.println(xl + ":" + prob);
                        System.out.println();
                    }
                }

                if(_printDetails)
                    System.out.println("");

            }while(mergeNumberAddOne(mergeNumber,upper));

            expectedXL = expectedXL / gtprob;
            xllist.add(expectedXL);
        }

        return xllist;
    }


    /**
     * Computes the coalescent histories of a gene tree under a given allele mapping
     * @param gt	the given gene tree
     * @param gtTaxa	the list of taxa in the gene tree
     * @param gtTaxa	the allele mapping
     *
     * @return	a list of coalescent histories
     */
    private List<int[]> computeHistories(Tree gt, List<String> gtTaxa, int mapping[]){
        List<int[]> histories = new ArrayList<int[]>();
        Map<String,String> aname2tname = new HashMap<String,String>();
        for(int i = 0; i<mapping.length; i++){
            aname2tname.put(gtTaxa.get(i), _stTaxa.get(mapping[i]));
        }
        calculateM(gt,_mulTree,aname2tname);
        Map<CEPair, Integer> ro = new HashMap<CEPair,Integer>();
        int[] his = new int[gt.getNodeCount()];
        Arrays.fill(his, -1);
        histories.add(his);
        enumCoalHistories(gt.getRoot(), ro, histories);
        ro.clear();
        return histories;
    }


    /**
     * Computes the probability of a coalescent history under a given allele mapping
     * @param mapping	the allele mapping
     * @param history   the coalescent history
     * @param countXL	whether this method is used for computing expected number of extra lineages; if so, the number of extra lineages will be computed in the meantime
     *
     * @return	the probability along with the number of extra lineages if countXL is set to true
     */
    private String computeProbability(int[] mapping, int[] history, boolean countXL){
        //return "0.5";
        double gtmaphisprob = 1;
        boolean first = true;
        int xl = 0;
        String hybrid_node_name;

        for(TNode b: _mulTree.postTraverse()){
            String nname = _tname2nname.get(b.getName());
            if(nname != null){
                if(_hname2tnodes.containsKey(nname)){ //todo: confirm These conditions are always equivalent
                    continue;
                }
            }

            int u = calculateU(_mulTree, b, mapping, history);
            if(u==0)continue;
            int c = calculateC(b,history);
            //double gij = gij(b.getParentDistance(),u,u-c);
            double gij = coalescence_from_scaled_exp(_branch_length_mean, u, u-c);

            if (b == _mulTree.getRoot()){
                // Coalescence in root always has probability 1
                gij = (u - c == 1) ? 1 : 0;
            }

            //System.out.println("\ngij params: " + b.getParentDistance() + ", " + u + ", " + (u - c) + "; gij old: " + gij1 + ", gij new: " + gij);
            long w = calculateW(b,c,history);
            long d = calculateD(u,c);
            double gamma = ((STINode<Double>)b).getData(); //todo: don't get the data
            gtmaphisprob *= gij * w/d; // * Math.pow(gamma, u);



            // update hybrid node values
            if (_hybrid_node_check.get(b.getName())) { //todo: maybe just check whether 2 lines down == null
                // need to store u for inheritance calculation
                 hybrid_node_name = _hybrid_node_to_netnode.get(b.getName());
                _hybrid_node_u_count.put(hybrid_node_name, _hybrid_node_u_count.get(hybrid_node_name) + u);
                gtmaphisprob *= fact(1, u + 1);
            }



            if(countXL && b.getParentDistance()!=0){
                xl += Math.max(0, u-c-1);
            }
            if(_printDetails){
                String prefix = "*";
                if(first){
                    prefix = "+";
                }
                if(gij!=1){
                    System.out.print(prefix+"g"+u+(u-c)+"("+b.getParentDistance()+")");
                    prefix = "*";
                    first = false;
                }
                if(d!=1){
                    System.out.print(prefix+"("+w+"/"+d+")");
                    prefix = "*";
                    first = false;
                }
                if(gamma!=1 && u!=0){
                    if(u!=1){
                        System.out.print(prefix+"("+gamma+")^"+u);
                    }
                    else{
                        System.out.print(prefix+"("+gamma+")");
                    }
                    first = false;
                }
            }
        }
        for(Map.Entry<String,List<TNode>> entry: _hname2tnodes.entrySet()){
            int sum_u = 0;
            int sum_c = 0;
            double prod_w =1;
            double distance = 0;
            for(TNode hnode: entry.getValue()){
                // distance = hnode.getParentDistance();
                int u = calculateU(_mulTree, hnode, mapping, history);
                int c = calculateC(hnode,history);
                //double gamma = ((STINode<Double>)hnode).getData();
                double w = calculateHW(hnode,c,history);
                //gtmaphisprob *= Math.pow(gamma, u); // could be done with sum_u
                if(_printDetails){
                    String prefix = "*";
                    if(first){
                        prefix = "+";
                    }
//                    if(gamma!=1 && u-c!=0){
//                        if(u-c!=1){
//                            System.out.print(prefix+"("+gamma+")^"+u);
//                        }
//                        else{
//                            System.out.print(prefix+"("+gamma+")");
//                        }
//                        first = false;
//                    }
                }
                sum_u += u;
                sum_c += c;
                prod_w *= w;
            }



            // update hybrid node values
            String one_of_the_multree_names = entry.getValue().get(1).getName(); // should always have at least 2 elements
            if (_hybrid_node_check.get(one_of_the_multree_names)) { //todo: maybe just check whether 2 lines down == null
                // need to store u for inheritance calculation
                hybrid_node_name = _hybrid_node_to_netnode.get(one_of_the_multree_names);
//                if (hybrid_node_name == null) {
//                    System.out.println(entry.toString());
//                    System.out.println(_hybrid_node_to_netnode.keySet());
//                }
                _hybrid_node_u_count.put(hybrid_node_name, _hybrid_node_u_count.get(hybrid_node_name) + sum_u);
                gtmaphisprob *= fact(1, sum_u + 1);
            }



            //double gij = gij(distance,sum_u,sum_u-sum_c);
            double gij = coalescence_from_scaled_exp(_branch_length_mean, sum_u, sum_u-sum_c);
            //System.out.println("\ngij params (below hybrid): " + distance + ", " + sum_u + ", " + (sum_u - sum_c) + "; gij old: " + gij1 + ", gij new: " + gij);
            long d = calculateD(sum_u,sum_c);
            prod_w *= fact(1,sum_c);
            gtmaphisprob *= gij*prod_w/d;
//            if(countXL && distance!=0){
//                xl += Math.max(0, sum_u-sum_c-1);
//            }
            if(_printDetails){
                String prefix = "*";
                if(first){
                    prefix = "+";
                }
                if(gij!=1){
//                    System.out.print(prefix+"g"+sum_u+(sum_u-sum_c)+"("+distance+")");
                    prefix = "*";
                    first = false;
                }
                if(d!=1){
                    System.out.print(prefix+"("+prod_w+"/"+d+")");
                    first = false;
                }
            }
        }

        // bottom half of inheritance probability
        for (String key: _hybrid_node_u_count.keySet()) {
            gtmaphisprob *= 6.0 / fact(1, _hybrid_node_u_count.get(key) + 3); // 6 = 1/Beta(2, 2)
            _hybrid_node_u_count.put(key, 0);
        }

        if(_printDetails)
            System.out.println(" = " + gtmaphisprob);
        if(countXL){
            return gtmaphisprob+"|"+xl;
        }
        else{
            return gtmaphisprob+"";
        }
    }


    /**
     * This function is to help enumerate allele mappings
     *
     * @return	false if it reaches the end
     */
    private boolean mergeNumberAddOne(List<int[]> mergeNumber, int[] upper){
        for(int i=0; i<mergeNumber.size(); i++){
            int[] partNumber = mergeNumber.get(i);
            int max = upper[i];
            for(int j=0; j<partNumber.length; j++){
                if(partNumber[j]==max){
                    partNumber[j] = 1;
                }
                else{
                    partNumber[j] = partNumber[j]+1;
                    return true;
                }
            }
            Arrays.fill(partNumber, 1);
        }
        return false;
    }


    /**
     * Converts a species network to a multi-labeled tree
     *
     * @param net 	the given network
     */
    private void networkToTree(Network<Double> net){
        removeBinaryNodes(net);
        _mulTree = new STITree<Double>();
        ((STINode<Double>)(_mulTree.getRoot())).setData(1.0);
        ((STINode<Double>)(_mulTree.getRoot())).setName("root");
        Queue<NetNode<Double>> source = new LinkedList<NetNode<Double>>();
        Queue<TMutableNode> dest = new LinkedList<TMutableNode>();
        source.offer(net.getRoot());
        dest.offer((TMutableNode) _mulTree.getRoot());
        int nameid = 1;

        _hybrid_node_check.put(_mulTree.getRoot().getName(), false);

        while(!source.isEmpty()){
            NetNode<Double> parent = source.poll();
            TMutableNode peer = dest.poll();
            for (NetNode<Double> child : parent.getChildren()) {
                TMutableNode copy;
                if (child.getName().equals(NetNode.NO_NAME)) {
                    child.setName("hnode" + (nameid++));
                }
                String name = child.getName();
                if(child.isNetworkNode()){
                    name = child.getName()+"TO"+parent.getName();
                    _hybrid_node_u_count.put(child.getName(), 0);
                }
                Integer amount = _nname2tamount.get(name); // so top-level hybrid nodes have count 1 but 2 copies.
                if(amount==null){
                    amount = 0;
                }
                _nname2tamount.put(name, ++amount);
                String newname = name + "_" + amount;
                copy = peer.createChild(newname);
                if(child.isLeaf()) {
                    _stTaxa.add(newname);
                }

                // Update the mappings for hybrid nodes
                _hybrid_node_check.put(copy.getName(), child.isNetworkNode()); // might be replacable with a check to ..._to_netnode
                if (child.isNetworkNode()) {
                    _hybrid_node_to_netnode.put(copy.getName(), child.getName());
                }

                // Update the distance and data for this child.
                double distance = child.getParentDistance(parent);
                if (distance == NetNode.NO_DISTANCE) {
                    copy.setParentDistance(0);
                }
                else {
                    copy.setParentDistance(distance);
                }

                double gamma = child.getParentProbability(parent);
                gamma = gamma==NetNode.NO_PROBABILITY?1.0:gamma;
                ((STINode<Double>)copy).setData(gamma);

                // Continue to iterate over the children of nn and tn.
                source.offer(child);
                dest.offer(copy);
            }
        }
    }


    /**
     * Collects all nodes under reticulation nodes which will be treated differently later when calculating probabilities
     * @param st the given species tree (MUL tree)
     */
    private void computeNodesUnderHybrid(Tree st){
        for(Map.Entry<String, Integer> entry: _nname2tamount.entrySet()){
            if (entry.getValue() > 1) {
                _hname2tnodes.put(entry.getKey(), new ArrayList<TNode>());
            }
        }

        for (TNode node : st.postTraverse()) {
            String name = _tname2nname.get(node.getName());

            if(name != null)
            {
                List<TNode> nodelist = _hname2tnodes.get(name);
                if(nodelist!=null){
                    nodelist.add(node);
                }
            }


        }
    }


    /**
     * Calculates the g_{ij} function.
     * @param	length	the branch length
     * @param 	i	the number of lineages in
     * @param	j	the number of lineages out
     *
     * @return	the resulting probability
     */
    private double gij(double length, int i, int j){
        if(length == TNode.NO_DISTANCE || length == -1){
            if(j == 1){
                return 1;
            }else{
                return 0;
            }
        }
        if(length==0){
            if(i==j)
                return 1;
            else
                return 0;
        }
        if(i==0){
            return 1;
        }
        double result = 0;
        for(int k=j; k<=i; k++){
            result += Math.exp(0.5*k*(1.0-k)*length)*(2.0*k-1.0)*Math.pow(-1,k-j)*fact(j,j+k-2)*fact(i-k+1,i)/(fact(1,j)*fact(1,k-j)*fact(i,i+k-1));
        }
        return result;
    }






    /**
     * Integrated version of g_{ij}
     *
     * Calculates the likelihood of observing u lineages coalescing into v,
     * assuming a truncated Exponential prior on length with support in (0, tau] and hyperparameter 1.
     */
    private double coalescence_from_truncated_exp(double tau, int u, int v) {

        double t1 = 1.0 / (1.0 - Math.exp(- tau));
        double t2 = 0.0;
        for (int j = v; j <= u; j++) {
            double f_term = calculate_f(u, v, j);
            double middle_term = 2.0 / (j * (j - 1) + 2.0);
            double last_term = 1.0 - Math.exp(-tau * (j * (j - 1) + 2) / 2);
            t2 += f_term * middle_term * last_term;
        }
        return t1 * t2;
    }

    /**
     * Replaces g_{ij} with a scaled exp prior
     * @param beta      parameter of the prior
     * @param u         number of lineages before coalescence / later in time
     * @param v         number of lineages after coalescence / earlier in time
     * @return          probability of this coalescence event
     */
    private double coalescence_from_scaled_exp(double beta, int u, int v) {
        double result = 0.0;
        for (int j = v; j <= u; j++) {
            double f_term = calculate_f(u, v, j);
            double scaled_term = 2.0 / (j * (j-1) * beta + 2);
            result += f_term * scaled_term;
//            if (f_term * scaled_term <= 0) {
//                System.out.println("f: " + f_term + ", scaled: " + scaled_term);
//            }
        }
        return result;
    }

    /**
     * Helper method for calculating coalescence probability
     */
    private double calculate_f(int u, int v, int j) {

        //refactor. based on g_{ij}

        // NOTE: the term (v+j-1) must be included. Do not wrap it into the fact(v,v+j-1) because this will give an incorrect factor of -1
        return (2.0*j-1.0) / (v+j-1) * Math.pow(-1,j-v) * fact(v,v+j-1) * fact(u-j+1,u) / (fact(1,v) * fact(1,j-v) * fact(u,u+j-1));
        //return (2.0*j-1.0) * Math.pow(-1,j-v) * fact(v,v+j-2) * (v+j-1) * fact(u-j,u) / (fact(1,v) * fact(1,j-v) * fact(u,u+j));
    }

    /**
     * Integrated version of gamma^u
     *
     * Calculates the likelihood of inheriting u1 lineages from first parent and u2 lineages from second parent,
     * assuming a Beta(2, 2) distribution on inheritance probability.
     */
    private double integrated_inheritance(int u1, int u2) {
        return fact(1, u1+1) * fact(1, u2+1) / (1.0 * fact(1, u1 + u2 + 3));
        // todo: some of the terms cancel
    }







    /**
     * The function is to calculate factorial
     * @param	start	the first number
     * @param 	end		the last number
     *
     * @return	the resulting factorial
     */
    private long fact(int start, int end){
        long result = 1;
        for(int i=start; i<=end; i++){
            result = result * i;
        }
        return result;
    }


    /**
     * The function is to calculate "N choose K"
     */
    private long choose(int N, int K) {
        long ret = 1;
        for (int k = 0; k < K; k++) {
            ret = ret * (N-k) / (k+1);
        }
        return ret;
    }


    /**
     * The function is to calculate the _R matrix for the given tree.
     * Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
     *
     * Note: used to efficiently find coalescent histories.
     */
    private boolean[][] calculateSorR(Tree tree){
        int nnode = tree.getNodeCount();

        boolean[][] matrix = new boolean[nnode][nnode];
        for(int i=0; i<nnode; i++){
            for(int j=0; j<nnode; j++){
                matrix[i][j] = false;
            }
        }
        Map<Integer, BitSet> map = new HashMap<Integer, BitSet>();

        for (TNode node : tree.postTraverse()) {
            BitSet bs = new BitSet(nnode);
            for(TNode child : node.getChildren()) {
                BitSet childBS = map.get(child.getID());

                bs.or(childBS);
            }
            int id = node.getID();
            for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i+1)) {
                matrix[id][i] = true;
            }
            bs.set(id);
            map.put(id, bs);

        }
        return matrix;
    }


    /**
     * The function is to calculate the _M matrix for the given tree.
     * Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
     *
     * @param gt the given gene tree
     * @param st the given species tree
     * @param allele2species mapping from the names of alleles to the names of the species, which is used when multiple alleles are sampled per species
     */
    private void calculateM(Tree gt, Tree st, Map<String,String> allele2species){
        int ngtnode = gt.getNodeCount();
        int nstnode = st.getNodeCount();
        _M = new boolean[nstnode][ngtnode];

        Map<Integer, BitSet> gtId2bs = new HashMap<Integer, BitSet>();
        List<Integer> rmlist = new ArrayList<Integer>();
        for (TNode node : gt.postTraverse()) {
            BitSet bs = new BitSet(_stTaxa.size());
            if (node.isLeaf()) {
                String name = allele2species.get(node.getName());
                bs.set(_stTaxa.indexOf(name));
                rmlist.add(node.getID());
            }
            else {
                for (TNode child : node.getChildren()) {
                    BitSet childCluster = gtId2bs.get(child.getID());
                    bs.or(childCluster);
                }
            }
            gtId2bs.put(node.getID(), bs);
        }

        for(Integer rmid: rmlist) {
            gtId2bs.remove(rmid);
        }
        rmlist.clear();

        Map<Integer, BitSet> stId2bs = new HashMap<Integer, BitSet>();
        for (TNode node : st.postTraverse()) {
            BitSet bs = new BitSet(_stTaxa.size());
            int stID = node.getID();
            if (node.isLeaf()) {
                String name = node.getName();
                bs.set(_stTaxa.indexOf(name));
            }
            else {
                for (TNode child : node.getChildren()) {
                    BitSet childCluster = stId2bs.get(child.getID());
                    bs.or(childCluster);
                }
            }
            stId2bs.put(stID, bs);

            for(Map.Entry<Integer, BitSet> entry : gtId2bs.entrySet()){
                int gtId = entry.getKey();
                BitSet gtBs = (BitSet) entry.getValue().clone();
                gtBs.and(bs);
                if(gtBs.equals(entry.getValue())){
                    rmlist.add(gtId);
                    _M[stID][gtId] = true;
                    for(int i=0; i< _S.length; i++){
                        if(_S[i][stID]){
                            _M[i][gtId] = true;
                        }
                    }
                }
            }
            for(int rmid: rmlist){
                gtId2bs.remove(rmid);
            }
        }
    }



    /**
     * The function is to enumerate coalescent histories.
     * Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
     */
    private void enumCoalHistories(TNode gtCluster, Map<CEPair,Integer> ro, List<int[]> histories){
        if(gtCluster.isLeaf()) {
            return;
        }
        // if this is the lowest cluster
        if(gtCluster.getLeafCount() <= 2) {
            List<int[]> temp = new ArrayList<int[]>();
            temp.addAll(histories);
            histories.clear();
            for(int edge=0; edge< _M.length; edge++){
                if(_M[edge][gtCluster.getID()]){
                    ro.put(new CEPair(gtCluster.getID(),edge),1);
                    for(int[] his: temp){
                        int[] newHis = his.clone();
                        newHis[gtCluster.getID()] = edge;
                        histories.add(newHis);
                    }
                }
            }
        } else {
            // compute children first
            for(TNode child : gtCluster.getChildren()) {
                enumCoalHistories(child, ro, histories);
            }

            // compute this cluster's score
            List<int[]> tempHis = new ArrayList<int[]>();
            tempHis.addAll(histories);
            histories.clear();

            for(int edge=0; edge< _M.length; edge++){
                if(!_M[edge][gtCluster.getID()]){
                    continue;
                }
                // compute the product of the sums
                CEPair cp = new CEPair();
                int sumProd = 1;
                List<int[]> childcoaledges = new ArrayList<int[]>();
                int[] coaledge = new int[_R.length];
                Arrays.fill(coaledge, -1);
                childcoaledges.add(coaledge);
                for(TNode child : gtCluster.getChildren()) {
                    if(sumProd == 0){
                        break;
                    }
                    // skip leaves since they aren't clusters
                    if(child.isLeaf()) {
                        continue;
                    }
                    List<int[]> tempcoaledges = new ArrayList<int[]>();
                    tempcoaledges.addAll(childcoaledges);
                    childcoaledges.clear();

                    int sum = 0;

                    for(int cedge=0; cedge< _M.length; cedge++){
                        if(_M[cedge][child.getID()]){
                            if(_S[edge][cedge] || edge==cedge){
                                cp.set(child.getID(),cedge);
                                sum += ro.get(cp);
                                for(int[] edges: tempcoaledges){
                                    int[] newedge = edges.clone();
                                    newedge[child.getID()] = cedge;
                                    childcoaledges.add(newedge);
                                }
                            }
                        }
                    }
                    tempcoaledges.clear();
                    sumProd *= sum;
                }
                for(int[] history: tempHis){
                    for(int[] childcoaledge: childcoaledges){
                        boolean add = true;
                        for(int i=0; i<childcoaledge.length; i++){
                            if(childcoaledge[i]==-1)continue;
                            if(history[i]!=childcoaledge[i]){
                                add = false;
                                break;
                            }
                        }
                        if(add){
                            int[] newHis = history.clone();
                            newHis[gtCluster.getID()] = edge;
                            histories.add(newHis);
                        }
                    }
                }
                childcoaledges.clear();
                ro.put(new CEPair(gtCluster.getID(),edge), Math.max(1, sumProd));
            }
            tempHis.clear();
        }
    }


    /**
     * The function is to calculate the number of lineages going into a branch for a given coalescent history
     *
     * @param	st		the multilabel species tree
     * @param	node	the node that the branch is incident into
     * @param	mapping		the mapping
     * @param	history		the given coalescent history of the gene tree
     */
    private int calculateU(Tree st, TNode node, int[] mapping, int[] history){
        int u = 0;
        for (int aMapping : mapping) {
            int mappingID = st.getNode(_stTaxa.get(aMapping)).getID();
            if (node.isLeaf()) {
                if (node.getID() == mappingID) {
                    u++;
                }
            } else {
                if (_S[node.getID()][mappingID]) {
                    u++;
                }
            }
        }
        for (int aHistory : history) {
            if (aHistory != -1) {
                if (_S[node.getID()][aHistory]) {
                    u--;
                }
            }
        }
        return u;
    }


    /**
     * The function is to calculate the number of coalescent events in a branch for a given coalescent history
     * @param node	the node that the branch is incident into
     * @param history	the coalescent history
     */
    private int calculateC(TNode node, int[] history){
        int c = 0;
        for (int aHistory : history) {
            if (aHistory == node.getID()) {
                c++;
            }
        }
        return c;
    }


    /**
     * The function is to calculate the number of all possible ordering coalescent events
     * @param	u	the number of lineages entering the branch
     * @param	c	the number of coalescent events happening on the branch
     */
    private long calculateD(int u, int c){
        long d = 1;
        if(c!=0){
            for(int i=1; i<=c; i++){
                d *= choose(u-i+1,2);
            }
        }
        return d;
    }


    /**
     * The function is to calculate the number of ways that coalescent events on a branch can occur consistently with the gene tree
     * @param node	the node that the branch is incident into
     * @param c	the number of coalescent events on the branch
     * @param history	coalescent history of the gene tree
     */
    private long calculateW(TNode node, int c, int[] history){
        long w = 1;
        if(c!=0){
            w = fact(1,c);
            for(int k=0; k<history.length; k++){
                if(history[k]==node.getID()){
                    int sum = 0;
                    for(int j=0; j<history.length; j++){
                        if(j==k || history[j]==-1)continue;
                        if((history[j]==node.getID() || _S[history[j]][node.getID()]) && _R[k][j]){
                            sum ++;
                        }
                    }
                    w *= 1.0/(1+sum);
                }
            }
        }
        return w;
    }


    /**
     * The function is to calculate the number of ways that coalescent events on a branch can occur consistently with the gene tree divided by c!
     * It is used for branches under hybridization events.
     * @param node	the node that the branch is incident into
     * @param c	the number of  coalescent events on the branch
     * @param history	coalescent history of the gene tree
     */
    private double calculateHW(TNode node, int c, int[] history){
        double w = 1;
        if(c!=0){
            for(int k=0; k<history.length; k++){
                if(history[k]==node.getID()){
                    int sum = 0;
                    for(int j=0; j<history.length; j++){
                        if(j==k || history[j]==-1)continue;
                        if((history[j]==node.getID() || _S[history[j]][node.getID()]) && _R[k][j]){
                            sum ++;
                        }
                    }
                    w /= 1+sum;
                }
            }
        }
        return w;
    }


    /**
     * The function is to print a matrix for debugging
     */
    private void printMatrix(boolean[][] matrix){
        for(int i=0; i<matrix.length; i++){
            for(int j=0; j<matrix[0].length; j++){
                if(matrix[i][j])
                    System.out.print(1 + "\t");
                else
                    System.out.print(0 + "\t");
            }
            System.out.println();
        }
        System.out.println();
    }


    /**
     * The function is to print a list of coalescent histories for debugging
     */
    private void printHistories(List<int[]> histories){
        System.out.println("total size:"+histories.size());
        for(int[] his: histories){
            System.out.print("[");
            for(int edge: his){
                System.out.print(edge+" ");
            }
            System.out.println("]");
        }
    }

    /**
     * The function is to print a coalescent history for debugging
     */
    private void printHistory(int[] history){
        System.out.print("[");
        for(int edge: history){
            System.out.print(edge+" ");
        }
        System.out.println("]");
    }


    /**
     * The function is to remove binary nodes of a species network
     */
    private void removeBinaryNodes(Network<Double> net)
    {
        // Find all binary nodes.
        List<NetNode<Double>> binaryNodes = new LinkedList<NetNode<Double>>();
        for (NetNode<Double> node : net.bfs()) {
            if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                binaryNodes.add(node);
            }
        }

        // Remove them.
        for (NetNode<Double> node : binaryNodes) {
            NetNode<Double> child = node.getChildren().iterator().next();	// Node's only child.
            if(child.getIndeg() != 1){
                continue;
            }
            NetNode<Double> parent = node.getParents().iterator().next();	// Node's only parent.
            double distance = node.getParentDistance(parent) + child.getParentDistance(node);
            double gamma1 = node.getParentProbability(parent)==NetNode.NO_PROBABILITY?1.0:node.getParentProbability(parent);
            double gamma2 = child.getParentProbability(node)==NetNode.NO_PROBABILITY?1.0:child.getParentProbability(node);
            double gamma =  gamma1 * gamma2;
            parent.removeChild(node);
            node.removeChild(child);
            parent.adoptChild(child, distance);
            child.setParentProbability(parent, gamma);
        }
    }


    /**
     * The class is to store cluster - edge pair.
     */
    private class CEPair {

        public int _clusterID;
        public int _edgeID;

        public CEPair(){}

        public CEPair(int cluster, int edge) {
            _edgeID = edge;
            _clusterID = cluster;
        }

        public void set(int cluster, int edge) {
            _edgeID = edge;
            _clusterID = cluster;
        }

        public int hashCode() {
            return _edgeID;
        }

        public boolean equals(Object o) {
            if(!(o instanceof CEPair)) {
                return false;
            }

            CEPair p2 = (CEPair) o;

            return (_clusterID == p2._clusterID) && (_edgeID == p2._edgeID);
        }

        public String toString(){
            return "edge:"+ _edgeID +"/node:"+ _clusterID;
        }
    }

}
