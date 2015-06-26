package edu.rice.cs.bioinfo.programs.phylonet.algos.network;


import edu.rice.cs.bioinfo.library.programming.MutableTuple;
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
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 3/13/12
 * Time: 11:31 AM
 * To change this template use File | Settings | File Templates.
 */

public class GeneTreeProbability {
    private boolean _printDetails = false;
    private List<String> _netTaxa;
    private List<String> _stTaxa;
    private boolean [][] _R, _M, _S;
    private Map<String,Integer> _nname2tamount;  //map the node name in the network to the number of corresponding nodes in the species tree
    private Map<String,String> _tname2nname;	 //map the node name in the species tree to the name of the corresponding node in the network
    private Map<String,List<TNode>> _hname2tnodes;  //map the name of hybrid node to the corresponding nodes in the species tree
    private Tree _mulTree;


    /**
     * Constructor that initialize the variables.
     */
    public GeneTreeProbability(){
        _netTaxa = new ArrayList<String>();
        _stTaxa = new ArrayList<String>();
        _nname2tamount = new TreeMap<String,Integer>();
        _hname2tnodes = new TreeMap<String,List<TNode>>();
        _tname2nname = new TreeMap<String,String>();
        _printDetails = false;
    }

    /**
     * kliu - paranoid.
     * Empty out state of this object. Just in case
     * calculateGTDistribution doesn't do a clean set of variables
     * over successive calls.
     */
    public void emptyState () {
	_netTaxa.clear();
	_stTaxa.clear();
	_nname2tamount.clear();
	_hname2tnodes.clear();
	_tname2nname.clear();
	_printDetails = false;
    }

    /**
     * The public function for calculating the probabilities.
     * @param	net 	the given network
     * @param 	gts		the given set of gene trees
     * @param	allele2species		the mapping from the names of allele to the names of the species. It is used for multiple alleles
     * @return	a list of probabilities corresponding to the list of gene trees.
     */
    public List<Double> calculateGTDistribution(Network<Double> net, List<Tree> gts, Map<String,String> allele2species, boolean toPrint){
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
                    System.out.print("Mapping: ");
                    for(int i=0; i<mapping.length; i++){
                        System.out.print(gtTaxa.get(i)+"->"+_stTaxa.get(mapping[i])+"\t");
                    }
                    System.out.println();
                }


                List<int[]> histories = computeHistories(gt, gtTaxa, mapping);

                double gtmapprob = 0;
                for(int[] history: histories){
                    gtmapprob += Double.parseDouble(computeProbability(mapping, history, false));
                }
                gtprob += gtmapprob;

                if(_printDetails){
                    System.out.println("Probability of this mapping: " + gtmapprob);
                    System.out.println();
                }


            }while(mergeNumberAddOne(mergeNumber,upper));

            if(_printDetails){
                System.out.println("Total probability of gene tree " + gt + " : " + gtprob);
            }
            problist.add(gtprob);
            totalprob += gtprob;
        }
        if(_printDetails){
            System.out.println();
            System.out.println("Total probability of all gene trees: "+totalprob);
        }
        return problist;
    }


    public List<Tuple3<int[],int[],Double>> findOptimalMapping(Network<Double> net, Tree gt, List<String> gtTaxa, Map<String,String> allele2species){
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

        List<Tuple3<int[],int[],Double>> results = new ArrayList<>();
        //for(Tree gt: gts){
            if(_printDetails){
                System.out.println("Gene tree " + gt+" :");
            }
            _R = calculateSorR(gt);
            //List<String> gtTaxa = Arrays.asList(gt.getLeaves());
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

            double maxProb = -1;
            //List<int[]> optimalMappings = new ArrayList<int[]>();
            //List<int[]> optimalHistories = new ArrayList<int[]>();
            MutableTuple<int[],int[]> optimalCoalescentHistory = new MutableTuple(null,null);
            //results.add(optimalCoalescentHistory);
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
                //TODO
                if(_printDetails){
                    for(int i=0; i<mapping.length; i++){
                        System.out.print(gtTaxa.get(i)+"->"+_stTaxa.get(mapping[i])+"\t");
                    }
                    System.out.println();
                }

                List<int[]> histories = computeHistories(gt, gtTaxa, mapping);

                for(int[] history: histories){
                    double prob = Double.parseDouble(computeProbability(mapping, history, false));
                    results.add(new Tuple3<int[], int[], Double>(mapping, history, prob));
                    /*
                    if(prob >= maxProb){
                        if(prob > maxProb){
                            maxProb = prob;
                            optimalMappings.clear();
                            optimalHistories.clear();
                        }
                        optimalMappings.add(mapping);
                        optimalHistories.add(history);                        
                    }
                    */

                    if(prob > maxProb){
                        maxProb = prob;
                        optimalCoalescentHistory.Item1 = mapping;
                        optimalCoalescentHistory.Item2 = history;
                    }

                }
                if(_printDetails)
                    System.out.println("");

            }while(mergeNumberAddOne(mergeNumber,upper));


               //if(_printDetails){
                   System.out.println(_mulTree);
                   System.out.println("Probability:" + maxProb);
                   //for(int i=0; i<optimalMappings.size(); i++){
                       int[] mapping = optimalCoalescentHistory.Item1;
                       int[] history = optimalCoalescentHistory.Item2;
                       System.out.println("Mapping:");
                       for(int j=0; j<mapping.length; j++){
                           System.out.print(gtTaxa.get(j)+"->"+_stTaxa.get(mapping[j])+"\t");
                       }
                       System.out.println();
                       System.out.println("History:");
                       for(int j=0; j<history.length; j++){
                           if(history[j] != -1){
                               System.out.println(gt.getNode(j).toString() + ":" + _mulTree.getNode(history[j]).getName());
                           }
                       }
                   //}
                   System.out.println();
               //}

            //allOptimalHistories.add(optimalHistories);
        //}

        return results;
    }

    /**
     * The public function for calculating the probabilities.
     * @param	net 	the given network
     * @param 	gts		the given set of gene trees
     * @param	allele2species		the mapping from the names of allels to the names of the species. It is used for multiple alleles
     * @return	a list of probabilities corresponding to the list of gene trees.
     */
    public List<Double> calculateExpectXL(Network<Double> net, List<Tree> gts, Map<String,String> allele2species){
        networkToTree(net);
        //System.out.println(st.toNewickWD());
        //System.exit(0);

        //String[] leaves = {"a","b1","b2","b3","c"};
        //gts = Trees.generateAllBinaryTrees(leaves);
        //System.out.println(gts.size());
        //System.exit(0);

        for(NetNode leaf: net.getLeaves()){
            _netTaxa.add(leaf.getName());
        }

        //
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

                //TODO

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

            //System.out.println();
            //TODO
            expectedXL = expectedXL / gtprob;
            //System.out.println(expectedXL);
            xllist.add(expectedXL);

        }

        return xllist;
    }


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


    private String computeProbability(int[] mapping, int[] history, boolean countXL){
        double gtmaphisprob = 1;
        boolean first = true;
        int xl = 0;

        for(TNode b: _mulTree.postTraverse()){
            String nname = _tname2nname.get(b.getName());
            if(nname != null){
                if(_hname2tnodes.containsKey(nname)){
                   continue;
                }
            }

            int u = calculateU(_mulTree, b, mapping, history);
            if(u==0)continue;
            int c = calculateC(b,history);
            double gij = gij(b.getParentDistance(),u,u-c);
            long w = calculateW(b,c,history);
            long d = calculateD(u,c);
            double gamma = ((STINode<Double>)b).getData();
            gtmaphisprob *= gij*w/d*Math.pow(gamma, u);
            if(countXL && b.getParentDistance()!=0){
                xl += Math.max(0, u-c-1);
            }
            if(_printDetails){
                //System.out.println("not h");
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
                distance = hnode.getParentDistance();
                int u = calculateU(_mulTree, hnode, mapping, history);
                int c = calculateC(hnode,history);
                double gamma = ((STINode<Double>)hnode).getData();
                double w = calculateHW(hnode,c,history);
                gtmaphisprob *= Math.pow(gamma, u);
                if(_printDetails){
                    String prefix = "*";
                    if(first){
                        prefix = "+";
                    }
                    if(gamma!=1 && u-c!=0){
                        if(u-c!=1){
                            System.out.print(prefix+"("+gamma+")^"+u);
                        }
                        else{
                            System.out.print(prefix+"("+gamma+")");
                        }
                        first = false;
                    }
                }
                sum_u += u;
                sum_c += c;
                prod_w *= w;
            }

            double gij = gij(distance,sum_u,sum_u-sum_c);
            long d = calculateD(sum_u,sum_c);
            prod_w *= fact(1,sum_c);
            gtmaphisprob *= gij*prod_w/d;
            if(countXL && distance!=0){
                xl += Math.max(0, sum_u-sum_c-1);
            }
            if(_printDetails){
                String prefix = "*";
                if(first){
                    prefix = "+";
                }
                if(gij!=1){
                    System.out.print(prefix+"g"+sum_u+(sum_u-sum_c)+"("+distance+")");
                    prefix = "*";
                    first = false;
                }
                if(d!=1){
                    System.out.print(prefix+"("+prod_w+"/"+d+")");
                    first = false;
                }
            }
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
     * The function is to convert a network to a multilabel tree.
     * @param	net 	the given network
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
        //long nameid = 0;
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
                }
                Integer amount = _nname2tamount.get(name);
                if(amount==null){
                    amount = 0;
                }
                _nname2tamount.put(name, ++amount);
                String newname = name + "_" + amount;
                copy = peer.createChild(newname);
                if(child.isLeaf()) {
                    _stTaxa.add(newname);
                }

                // Update the distance and data for this child.
                double distance = child.getParentDistance(parent);
                if (distance == NetNode.NO_DISTANCE) {
                    //copy.setParentDistance(TNode.NO_DISTANCE);
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
                //index ++;
            }
        }
    }


    /**
     * The function is to collect all nodes under hybridization so that they can be treated differently when calculating probabilities
     * @param	st	a tree
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
     * The function is to calculate the g_{ij} function.
     * @param	length	the branch length
     * @param 	i	the number of lineages in
     * @param	j	the number of lineages out
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
     * The function is to calculate factorial
     * @param	start	the first number
     * @param 	end		the last number
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
        //BigInteger ret = BigInteger.ONE;
        long ret = 1;
        for (int k = 0; k < K; k++) {
            //ret = ret.multiply(BigInteger.valueOf(N-k))
            //.divide(BigInteger.valueOf(k+1));
            ret = ret * (N-k) / (k+1);
        }
        //return ret.longValue();
        return ret;
    }


    /**
     * The function is to calculate the _R matrix for the given tree.
     * Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
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
     * The function is to calculate the number of lineages going into a branch
     * @param	st		the multilabel species tree
     * @param	node	the node that the branch is incident into
     * @param	mapping		the mapping
     * @param	history		the coalescent history of the gene tree
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
     * The function is to calculate the number of coalescent events in a branch
     * @param	node	the node that the branch is incident into
     * @param	history		the coalescent history
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
     * The function is to calculate the number of possible ordering coalescent events
     * @param	u	the number of lineages entering the branch
     * @param	c	the number of coalescent events
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
     * @param	c	the number of coalescent events on the branch
     * @param	history		coalescent history of the gene tree
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
     * @param	c	the number of  coalescent events on the branch
     * @param	history		coalescent history of the gene tree
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
     * The function is to _printDetails matrix for debugging
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
     * The function is to _printDetails a list of coalescent histories for debugging
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
     * The function is to _printDetails coalescent histories for debugging
     */
    private void printHistory(int[] history){
        System.out.print("[");
        for(int edge: history){
            System.out.print(edge+" ");
        }
        System.out.println("]");
    }

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
