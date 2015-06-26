package edu.rice.cs.bioinfo.programs.phylonet.algos.network;

import java.util.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.*;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 5/10/12
 * Time: 3:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneTreeWithBranchLengthProbability {
    boolean _printDetail;
    private List<String> _netTaxa;
    private List<String> _stTaxa;
    private boolean [][] _S;
    private ArrayList<NodeInfo> _stNodesList;
    private HashMap<String,Integer> _nname2tamount;   //map the node name in the network to the number of corresponding nodes in the species tree
    private HashMap<String,String> _tname2nname;	 //map the node name in the species tree to the name of the corresponding node in the network
    private HashMap<String,List<TNode>> _hname2tnodes;  //map the name of hybrid node to the corresponding nodes in the species tree
    private Tree _st;
    private double[] _stAges;

    public GeneTreeWithBranchLengthProbability(){
        _printDetail = false;
        _netTaxa = new ArrayList<String>();
        _stTaxa = new ArrayList<String>();
        _stNodesList = new ArrayList<NodeInfo>();
        _nname2tamount = new HashMap<String,Integer>();   //map the node name in the network to the number of corresponding nodes in the species tree
        _tname2nname = new HashMap<String,String>();	 //map the node name in the species tree to the name of the corresponding node in the network
        _hname2tnodes = new HashMap<String,List<TNode>>();  //map the name of hybrid node to the corresponding nodes in the species tree

    }

    public void setPrintDetails(boolean p){
        _printDetail = p;
    }



    /**
     * The public function for calculating the probabilities.
     * @param	net 	the given network
     * @param 	gts		the given set of gene trees
     * @param	allele2species		the mapping from the names of allels to the names of the species. It is used for multiple alleles
     * @return	a list of probabilities corresponding to the list of gene trees.
     */
    public List<Double> calculateGTDistribution(Network<Double> net, List<Tree> gts, Map<String,String> allele2species){
        networkToTree(net);
        //System.out.println(_st.toNewickWD());
        //System.exit(0);
        //String[] leaves = {"A","B","C","D","E","F","G","H","I","J"};
        //gts = Trees.generateAllBinaryTrees(leaves);
        //System.out.println(gts.size());
        //System.exit(0);

        for(NetNode leaf: net.getLeaves()){
            _netTaxa.add(leaf.getName());
        }
        for(Map.Entry<String, Integer> entry: _nname2tamount.entrySet()){
            for(int i=1; i<=entry.getValue(); i++){
                String name = entry.getKey();
                _tname2nname.put(name + "_" + i, name);
            }
        }
        _S = calculateS(_st);
        generateSTNodesList();
        computeNodesUnderHybrid(_st);
        List<Double> probs = performCalculating(gts, allele2species);
        return probs;
    }



    /**
     * The actual calculation of the gene tree probabilities in the network
     * @param 	gts		the given set of gene trees
     * @param	allele2species		the mapping from the names of allels to the names of the species. It is used for multiple alleles
     * @return	a list of probabilities corresponding to the list of gene trees.
     */

    private List<Double> performCalculating(List<Tree> gts, Map<String,String> allele2species){
        List<Double> problist = new ArrayList<Double>();
        for(Tree gt: gts){
            //System.out.println(gt.toNewick()+": ");
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

            double[] gtNodeAges = getGTNodeAges(gt);
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

                if(_printDetail){
                    System.out.println();
                    for(int i=0; i<mapping.length; i++){
                        System.out.print(gtTaxa.get(i)+"->"+ _stTaxa.get(mapping[i])+"\t");
                    }
                    System.out.println();
                }
                Map<String,String> aname2tname = new HashMap<String,String>();
                for(int i = 0; i<mapping.length; i++){
                    aname2tname.put(gtTaxa.get(i), _stTaxa.get(mapping[i]));
                }

                int[] history = findCoalescentHistory(gt, gtNodeAges, aname2tname);
                if(history == null){
                    if(_printDetail){
                        System.out.println("0");
                        System.out.println();
                    }
                    continue;
                }

                double gtmapprob = 1;
                for(TNode b: _st.postTraverse()){
                    String nname = _tname2nname.get(b.getName());
                    if(_hname2tnodes.containsKey(nname)){
                        continue;
                    }

                    int u = calculateU(b, mapping, history);
                    if(u==0)continue;
                    List<Integer> nodelist = new ArrayList<Integer>();
                    int c = calculateC(b,history, nodelist);
                    double bprob = calculateBranchCoalProb(b, gtNodeAges, nodelist, u);
                    double gamma = ((STINode<Double>)b).getData();
                    gtmapprob *= bprob*Math.pow(gamma, u);
                    if(_printDetail){
                        if(gamma!=1 && u-c!=0){
                            if(u-c!=1){
                                System.out.print("*("+gamma+")^"+u);
                            }
                            else{
                                System.out.print("*("+gamma+")");
                            }
                        }
                    }
                }
                for(Map.Entry<String,List<TNode>> entry: _hname2tnodes.entrySet()){
                    int sumU = 0;
                    int sumC = 0;
                    List<Integer> nodelist = new ArrayList<Integer>();
                    TNode node = null;
                    for(TNode hnode: entry.getValue()){
                        if(node == null){
                            node = hnode;
                        }
                        int u = calculateU(hnode, mapping, history);
                        int c = calculateC(hnode, history, nodelist);
                        double gamma = ((STINode<Double>)hnode).getData();
                        gtmapprob *= Math.pow(gamma, u);
                        if(_printDetail){
                            if(gamma!=1 && u-c!=0){
                                if(u-c!=1){
                                    System.out.print("*("+gamma+")^"+u);
                                }
                                else{
                                    System.out.print("*("+gamma+")");
                                }
                            }
                        }
                        sumU += u;
                        sumC += c;
                    }

                    gtmapprob *= calculateBranchCoalProb(node, gtNodeAges, nodelist, sumU);
                }
                gtprob += gtmapprob;
                if(_printDetail)
                    System.out.println();

            }while(mergeNumberAddOne(mergeNumber,upper));


            if(_printDetail){
                System.out.println();
            }
            //System.out.println(gtprob);
            problist.add(gtprob);
        }
        //System.out.println("total probability:"+totalprob);
        return problist;
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
        _st = new STITree<Double>();
        ((STINode<Double>)(_st.getRoot())).setData(1.0);
        Queue<NetNode<Double>> source = new LinkedList<NetNode<Double>>();
        Queue<TMutableNode> dest = new LinkedList<TMutableNode>();
        source.offer(net.getRoot());
        dest.offer((TMutableNode) _st.getRoot());
        long nameid = System.currentTimeMillis();
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
                ((STINode<Double>)copy).setData(gamma);

                // Continue to iterate over the children of nn and tn.
                source.offer(child);
                dest.offer(copy);
            }
        }
    }


    /**
     * The function is to collect all nodes under hybridization so that they can be treated differently when calculating probabilities
     * @param	st	a tree
     */
    private void computeNodesUnderHybrid(Tree st){
        for(Map.Entry<String, Integer> entry: _nname2tamount.entrySet()){
            if(entry.getValue() > 1) {
                _hname2tnodes.put(entry.getKey(), new ArrayList<TNode>());
            }
        }

        for (TNode node : st.postTraverse()) {
            String name = _tname2nname.get(node.getName());
            List<TNode> nodelist = _hname2tnodes.get(name);
            if(nodelist!=null){
                nodelist.add(node);
            }
        }
    }


    private int[] findCoalescentHistory(Tree gt, double[] gtNodeAges, Map<String,String> allele2species){
        int[] history = new int[gt.getNodeCount()];
        Arrays.fill(history, -1);
        Map<Integer, BitSet> gtId2bs = new HashMap<Integer, BitSet>();

        for (TNode node : gt.postTraverse()) {
            BitSet gtBs = new BitSet(_stTaxa.size());
            if (node.isLeaf()) {
                String name = allele2species.get(node.getName());
                gtBs.set(_stTaxa.indexOf(name));
            }
            else {
                for (TNode child : node.getChildren()) {
                    BitSet childCluster = gtId2bs.get(child.getID());
                    gtBs.or(childCluster);
                }
                boolean found = false;
                double gtAge = gtNodeAges[node.getID()];
                for(NodeInfo ni: _stNodesList){
                    if(ni._bs.cardinality() >= gtBs.cardinality() && gtAge >= _stAges[ni._id]){
                        BitSet temp = (BitSet) gtBs.clone();
                        temp.and(ni._bs);
                        if(temp.equals(gtBs)){
                            if(ni._parentid == -1 || gtAge <= _stAges[ni._parentid]){
                                history[node.getID()] = ni._id;
                                found = true;
                                break;
                            }
                        }
                    }
                }
                if(!found){
                    return null;
                }
            }
            gtId2bs.put(node.getID(), gtBs);
        }
        //printHistory(history);
        return history;
    }

    /**
     * The function is to calculate the _R matrix for the given tree.
     * Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
     */

    private boolean[][] calculateS(Tree tree){
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


    private void generateSTNodesList(){
        _stAges = new double[_st.getNodeCount()];
        Map<Integer, BitSet> id2bs = new HashMap<Integer, BitSet>();
        for (TNode node : _st.postTraverse()) {
            BitSet bs = new BitSet(_stTaxa.size());
            int stId = node.getID();
            double age = Double.MAX_VALUE;
            if (node.isLeaf()) {
                String name = node.getName();
                bs.set(_stTaxa.indexOf(name));
                age = 0.0;
            }
            else {
                for (TNode child : node.getChildren()) {
                    age = Math.min(age, _stAges[child.getID()]+child.getParentDistance());
                    BitSet childCluster = id2bs.get(child.getID());
                    bs.or(childCluster);
                }
            }

            id2bs.put(stId, bs);
            _stAges[stId] = age;
            if(node.isRoot()) {
                _stNodesList.add(new NodeInfo(stId, bs, -1));
            } else {
                _stNodesList.add(new NodeInfo(stId, bs, node.getParent().getID()));
            }
        }

    }

    private double[] getGTNodeAges(Tree gt){
        double[] ages = new double[gt.getNodeCount()];

        for (TNode node : gt.postTraverse()) {
            Double gtAge = 0.0;
            if (!node.isLeaf()){
                double minTime = Double.MAX_VALUE;
                for(TNode child: node.getChildren()){
                    minTime = Math.min(ages[child.getID()]+child.getParentDistance()*2, minTime);
                }
                gtAge = minTime;
            }
            ages[node.getID()] = gtAge;
        }

        return ages;
    }



    /**
     * The function is to calculate the number of lineages going into a branch
     * @param	node	the node that the branch is incident into
     * @param	mapping		the mapping
     * @param	history		the coalescent history of the gene tree
     */
    private int calculateU(TNode node, int[] mapping, int[] history){
        int u = 0;
        for (int aMapping : mapping) {
            int mappingId = _st.getNode(_stTaxa.get(aMapping)).getID();
            if (node.isLeaf()) {
                if (node.getID() == mappingId) {
                    u++;
                }
            } else {
                if (_S[node.getID()][mappingId]) {
                    u++;
                }
            }
        }
        for(int aHistory : history){
            if(aHistory!=-1){
                if(_S[node.getID()][aHistory]){
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
    private int calculateC(TNode node, int[] history, List<Integer> nodelist){
        int c = 0;
        for(int i=0; i<history.length; i++){
            if(history[i]==node.getID()){
                nodelist.add(i);
                c++;
            }
        }
        return c;
    }


    private double calculateBranchCoalProb(TNode node, double[] gtNodeAges, List<Integer> nodelist, int u){
        double prob = 1;
        List<Double> coalTimes = new ArrayList<Double>();
        //System.out.println(node);
        double lowTao = _stAges[node.getID()];
        for(int id: nodelist){
            double gtAge = gtNodeAges[id];
            int j;
            for(j=0; j<coalTimes.size(); j++){
                if(gtAge < coalTimes.get(j)){
                    coalTimes.add(j, gtAge);
                    break;
                }
            }
            if(j == coalTimes.size()){
                coalTimes.add(gtAge);
            }
        }
        for(Double t: coalTimes){
            prob *= Math.exp((-1)*(t-lowTao)*u*(u-1)/2);
            if(_printDetail){
                System.out.print("*exp(-" + u*(u-1) + "*" + (t-lowTao) + "/2)");
            }
            lowTao = t;
            u--;
        }

        if(u != 1){
            //System.out.println(node.isRoot());
            double highTao = _stAges[node.getParent().getID()];
            prob *= Math.exp((-1)*(highTao-lowTao)*u*(u-1)/2);
            if(_printDetail){
                System.out.print("*exp(-" + u*(u-1) + "*" + (highTao-lowTao) + "/2)");
            }
        }
        //System.out.println();
        return prob;
    }


    /**
     * The function is to _printDetail matrix for debugging
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
     * The function is to _printDetail coalescent histories for debugging
     */
    private void printHistory(int[] history){
        System.out.print("[");
        for(int edge: history){
            System.out.print(edge+" ");
        }
        System.out.println("]");
    }


    private class NodeInfo{
        public int _id;
        public BitSet _bs;
        public int _parentid;

        public NodeInfo(int id, BitSet bs, int pid){
            _id = id;
            _bs = bs;
            _parentid = pid;
        }

        public String toString(){
            return _id + ":" + _bs;
        }
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
            double gamma = node.getParentProbability(parent) * child.getParentProbability(node);
            parent.removeChild(node);
            node.removeChild(child);
            parent.adoptChild(child, distance);
            child.setParentProbability(parent, gamma);
        }
    }
}
