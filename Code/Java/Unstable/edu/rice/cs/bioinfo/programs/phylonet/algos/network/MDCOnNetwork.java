package edu.rice.cs.bioinfo.programs.phylonet.algos.network;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 5/9/12
 * Time: 1:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class MDCOnNetwork {
    boolean _printDetail = false;

    public void setPrintDetails(boolean p){
        _printDetail = p;
    }
    /**
     * The public function for calculating the probabilities.
     * @param	net 	the given network
     * @param 	gts		the given set of gene trees
     * @param	alleles2species		the mapping from the names of allels to the names of the species. It is used for multiple alleles
     * @return	a list of probabilities corresponding to the list of gene trees.
     */
    public List<Integer> countExtraCoal(Network net, List<MutableTuple<Tree,Double>> gts, Map<String, String> alleles2species){
        List<Integer> xlList = new ArrayList<Integer>();
        Map<String,Integer> nname2tamount = new HashMap<String,Integer>();
        Tree superst = networkToTree(net, nname2tamount);
        //System.out.println(superst.toNewickWD());
        Map<String,String> tname2nname = new HashMap<String,String>();

        for(Map.Entry<String, Integer> entry: nname2tamount.entrySet()){
            if(entry.getValue()>1){
                for(int i=1; i<=entry.getValue(); i++){
                    tname2nname.put(entry.getKey()+"_"+i, entry.getKey());
                }
            }
        }

        for(MutableTuple<Tree,Double> tuple: gts){
            Tree gt = tuple.Item1;
            List<String> gtTaxa = Arrays.asList(gt.getLeaves());
            if(alleles2species == null){
                alleles2species = new HashMap<String,String>();
                for(String taxon: gtTaxa){
                    alleles2species.put(taxon, taxon);
                }
            }

            List<String> netTaxa = new ArrayList<String>();
            List<List<String>> allelesList = new ArrayList<List<String>>();
            //List<Integer> alleleNum = new ArrayList<Integer>();
            List<Integer> upper = new ArrayList<Integer>();
            for(String gtleaf: gtTaxa){
                String nleaf = alleles2species.get(gtleaf);
                int index = netTaxa.indexOf(nleaf);
                if(index==-1){
                    netTaxa.add(nleaf);
                    List<String> alleles = new ArrayList<String>();
                    alleles.add(gtleaf);
                    allelesList.add(alleles);
                    upper.add(nname2tamount.get(nleaf));
                }
                else{
                    allelesList.get(index).add(gtleaf);
                }
            }

            List<int[]> mergeNumber = new ArrayList<int[]>();
            for(List<String> alleles: allelesList){
                int[] first = new int[alleles.size()];
                Arrays.fill(first, 1);
                mergeNumber.add(first);
            }
            int minCoal = Integer.MAX_VALUE;
            do{
                Map<String, String> mapping = new HashMap<String,String>();
                for(int i=0; i<netTaxa.size(); i++){
                    String baseName = netTaxa.get(i);
                    List<String> alleles = allelesList.get(i);
                    int[] subscribes = mergeNumber.get(i);
                    for(int j=0; j<alleles.size(); j++){
                        mapping.put(alleles.get(j), baseName+"_"+subscribes[j]);
                    }
                }
                int coal = countExtraCoal(gt, superst, mapping, tname2nname);
                if(_printDetail){
                    System.out.println(mapping + ":  " + coal);
                }
                minCoal = Math.min(coal, minCoal);

            }while(mergeNumberAddOne(mergeNumber,upper));

            //System.out.println(gt +": " + min_coal + " extra lineages");
            xlList.add(minCoal);
        }
        return xlList;
    }

    private boolean mergeNumberAddOne(List<int[]> mergeNumber, List<Integer> upper){
        for(int i=0; i<mergeNumber.size(); i++){
            int[] partNumber = mergeNumber.get(i);
            int max = upper.get(i);
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

    private Tree networkToTree(Network<Double> net, Map<String,Integer> nname2tamount){
        removeBinaryNodes(net);
        Tree st = new STITree<Double>();
        Queue<NetNode<Double>> source = new LinkedList<NetNode<Double>>();
        Queue<TMutableNode> dest = new LinkedList<TMutableNode>();
        source.offer(net.getRoot());
        dest.offer((TMutableNode)st.getRoot());
        long nameid = System.currentTimeMillis();
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
                Integer amount = nname2tamount.get(name);
                if(amount==null){
                    amount = 0;
                }
                nname2tamount.put(name, ++amount);
                String newname = name + "_" + amount;
                copy = peer.createChild(newname);

                double distance = child.getParentDistance(parent);
                if (distance == NetNode.NO_DISTANCE) {
                    //copy.setParentDistance(TNode.NO_DISTANCE);
                    copy.setParentDistance(0);
                }
                else {
                    copy.setParentDistance(distance);
                }

                source.offer(child);
                dest.offer(copy);
            }
        }
        //Trees.removeBinaryNodes((MutableTree)st);
        return st;
    }

    private int countExtraCoal(Tree gt,Tree st, Map<String, String> taxonMap, Map<String,String> tname2nname){
        int sum = 0;
        String[] stTaxa = st.getLeaves();

        Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
        Map<String, Integer> nname2xl = new HashMap<String, Integer>();

        for (TNode node : st.postTraverse()) {
            BitSet bs = new BitSet();
            if (node.isLeaf()) {
                for (int i = 0; i < stTaxa.length; i++) {
                    if (node.getName().equals(stTaxa[i])) {
                        bs.set(i);
                        break;
                    }
                }
                map.put(node, bs);
            }
            else {
                for (TNode child : node.getChildren()) {
                    BitSet childCluster = map.get(child);
                    bs.or(childCluster);
                }
                map.put(node, bs);
            }

            if(node.getChildCount()==1 && node.getParentDistance()==0){
                continue;
            }

            STITreeCluster c = new STITreeCluster(stTaxa);
            c.setCluster(bs);

            if(c.getClusterSize()<stTaxa.length){
                int el = getClusterCoalNum(gt, c, taxonMap);
                //System.out.println(c + ":" + el);
                String tname = node.getName();
                if(tname!=null && tname2nname.containsKey(tname)){
                    String nname = tname2nname.get(tname);
                    Integer xl = nname2xl.get(nname);
                    if(xl==null){
                        xl = 0;
                    }
                    xl += el;
                    nname2xl.put(nname, xl);
                }
                else{
                    sum += Math.max(0, el-1);
                }
            }
        }

        for(Map.Entry<String, Integer> entry: nname2xl.entrySet()){
            sum += Math.max(0, entry.getValue()-1);
        }

        return sum;
    }

    private int getClusterCoalNum(Tree tr, STITreeCluster cluster, Map<String, String> taxonMap) {
        Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
        List<String> taxa = new LinkedList<String>();	// List of species taxa.

        Collections.addAll(taxa, cluster.getTaxa());

        int count = 0;
        for (TNode node : tr.postTraverse()) {
            if (node.isLeaf()) {
                String stTaxon = taxonMap.get(node.getName());	// Get the corresponding species name.
                int index = taxa.indexOf(stTaxon);
                BitSet bs = new BitSet(taxa.size());
                bs.set(index);
                if (cluster.containsCluster(bs)) {
                    count++;
                }
                map.put(node, bs);
            }
            else {
                BitSet bs = new BitSet(taxa.size());
                int intersect = 0;
                int childCount = node.getChildCount();
                for (TNode child : node.getChildren()) {
                    BitSet v = map.get(child);
                    bs.or(v);
                    if(childCount>2){
                        if(cluster.containsCluster(v)){
                            intersect ++;
                        }
                    }
                }

                if (cluster.containsCluster(bs)) {
                    count -= node.getChildCount();
                    count++;
                }
                else if(intersect>1){
                    count -= intersect;
                    count ++;
                }

                map.put(node, bs);
            }
        }
        return count;
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
