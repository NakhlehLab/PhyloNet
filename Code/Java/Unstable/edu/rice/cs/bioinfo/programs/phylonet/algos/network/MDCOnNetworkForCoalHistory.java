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
public class MDCOnNetworkForCoalHistory {
    boolean _printDetail = false;

    public void setPrintDetails(boolean p){
        _printDetail = p;
    }




    public int countExtraCoal(Tree mulTree, Tree gt, Map<String, String> alleleMapping, Map<String,String> tname2nname, Map<TNode, Integer> coalHistoryMapping, boolean[][] mulTreeMatrix){
        int sum = 0;
        String[] stTaxa = mulTree.getLeaves();

        Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
        Map<String, Integer> nname2xl = new HashMap<String, Integer>();

        for (TNode node : mulTree.postTraverse()) {
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


            STITreeCluster c = new STITreeCluster(stTaxa);
            c.setCluster(bs);

            int el = getClusterCoalNum(gt, c, node.getID(), alleleMapping, coalHistoryMapping, mulTreeMatrix);
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
                //System.out.println(node.toString() + ":" + Math.max(0, el-1));
            }


        }

        for(Map.Entry<String, Integer> entry: nname2xl.entrySet()){
            sum += Math.max(0, entry.getValue()-1);
            //System.out.println(entry.getKey() + ":" + Math.max(0, entry.getValue()-1));
        }

        return sum;
    }

    private int getClusterCoalNum(Tree tr, STITreeCluster cluster, int mulTreeNodeID, Map<String, String> alleleMapping, Map<TNode, Integer> coalHistoryMapping, boolean[][] mulTreeMatrix) {
        Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
        List<String> taxa = new LinkedList<String>();	// List of species taxa.

        Collections.addAll(taxa, cluster.getTaxa());

        int count = 0;
        for (TNode node : tr.postTraverse()) {
            if (node.isLeaf()) {
                String stTaxon = alleleMapping.get(node.getName());	// Get the corresponding species name.
                int index = taxa.indexOf(stTaxon);
                BitSet bs = new BitSet(taxa.size());
                bs.set(index);
                if (cluster.containsCluster(bs)) {
                    count++;
                }
                map.put(node, bs);
            }
            else {
                int coalNodeID = coalHistoryMapping.get(node);
                BitSet bs = new BitSet(taxa.size());
                //int intersect = 0;
                //int childCount = node.getChildCount();
                for (TNode child : node.getChildren()) {
                    BitSet v = map.get(child);
                    bs.or(v);
                    /*
                    if(childCount>2){
                        if(cluster.containsCluster(v)){
                            intersect ++;
                        }
                    }
                    */
                }

                if (cluster.containsCluster(bs) && mulTreeMatrix[mulTreeNodeID][coalNodeID]) {
                    count -= node.getChildCount();
                    count++;
                }
                /*
                else if(intersect>1){
                    count -= intersect;
                    count ++;
                }
                */

                map.put(node, bs);
            }
        }
        return count;
    }



}
