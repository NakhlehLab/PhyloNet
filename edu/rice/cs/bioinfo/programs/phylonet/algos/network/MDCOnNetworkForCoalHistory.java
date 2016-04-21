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
 * Created with Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is to count the number of extra lineages given a gene tree and its coalescent history within the branches of a multree
 *
 * See "Parsimonious Inference of Hybridization in the Presence of Incomplete Lineage Sorting", Systematic Biology, 2013.
 */
public class MDCOnNetworkForCoalHistory {
    boolean _printDetail = false;


    /**
     * This function is for setting the option of printing details
     */
    public void setPrintDetails(boolean p){
        _printDetail = p;
    }



    /**
     * This is the main function for counting the number of extra lineages given a gene tree and its coalescent history within the branches of a multree
     *
     * @param mulTree 	            the multree
     * @param gt		            the gene tree
     * @param alleleMapping	        the mapping from alleles to the species they are sampled from
     * @param tname2nname           the mapping from leaves in multree to leaves in original network
     * @param coalHistoryMapping    the coalescent history
     * @param mulTreeMatrix         a matrix that stores the ancestral relationships among nodes in multree
     *
     * @return	the number of extra lineages
     */
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

        for(Map.Entry<String, Integer> entry: nname2xl.entrySet()){
            sum += Math.max(0, entry.getValue()-1);
        }

        return sum;
    }



    /**
     * This function is to count the minimal number of extra lineages given a gene tree and a cluster in the multree
     *
     * @param tr 	                the gene tree
     * @param cluster		        the cluster in multree
     * @param mulTreeNodeID	        the node that the cluster is induced from
     * @param alleleMapping         the mapping from alleles to the species they are sampled from
     * @param coalHistoryMapping    the coalescent history
     * @param mulTreeMatrix         a matrix that stores the ancestral relationships among nodes in multree
     *
     * @return	the number of extra lineages
     */
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
                for (TNode child : node.getChildren()) {
                    BitSet v = map.get(child);
                    bs.or(v);
                }

                if (cluster.containsCluster(bs) && mulTreeMatrix[mulTreeNodeID][coalNodeID]) {
                    count -= node.getChildCount();
                    count++;
                }

                map.put(node, bs);
            }
        }
        return count;
    }



}
