package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 8/27/18
 * Time: 11:47 AM
 * To change this template use File | Settings | File Templates.
 */
public class TreeEmbedding {
    // edge (parent, child) -> list of clusters
    public Map<Tuple<NetNode, NetNode>, Collection<STITreeCluster>> embedding;
    public double prob = 1.0;
    public double probsum = 1.0;
    private String gtTaxa_[];
    private Map<TNode, STITreeCluster> clusterMap_;

    public TreeEmbedding() {
    }

    public TreeEmbedding(TreeEmbedding src) {
        setTo(src);
    }

    public TreeEmbedding(String gtTaxa[], Map<TNode, STITreeCluster> clusterMap) {
        gtTaxa_ = gtTaxa.clone();
        clusterMap_ = new HashMap<>(clusterMap);
        embedding = new HashMap<>();
    }

    @Override
    public TreeEmbedding clone() {
        TreeEmbedding newEmbedding = new TreeEmbedding(gtTaxa_, clusterMap_);
        newEmbedding.setTo(this);
        return newEmbedding;
    }

    public void setDirection(TNode gtnode, NetNode parent, NetNode child) {
        Tuple<NetNode, NetNode> edge = new Tuple(parent, child);
        if(!embedding.containsKey(edge)) {
            embedding.put(edge, new HashSet<>());
        }
        embedding.get(edge).add(clusterMap_.get(gtnode));
    }

    public void merge(TreeEmbedding src) {
        prob *= src.prob;
        probsum *= src.probsum;

        for(Tuple<NetNode, NetNode> edge : src.embedding.keySet()) {
            if(!embedding.containsKey(edge)) {
                embedding.put(edge, new HashSet<>());
            }
            embedding.get(edge).addAll(src.embedding.get(edge));
        }
    }

    public void setTo(TreeEmbedding src) {
        prob = src.prob;
        probsum = src.probsum;
        embedding = new HashMap<>();
        embedding.putAll(src.embedding);
        gtTaxa_ = src.gtTaxa_.clone();
        clusterMap_ = new HashMap<>(src.clusterMap_);
    }

    @Override
    public int hashCode() {
        int ret = Arrays.hashCode(this.gtTaxa_);
        //for(Tuple<NetNode, NetNode> tuple : this.embedding.keySet()) {
        //    ret = ret ^ tuple.hashCode();
        //}
        return ret;
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof TreeEmbedding))
            return false;
        if (obj == this)
            return true;

        TreeEmbedding rhs = (TreeEmbedding) obj;

        if(!Arrays.equals(this.gtTaxa_, rhs.gtTaxa_)) {
            return false;
        }

        for(Tuple<NetNode, NetNode> tuple : this.embedding.keySet()) {
            if(!rhs.embedding.containsKey(tuple)) return false;
            /*for(STITreeCluster cl : rhs.embedding.get(tuple)) {
                if(!this.embedding.get(tuple).contains(cl))
                    System.out.println(cl.hashCode());
                    System.out.println(this.embedding.get(tuple).iterator().next().hashCode());
                    return false;
            }*/
            if(!rhs.embedding.get(tuple).containsAll(this.embedding.get(tuple)))
                return false;
            if(!this.embedding.get(tuple).containsAll(rhs.embedding.get(tuple))) return false;
        }
        for(Tuple<NetNode, NetNode> tuple : rhs.embedding.keySet()) {
            if(!this.embedding.containsKey(tuple)) return false;
            if(!rhs.embedding.get(tuple).containsAll(this.embedding.get(tuple))) return false;
            if(!this.embedding.get(tuple).containsAll(rhs.embedding.get(tuple))) return false;
        }

        return true;

    }

    public boolean initialized() {
        if(gtTaxa_ == null) return false;
        return true;
    }

    public void clear() {
        gtTaxa_ = null;
    }

}
