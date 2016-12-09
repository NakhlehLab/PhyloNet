package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.distance.Distance;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/**
 * Build UPGMA tree given distances
 * Created by wendingqiao on 1/26/16.
 */
public class UPGMATree {

    private Tree _upgmaTree;
    private Distance _distance;
    private Cluster _cluster;
    private int id = 0;

    class Cluster implements Comparable<Cluster>{
        public Cluster _c1;
        public Cluster _c2;
        public double _height;
        public Set<String> _taxa = new HashSet<>();

        // constructor for leaves
        public Cluster(String taxon) {
            _c1 = null;
            _c2 = null;
            _height = Utils.DEFAULT_TREE_LEAF_HEIGHT;
            _taxa.add(taxon);
        }
        // constructor for internal nodes
        public Cluster(Cluster c1, Cluster c2) {
            _c1 = c1;
            _c2 = c2;
            _height = _c1.getDistance(_c2);
            _taxa.addAll(_c1._taxa);
            _taxa.addAll(_c2._taxa);
        }
        @Override
        public int compareTo(Cluster o) {
            return _height > o._height ? 1 : -1;
        }
        public void resetCluster1(Cluster c1) {
            _c1 = c1;
            _height = _c1.getDistance(_c2);
        }
        public void resetCluster2(Cluster c2) {
            _c2 = c2;
            _height = _c1.getDistance(_c2);
        }
        private double getDistance(Cluster c) {
            double dist = 0.0;
            for(String taxon1 : _taxa) {
                for(String taxon2 : c._taxa) {
                    dist += _distance.getDistance(taxon1, taxon2);
                }
            }
            return dist / _taxa.size() / c._taxa.size();
        }
    }

    public UPGMATree(Distance distance) {
        this._distance = distance;
        Set<String> taxa = _distance.getTaxa();

        // build clusters list of leaves
        List<Cluster> clusters = new ArrayList<>();
        for(String taxon : taxa) {
            clusters.add(new Cluster(taxon));
        }
        int size = clusters.size();
        // build cluster pair priority queue
        PriorityQueue<Cluster> pq = new PriorityQueue<>(size * (size-1) / 2);
        for(int i = 0; i < clusters.size(); i++) {
            for (int j = i + 1; j < clusters.size(); j++) {
                pq.add(new Cluster(clusters.get(i), clusters.get(j)));
            }
        }
        // find the Cluster pair with minimum distance
        while(pq.size() > 1) {
            size--;
            Cluster cp = pq.poll();
            PriorityQueue<Cluster> temp = new PriorityQueue<>(size * (size-1) / 2);
            while(!pq.isEmpty()) {
                Cluster c = pq.poll();
                if(c._c1.equals(cp._c2) || c._c2.equals(cp._c2)) continue;
                if(c._c1.equals(cp._c1)) {
                    c.resetCluster1(cp);
                } else if (c._c2.equals(cp._c1)) {
                    c.resetCluster2(cp);
                }
                temp.add(c);
            }
            pq = temp;
        }
        _cluster = pq.poll();
        _upgmaTree = new STITree<>("I"+id++, true);
        buildTree(_cluster, (STINode) _upgmaTree.getRoot());
        setBranchLengths();
    }

    // This constructor should be used for debug only
    public UltrametricTree getUltrametricTree() {
        return new UltrametricTree(_upgmaTree);
    }

    public Tree getTree() {
        return _upgmaTree;
    }

    private double buildTree(Cluster cluster, STINode node) {
        if (cluster._c1 != null && cluster._c2 != null) {
            String name1 = cluster._c1._taxa.size() == 1 ? cluster._c1._taxa.iterator().next() : "I" + id++;
            STINode child1 = node.createChild(name1);
            String name2 = cluster._c2._taxa.size() == 1 ? cluster._c2._taxa.iterator().next() : "I" + id++;
            STINode child2 = node.createChild(name2);
            double childHeight = Math.max(buildTree(cluster._c1, child1), buildTree(cluster._c2, child2));
            if(cluster._height < childHeight) {
                cluster._height = childHeight * Utils.TREE_INTI_SCALE;
            }
        }
        node.setNodeHeight(cluster._height);
        return cluster._height;
    }

    private void setBranchLengths() {
        for(TNode node : _upgmaTree.getNodes()) {
            if(node.isRoot()) continue;
            double bl = node.getParent().getNodeHeight() - node.getNodeHeight();
            node.setParentDistance(bl);
        }
    }

}
