package edu.rice.cs.bioinfo.programs.phylonet.structs.network.util;


import edu.rice.cs.bioinfo.programs.phylonet.algos.bipartitematching.HungarianBipartiteMatcher;

/**
 * This class provides a simple way to represent the bipartite graph. The class Networks
 * use this class to compute the tree-based distance between two networks.
 */
public class BipartiteGraph {
    /**
     * This constructor initializes a bipartite graph with <code>lsize</code> nodes on
     * the left and <code>rsize</code> nodes on the right.
     */
    public BipartiteGraph(int lsize, int rsize)
    {
        assert(lsize > 0 && rsize > 0);

        _lsize = lsize;
        _rsize = rsize;

        _weights = new double[lsize][rsize];
        for (int l = 0; l < lsize; l++) {
            for (int r = 0; r < rsize; r++) {
                _weights[l][r] = NO_EDGE;
            }
        }
    }

    /**
     * This function adds an edge to the bipartite graph.
     *
     * @param l, r: indexes of the left and right endpoints for this edge.
     * @param w: the weight of this edge. We require it to be positive.
     */
    public void addEdge(int l, int r, double w)
    {
        assert(0 <= l && l < _lsize);
        assert(0 <= r && r < _rsize);
        assert(w >= 0);

        _weights[l][r] = w;
    }

    /**
     * This function returns the weight of a minimum edge cover of this bipartite graph.
     *
     * @return: The weight of a minimum edge cover.
     */
    public double getMinEdgeCoverWeight()
    {
        computeMinEdgeCover();

        double weight = _hbm.getMatchingWeight();
        weight = _hbm.getMatching().length * _invert_weight - weight;

        return weight / 2.0;
    }

    /**
     * Computes and returns the size of a minimum edge cover of this bipartite graph.
     *
     * @return: The number of edges in a minimum edge cover.
     */
    public int getMinEdgeCoverSize()
    {
        computeMinEdgeCover();

        int matching[][] = _hbm.getMatching();
        int size = matching.length;

        for (int i = 0; i < matching.length; i++) {
            int l = matching[i][0];
            int r = matching[i][1];

            if (l >= _lsize && r < _lsize) {
                size--;	// Edges in the copy of the graph is removed from the cover.
            }
        }

        return size;
    }

    /**
     * This function initializes the bipartite matcher to find the minimum-weight edge
     * cover. Refer to Chapters 17 and 19 in the book "Combinatorial Optimzation" by A. Schrijver
     * to see how the Hungarian algorithm is used to find the minimum edge cover.
     */
    private void initBipartiteMatcher()
    {
        _invert_weight = 2 * getMaxEdgeWeight() + 1.0;

        // Initialize the original and disjoint-copy parts for the new bipartite graph.
        _hbm = new HungarianBipartiteMatcher(_lsize + _rsize, _lsize + _rsize);

        for (int l = 0; l < _lsize; l++) {
            for (int r = 0; r < _rsize; r++) {
                if (_weights[l][r] != NO_EDGE) {
                    _hbm.addEdge(l, r + _lsize, _invert_weight - _weights[l][r]);	// Original part.
                    _hbm.addEdge(r + _lsize, l, _invert_weight - _weights[l][r]);	// Disjoint-copy part.
                }
            }
        }

        // Add edges for corresponding nodes.
        for (int l = 0; l < _lsize; l++) {
            double min = getMinIncidentEdgeWeight(l, true);

            _hbm.addEdge(l, l, _invert_weight -  2 * min);
        }
        for (int r = 0; r < _rsize; r++) {
            double min = getMinIncidentEdgeWeight(r, false);

            _hbm.addEdge(r + _lsize, r + _lsize, _invert_weight -  2 * min);
        }
    }

    /**
     * Computes the minimum edge covers for this bipartite graph.
     */
    private void computeMinEdgeCover()
    {
        initBipartiteMatcher();
        _hbm.findMatching();
    }

    /**
     * This function returns the minimum weight of all edges incident with a node indexed by
     * <code>index</index>.
     *
     * @param index: The index of the node, either in the left side or the right side.
     * @param left: <code>true</code> if this node is on the left side; <code>false</code> otherwise.
     */
    private double getMinIncidentEdgeWeight(int index, boolean left)
    {
        double min = Double.POSITIVE_INFINITY;

        if (left) {	// Find the min weight for a node on the left side.
            for (int r = 0; r < _rsize; r++) {
                if (min > _weights[index][r] && _weights[index][r] != NO_EDGE) {
                    min = _weights[index][r];
                }
            }
        }
        else {	// Find the min weight for a node on the right side.
            for (int l = 0; l < _lsize; l++) {
                if (min > _weights[l][index] && _weights[l][index] != NO_EDGE) {
                    min = _weights[l][index];
                }
            }
        }

        return min;
    }

    /**
     * This function returns the index of the edge with minimum weight of all edges incident
     * with a node indexed by <code>index</index>.
     *
     * @param index: The index of the node, either in the left side or the right side.
     * @param left: <code>true</code> if this node is on the left side; <code>false</code> otherwise.
     */
    private int getMinIncidentEdgeIndex(int index, boolean left)
    {
        double min = Double.POSITIVE_INFINITY;
        int i = -1;	// Indicating that this node has no incident edges.

        if (left) {	// Find the min weight for a node on the left side.
            for (int r = 0; r < _rsize; r++) {
                if (min > _weights[index][r] && _weights[index][r] != NO_EDGE) {
                    min = _weights[index][r];
                    i = r;
                }
            }
        }
        else {	// Find the min weight for a node on the right side.
            for (int l = 0; l < _lsize; l++) {
                if (min > _weights[l][index] && _weights[l][index] != NO_EDGE) {
                    min = _weights[l][index];
                    i = l;
                }
            }
        }

        return i;
    }

    /**
     * This function returns the maximum weight of all edges in the bipartite graph.
     */
    private double getMaxEdgeWeight()
    {
        double max = Double.NEGATIVE_INFINITY;

        for (int l = 0; l < _lsize; l++) {
            for (int r = 0; r < _rsize; r++) {
                if (max < _weights[l][r] && _weights[l][r] != NO_EDGE) {
                    max = _weights[l][r];
                }
            }
        }

        return max;
    }




    // Data members
    public final static double NO_EDGE = Double.NEGATIVE_INFINITY;

    private int _lsize;		// The number of nodes on the left side.
    private int _rsize;		// The number of nodes on the right side.
    private double _weights[][];			// Holds the edge weights.

    private HungarianBipartiteMatcher _hbm;	// The Hungarian matcher used to find the max. matching.
    private double _invert_weight;			// Hold a "big enough" weight, allowing the matcher to find min. matching.
}
