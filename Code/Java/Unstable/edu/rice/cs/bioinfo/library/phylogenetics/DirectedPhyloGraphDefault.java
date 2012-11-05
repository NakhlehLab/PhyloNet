package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.uci.ics.jung.graph.DirectedSparseGraph;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/25/12
 * Time: 2:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class DirectedPhyloGraphDefault<N>  extends DirectedGraphToGraphAdapter<N,PhyloEdge<N>> implements PhyloGraph<N>
{
    public DirectedPhyloGraphDefault()
    {
        super( new DirectedSparseGraph<N, PhyloEdge<N>>(), new Func1<PhyloEdge<N>, Tuple<N, N>>() {
            public Tuple<N, N> execute(PhyloEdge<N> edge) {
                return edge.NodesOfEdge;
            }
        });
    }
}