package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.UndirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/25/12
 * Time: 2:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class UndirectedPhyloGraphDefault<N>  extends UndirectedGraphToGraphAdapter<N,PhyloEdge<N>> implements PhyloGraph<N>
{
    public UndirectedPhyloGraphDefault()
    {
        super( new UndirectedSparseGraph<N, PhyloEdge<N>>(),
                new Func1<PhyloEdge<N>, Tuple<N, N>>() {
                    public Tuple<N, N> execute(PhyloEdge<N> edge) {
                        return new Tuple<N, N>(edge.Source, edge.Destination);
                    }
                });
    }
}