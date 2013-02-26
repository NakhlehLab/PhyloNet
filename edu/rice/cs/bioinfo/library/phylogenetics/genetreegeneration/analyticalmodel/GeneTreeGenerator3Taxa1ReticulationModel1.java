package edu.rice.cs.bioinfo.library.phylogenetics.genetreegeneration.analyticalmodel;

import edu.rice.cs.bioinfo.library.phylogenetics.*;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.rice.cs.bioinfo.library.programming.Predicate2;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.math.BigInteger;
import java.util.Iterator;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/4/12
 * Time: 1:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneTreeGenerator3Taxa1ReticulationModel1<N,E,BL,HP> implements GeneTreeGenerator<N,E>
{
    private final Func1<E,BL> _getBranchLength;

    private final Predicate2<BL,BL> _areSameBranchLength;

    private final Predicate1<E> _isHalfProb;

    public GeneTreeGenerator3Taxa1ReticulationModel1(Func1<E, BL> getBranchLength, Predicate2<BL, BL> areSameBranchLength, Predicate1<E> isHalfProb) {
        _getBranchLength = getBranchLength;
        _areSameBranchLength = areSameBranchLength;
        _isHalfProb = isHalfProb;
    }

    public boolean canGenerateForNetwork(GraphReadOnly<N,E> network) {

        GetDirectSuccessors<N,E> getDirectSuccessors = new GetDirectSuccessors<N,E>();
        IsNetworkNode<N,E> isNetworkNode = new IsNetworkNode<N,E>();

        if(!network.isRooted())
        {
            return false;
        }

        if(IterableHelp.countInt(network.getNodes()) != 7)
        {
            return false;
        }

        Iterable leafs = new GetLeafs().execute(network);

        if(IterableHelp.countInt(leafs) != 3)
        {
            return false;
        }

        N root = new FindRoot<N>().execute(network);
        Iterable<N> rootChildren = getDirectSuccessors.execute(network, root);

        if(IterableHelp.countInt(rootChildren) != 2)
        {
            return false;
        }

        Iterator<N> rootChildenIt = rootChildren.iterator();
        N rootChild1 = rootChildenIt.next();
        N rootChild2 = rootChildenIt.next();

        N rootChild1HybridChild = null;
        int rootChild1ChildCount = 0;

        for(N child : getDirectSuccessors.execute(network, rootChild1))
        {
            rootChild1ChildCount++;
            if(isNetworkNode.execute(child, network))
            {
                rootChild1HybridChild = child;
            }
        }

        N rootChild2HybridChild = null;
        int rootChild2ChildCount = 0;

        for(N child : getDirectSuccessors.execute(network, rootChild2))
        {
            rootChild2ChildCount++;
            if(isNetworkNode.execute(child, network))
            {
                rootChild2HybridChild = child;
            }
        }

        if(rootChild1HybridChild == null || rootChild2HybridChild == null || rootChild1ChildCount != 2 || rootChild2ChildCount == 2 ||
           !rootChild1HybridChild.equals(rootChild2HybridChild))
        {
            return false;
        }

        N hybridNode = rootChild1HybridChild;

        E rootToRootChild1Branch = network.getEdge(root, rootChild1);
        E rootToRootChild2Branch = network.getEdge(root, rootChild2);

        if(!_areSameBranchLength.execute(_getBranchLength.execute(rootToRootChild1Branch),
                                         _getBranchLength.execute(rootToRootChild2Branch)))
        {
            for(N hybridParent : new GetDirectPredecessors<N,E>().execute(network, hybridNode))
            {
                E hybridParentEdge = network.getEdge(hybridParent, hybridNode);
                if(!_isHalfProb.execute(hybridParentEdge))
                {
                    return false;
                }
            }
        }

        return true;

    }

    public Iterable<GraphReadOnly> generateGeneTrees(GraphReadOnly network, BigInteger numGeneTrees) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }


}
