package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.tree;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;

/**
 * Created by wendingqiao on 2/15/16.
 * Abstract operator for trees
 */
public abstract class TreeOperator extends Operator {

    protected UltrametricTree _tree;
    protected boolean _violate; // if the operation would violate the temporal constraint

    public TreeOperator(UltrametricTree tree)  {
        this._tree = tree;
        this._violate = false;
    }

    public Utils.MOVE_TYPE getCategory() {
        return Utils.MOVE_TYPE.TREE;
    }

    public boolean mayViolate() {
        return  _violate;
    }

    public void optimize (double logAlpha) {}

    protected TNode getOtherChild(TNode parent, TNode child) {
        for(TNode c : parent.getChildren()) {
            if(c.equals(child)) continue;
            return c;
        }
        throw new IllegalArgumentException("Cannot find the other child");
    }

    protected TNode getHigherChild(TNode parent) {
        TNode child = null;
        double height = -1;
        for(TNode c : parent.getChildren()) {
            if(_tree.getNodeHeight(c) > height) {
                child = c;
                height = _tree.getNodeHeight(c);
            }
        }
        return child;
    }

    protected void exchangeParents(TNode node1, TNode node2) {
        TNode p1 = node1.getParent();
        TNode p2 = node2.getParent();
        ((STINode) node1).setParent2((STINode) p2);
        ((STINode) node2).setParent2((STINode) p1);
        _tree.setNodeHeight(p1, _tree.getNodeHeight(p1));
        _tree.setNodeHeight(p2, _tree.getNodeHeight(p2));
    }

    protected TNode getChild(TNode parent, int idx) {
        int count = 0;
        for(TNode c : parent.getChildren()) {
            if(count++ == idx) return c;
        }
        return null;
    }

}
