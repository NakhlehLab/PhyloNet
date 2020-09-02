package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs;

/**
 * Data structure that associates with each node of the model species tree.
 * Stores:
 * 1. population size of the branch incident into this node
 * 2. index of this node used for generating MS command
 *
 * Created by Xinhao Liu on 11/1/19.
 */
public class TreeNodeInfo {
    private int _popSize;
    private int _prevPopSize = -1; //TODO: what is this

    private int _index; // The index of this node for generating MS command
    private int _label; // Integer label of this node for generating coalescent history
    public static int NO_LABEL = Integer.MIN_VALUE;

    public TreeNodeInfo(int popSize) {
        this._popSize = popSize;
        this._label = NO_LABEL;
    }

    public void setPopSize(int newPopSize) {
        this._popSize = newPopSize;
    }

    public int getPopSize() {
        return _popSize;
    }

    //TODO: what is this
    public void storePopSize() {
        this._prevPopSize = this._popSize;
    }

    public void recoverPopSize() {
        this._popSize = this._prevPopSize;
        this._prevPopSize = -1;
    }

    public void setIndex(int newIndex) {
        this._index = newIndex;
    }

    public int getIndex() {
        return _index;
    }

    public void setLabel(int newLabel) {
        this._label = newLabel;
    }

    public int getLabel() {
        return _label;
    }

}
