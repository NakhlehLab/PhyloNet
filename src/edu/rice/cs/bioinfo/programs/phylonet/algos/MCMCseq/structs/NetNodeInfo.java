package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs;

/**
 * Data structure that associates with each node of the network.
 * Created by wendingqiao on 2/25/16.
 */
public class NetNodeInfo {

    private double _height; // the time of the node
    private int _index; // the index of the node for likelihood calculation
    private double _prevHeight = -1;

    public NetNodeInfo(double height) {
        this._height = height;
        this._index = -1;
    }

    public NetNodeInfo(double height, int index, double prevHeight) {
        this._height = height;
        this._index = index;
        this._prevHeight = prevHeight;
    }

    public void setHeight(double newHeight) {
        this._height = newHeight;
    }

    public double getHeight() {
        return _height;
    }

    public void setIndex(int index) {
        this._index = index;
    }

    public int getIndex() {
        if(this._index == -1) {
            throw new IllegalArgumentException("Cannot access index value before initialization.");
        }
        return this._index;
    }

    public void storeHeight() {
        this._prevHeight = this._height;
    }

    public void recoverHeight() {
        this._height = this._prevHeight;
        this._prevHeight = -1;
    }

    public double getPrevHeight() {
        return this._prevHeight;
    }
}
