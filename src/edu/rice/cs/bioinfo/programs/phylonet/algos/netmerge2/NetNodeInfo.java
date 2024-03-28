package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge2;
/*
 * @ClassName:   NetNode
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        8/8/23 10:22 PM
 */


//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.structs.NetNodeInfo;

import java.util.HashMap;
import java.util.Map;

public class NetNodeInfo extends edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo{

    private double _height; // the time of the node
    private int _index; // the index of the node for likelihood calculation
//    private double _prevHeight = -1;
    private Map<String, Double> _NJDistances = null;
    private double _NJXSub = 0;

    public NetNodeInfo(double height) {
        super();
        this._height = height;
        this._index = -1;
        this._NJDistances = new HashMap<>();

    }


    public NetNodeInfo(){
        super();
        this._NJDistances = new HashMap<>();
    }


    public void setNJDistances(String leaf, double dist){
        _NJDistances.put(leaf, dist);

    }

    public void emptyNJDistances(){
        _NJDistances = new HashMap<>();

    }

    public void zeroNJXSub(){
        _NJXSub = 0;
    }

    public void incremetNJXSub(double d){
        _NJXSub += d;
    }

    public double getNJXsub(){
        return this._NJXSub;
    }

    public Map<String, Double> getNjDistances(){
        return this._NJDistances;
    }

    public void setNJXsub(double value){
        this._NJXSub = value;
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

//    public void storeHeight() {
//        this._prevHeight = this._height;
//    }
//
//    public void recoverHeight() {
//        this._height = this._prevHeight;
//        this._prevHeight = -1;
//    }

//    public double getPrevHeight() {
//        return this._prevHeight;
//    }

}
