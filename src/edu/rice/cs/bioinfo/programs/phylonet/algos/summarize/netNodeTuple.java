package edu.rice.cs.bioinfo.programs.phylonet.algos.summarize;
/*
 *@ClassName: netNodeTuple
 *@Description: This class is the information need by tripartition2
 *@Author: Zhen Cao
 *@Date:  2019-07-30 11:50
 *@Version: 1.0
 */

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

import java.util.List;
import java.util.Set;

public class netNodeTuple{
    NetNode _node;
    Set<Set<String>> _siblingSet;
    Set<String> _retiLeaves;

    public netNodeTuple(NetNode node, Set<Set<String>> siblingSet, Set<String> retiLeaves){
        _node = node;
        _siblingSet = siblingSet;
        _retiLeaves = retiLeaves;
    }

    public NetNode getNode(){
        return _node;
    }

    public Set<String> getRetiLeaves(){
        return _retiLeaves;
    }

    public Set<Set<String>> getSiblingSet(){
        return _siblingSet;
    }

    public boolean equals(netNodeTuple t2){
        boolean equal = false;

        if (_retiLeaves.equals(t2.getRetiLeaves()) && _siblingSet.equals(t2.getSiblingSet())) {
            return true;

        }

        return false;
    }

    public String toString(){
        String s = _node.getName();
        s += ":(";
        s += String.join(",",_retiLeaves);
        s += "):";

        for(Set sl: _siblingSet){
            s += "(";
            s += String.join(",", sl);
            s += ")|";

        }

        s = s.substring(0, s.length()-1);
        s += ";";
        return s;
    }

}
