package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.Random;

/**
 * Created by yunyu on 10/29/14.
 */
public abstract class EdgeParameterChange extends NetworkRearrangementOperation{
    double _windowSize = 0.1;

    public void setWindowSize(double windowSize){
        _windowSize = windowSize;
    }

    protected double drawRandomParameter(double originalValue, double lowerBound, double upperBound){
        double random;
        lowerBound = Math.max(lowerBound, originalValue - _windowSize);
        upperBound = Math.min(upperBound, originalValue + _windowSize);

        //System.out.println(originalValue + " in (" + lowerBound + "," + upperBound+")");
        do{
            random = originalValue + (Math.random()*2-1)*_windowSize;
        }while(random<lowerBound || random>upperBound);

        return random;
    }

    public void setParameters(Network network, Tuple<NetNode,NetNode> targetEdge, double targetEdgeBrlen, double targetEdgeInheriProb){
        super.setParameters(network, targetEdge, targetEdgeBrlen, targetEdgeInheriProb, null, null, null, null, null, null);
    }

}
