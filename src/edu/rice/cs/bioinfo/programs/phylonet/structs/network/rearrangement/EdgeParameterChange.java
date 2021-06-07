package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.Random;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of NetworkRearrangementOperation.
 * It handles the case of changing parameter of an edge
 */
public abstract class EdgeParameterChange extends NetworkRearrangementOperation{
    double _windowSize = 0.1;


    /**
     * This function is to set the size of the window when changing the parameters
     */
    public void setWindowSize(double windowSize){
        _windowSize = windowSize;
    }


    /**
     * This function is to randomly draw a parameter within the window which is centered at the original value
     *
     * @param originalValue  the original value
     * @param lowerBound     the lower bound
     * @param upperBound     the upper bound
     */
    protected double drawRandomParameter(double originalValue, double lowerBound, double upperBound){
        double random;
        lowerBound = Math.max(lowerBound, originalValue - _windowSize);
        upperBound = Math.min(upperBound, originalValue + _windowSize);
        do{
            random = originalValue + (Math.random()*2-1)*_windowSize;
        }while(random<lowerBound || random>upperBound);

        return random;
    }


    /**
     * This function is set the parameters
     *
     * @param network               the species network to be rearranged
     * @param targetEdge            the target edge that is involved in this rearrangement
     * @param targetEdgeBrlen       the branch length of <code>targetEdge</code>
     * @param targetEdgeInheriProb  the inheritance probability of <code>targetEdge</code>
     */
    public void setParameters(Network network, Tuple<NetNode,NetNode> targetEdge, double targetEdgeBrlen, double targetEdgeInheriProb){
        super.setParameters(network, targetEdge, targetEdgeBrlen, targetEdgeInheriProb, null, null, null, null, null, null);
    }

}
