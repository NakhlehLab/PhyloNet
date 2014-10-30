package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

/**
 * Created by yunyu on 10/29/14.
 */
public abstract class EdgeParameterChange extends NetworkRearrangementOperation{
    double _windowSize = 0.1;
    double _maxBrlen = 6;

    public void setWindowSize(double windowSize){
        _windowSize = windowSize;
    }

    protected double drawRandomParameter(double originalValue, double lowerBound, double upperBound){
        double random;
        lowerBound = Math.max(lowerBound, originalValue - _windowSize);
        upperBound = Math.min(upperBound, originalValue + _windowSize);
        do{
            random = originalValue + (Math.random()*2-1)*_windowSize;
        }while(random<lowerBound || random>upperBound);
        return random;
    }

}
