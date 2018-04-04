package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.summary;

import edu.rice.cs.bioinfo.library.programming.Tuple;

/**
 * Created by wendingqiao on 3/28/16.
 */
public class OperatorLog {

    private Tuple<String,Double> MAP;
    private double meanValue;

    public OperatorLog(Tuple<String,Double> t) {
        this.MAP = t;
        this.meanValue += t.Item2;
    }

    public void addSample(Tuple<String,Double> t) {
        if(this.MAP.Item2 < t.Item2) this.MAP = t;
        this.meanValue += t.Item2;
    }

    public String toString(int num) {
        return this.MAP.Item2 + ":" + this.MAP.Item1 + " Ave=" + getMean(num);
    }

    private double getMean(int num) {
        return this.meanValue / (double)num;
    }
}