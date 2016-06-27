package edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling;

import java.util.List;

/**
 * Created by yunyu on 10/1/15.
 * D: sample
 * E: observations
 */
public abstract class  GibbsSamplingBase <H, E extends Object> {
    protected int _numIterations = 100000;
    protected int _burnIn = 100;
    protected int _sampleInterval = 100;
    protected List<H> sampledHypothesis;
    protected E _observations;

    protected GibbsSamplingBase(){};

    protected GibbsSamplingBase(int numIterations, int burnIn, int sampleInterval){
        _numIterations = numIterations;
        _burnIn = burnIn;
        _sampleInterval = sampleInterval;
    }

    protected abstract H generateInitialSample();

    protected abstract H getNextSample(H currentHypothesis);

    protected abstract void updateNextValueFromConditionalPosterior(H hypothesis, int index);

    protected abstract double computePrior(H hypothesis);

    protected abstract double computeLikelihood(H hypothesis);
}
