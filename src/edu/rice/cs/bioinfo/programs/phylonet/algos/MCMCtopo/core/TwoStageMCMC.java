package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core;

/**
 * Created by dw20 on 5/12/17.
 */
public class TwoStageMCMC extends MCMCMC {

    protected double _logPseudoLikelihood;

    public TwoStageMCMC(MC3Organizer organizer, State start, long sampleFrequency, double k, long seed) {
        super(organizer, start, sampleFrequency, k, seed);
        _logPseudoLikelihood = ((TwoStageState) _state).calculatePseudoLikelihood();
    }

    @Override
    public boolean acceptNewState(double logThisPrior, double logQ) {
        return false;
    }

}
