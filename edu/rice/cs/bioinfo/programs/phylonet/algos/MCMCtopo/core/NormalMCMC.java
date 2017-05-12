package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

/**
 * Created by dw20 on 5/12/17.
 */
public class NormalMCMC extends MCMCMC {

    public NormalMCMC(MC3Organizer organizer, State start, long sampleFrequency, double k, long seed) {
        super(organizer, start, sampleFrequency, k, seed);
    }

    @Override
    public boolean acceptNewState(double logThisPrior, double logQ) {
        double logLikelihoodNext = _state.calculateLikelihood();
        double logPriorNext = logThisPrior;
        double logNext = logLikelihoodNext + logPriorNext;
        _logAlpha = (logNext - _logPost) / _temperature + logQ;

        if(_logAlpha >= Math.log(_random.nextDouble())) {
            _logLikelihood = logLikelihoodNext;
            _logPrior = logPriorNext;
            _logPost = logNext;
            return true;
        }
        return false;
    }
}
