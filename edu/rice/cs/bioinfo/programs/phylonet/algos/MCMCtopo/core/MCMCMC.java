package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * The class of core MCMCMC chain
 *
 * Created by wendingqiao on 10/17/14.
 */
public abstract class MCMCMC {

    protected MC3Organizer _organizer;

    protected long _sampleFrequency;

    protected double _k;

    protected Random _random;

    protected double _temperature;

    protected boolean _main = false;

    public void setTemperature(double temp) {
        this._temperature = temp;
    }

    public void setMain(boolean m) {
        this._main = m;
    }

    public double getTemperature() {
        return this._temperature;
    }

    public boolean getMain() {
        return this._main;
    }

    /**
     * Metropolis Hastings ratio.
     */
    protected double _logAlpha;

    /**
     *  State posterior value
     */
    protected double _logPost;

    /**
     * State Likelihood value
     */
    protected double _logLikelihood;

    /**
     * State Prior value;
     */
    protected double _logPrior;


    public double getPosterior() {
        return _logPost;
    }

    /**
     * Current state.
     */
    protected State _state;


    /**
     * GTT Constructor use defined settings.
     */
    public MCMCMC(MC3Organizer organizer,
                  State start,
                  long sampleFrequency,
                  double k,
                  long seed
    ) {
        this._organizer = organizer;
        this._state = start;
        this._sampleFrequency = sampleFrequency;
        this._k = k;
        this._random = new Random(seed);

        this._logLikelihood = _state.calculateLikelihood();
        this._logPrior = _state.calculatePrior(k);
        this._logPost = _logLikelihood + _logPrior;
    }

    /**
     * Main GTT loop.
     */
    public void run(int iteration, boolean doSample) {
        initialSystemOut(iteration == 1 && doSample);
        for (int i = 0; i < _sampleFrequency; i++) {
            String prev = _state.getNetwork().toString();
            boolean ac = false;
            double logQ = _state.propose();
            double logThisPrior = Double.MIN_VALUE;
            boolean proceed2likelihood = (logQ != Double.MIN_VALUE) &&
                    ((logThisPrior = _state.calculatePrior(_k)) != Double.MIN_VALUE);

            if(!proceed2likelihood || !(ac = acceptNewState(logThisPrior, logQ))) {
                // nullify the operation, or prior = 0
                // means the proposal or the network is not valid (e.g. has cycle)
                // undo may not be valid
                _state.setNetwork(Networks.readNetwork(prev));
            }

            if(_main) {
                if(ac) _organizer.incrementAC();
                _organizer.addInfo(ac, _state.getOperation());
                _organizer.addLog(i, ac, _logPost);
                List<Double> list = new ArrayList<>();
                list.add(_logPost);
                list.add(_logLikelihood);
                list.add(_logPrior);
            }
        }
        sampleSystemOut(doSample, iteration);
    }

    public abstract boolean acceptNewState(double logThisPrior, double logQ);

    public String report() {
        return _temperature + " : " + _state.toString();
    }

    private void initialSystemOut(boolean print) {
        if(!print || !_main) return;
        System.out.printf("%d;\t%2.5f;\t%2.5f;\t%2.5f;\t%2.5f;\t%2.5f;\t%d;\n",
                0, _logPost, 0.0,
                _logLikelihood, _logPrior, 0.0,
                _state.numOfReticulation());
        System.out.println(_state.toString());
    }

    private void sampleSystemOut(boolean print, int iteration) {
        if(!print || !_main) return;
        _organizer.addSample(_state.storeState(this._logPost));
        double essPost = _organizer.addPosteriorESS(_logPost);
        double essPrior = _organizer.addPriorESS(_logPrior);
        System.out.printf("%d;\t%2.5f;\t%2.5f;\t%2.5f;\t%2.5f;\t%2.5f;\t%d;\n",
                iteration, _logPost, essPost,
                _logLikelihood, _logPrior, essPrior,
                _state.numOfReticulation());
        System.out.println(_state.toString());
    }

}
