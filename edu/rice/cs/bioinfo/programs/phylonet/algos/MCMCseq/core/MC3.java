package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;

/**
 * Created by wendingqiao on 10/17/14.
 */
public class MC3 {

    private MC3Core _core;

    private State _state;

    private double _temperature;
    private boolean _main;
    private int _sampleFrequency;

    private double _logAlpha;
    private double _logPost;
    private double _logLikelihood;
    private double _logPrior;

    public MC3(MC3Core core,
               State start,
               double temp,
               boolean main,
               int sampleFrequency
    ) {
        this._core = core;
        this._state = start;
        this._temperature = temp;
        this._main = main;
        this._sampleFrequency = sampleFrequency;
        this._logLikelihood = _state.calculateLikelihood();
        this._logPrior = _state.calculatePrior();
        this._logPost = this._logLikelihood + this._logPrior;
    }

    /**
     * Main GTT loop.
     */
    public void run(int iteration, boolean doSample) {

        if(iteration == 1 && doSample) {
            if(_main) {
                System.out.printf("%d;    %2.5f;    %2.5f;    %2.5f;   %2.5f;    %2.5f;    %d;\n",
                        0, _logPost, 0.0,
                        _logLikelihood, _logPrior, 0.0,
                        _state.numOfReticulation());
                System.out.println(_state.toString());
            }
        }

        boolean accept;
        double logHastings;
        String op;

        for (int i = 0; i < _sampleFrequency; i++) {
            accept = false;

            logHastings = _state.propose();
            op = _state.getOperation().getName();

            if(logHastings != Utils.INVALID_MOVE) {

                if(Utils.DEBUG_MODE && !_state.isValidState()) {
                    throw new RuntimeException("INVALID state after operation and validation!!!"
                            + op + "\n" + _state.toString() + "\n" + _state.getNetwork());
                }

                double logPriorNext = _state.calculatePrior();
                double logLikelihoodNext = _state.calculateLikelihood();
                double logNext = logLikelihoodNext + logPriorNext;
                _logAlpha = (logNext - _logPost) / _temperature + logHastings;

                _state.getOperation().optimize(_logAlpha);

                if( _logAlpha >= Math.log(Randomizer.getRandomDouble()) ) {
                    _logLikelihood = logLikelihoodNext;
                    _logPrior = logPriorNext;
                    _logPost = logNext;
                    accept = true;
                    _state.accept(_logAlpha);
                } else {
                    _state.undo(_logAlpha);
                }
            } else {
                _state.undo(Utils.INVALID_MOVE);
            }
            if(Utils.DEBUG_MODE && !_state.isValidState()) {
                throw new RuntimeException("INVALID state!!!"
                        + op + "\n" + _state.toString() + "\n" + _state.getNetwork());
            }

            if(_main) {
                _core.addInfo(accept, op);
            }
        }
        if(doSample) {
            if(_main) {
                double essPost = _core.addPosteriorESS(_logPost);
                double essPrior = _core.addPriorESS(_logPrior);
                System.out.printf("%d;    %2.5f;    %2.5f;    %2.5f;   %2.5f;    %2.5f;    %d;\n",
                        iteration, _logPost, essPost,
                        _logLikelihood, _logPrior, essPrior,
                        _state.numOfReticulation());
                _core.addSample(_state.toList());
                System.out.println(_state.toString());
            }
        }
    }

    /**
     * report state
     * @return    state info
     */
    public String report() {
        return this._temperature + " : " + _state.toString();
    }

    /********** getters and setters ***********/

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

    public double getPosterior() {
        return _logPost;
    }

}