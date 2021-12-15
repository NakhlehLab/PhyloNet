package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
    private double _essPost;

    private long _lastSampleTime;

    private int examed = 0;

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
        this._lastSampleTime = System.currentTimeMillis();
    }

    /**
     * Main GTT loop.
     */
    public void run(int iteration, boolean doSample) {
        boolean accept;
        double logHastings;
        String op;

        for (int i = 0; i < _sampleFrequency; i++) {
            accept = false;
            _state.getUltrametricNetworkObject().checkLogDensity();
            logHastings = _state.propose();
            examed++;
            op = _state.getOperation().getName();

            if(Utils.DEBUG_MODE) {
                System.out.println(_temperature + " Proposed:"+ op + "\n" + _state.toString() + "\n" + _state.getNetwork() + "\n" + Networks.getTopologyString(_state.getNetworkObject()));
                System.out.println("logHasting: " + logHastings);

                Map<Tree, Integer> count = new HashMap<>();
                for(UltrametricTree gt : _state._geneTrees) {
                    boolean exist = false;
                    for(Tree key : count.keySet()) {

                        if(Trees.haveSameRootedTopology(key, gt.getTree())) {
                            exist = true;
                            count.put(key, count.get(key) + 1);
                            break;
                        }
                    }

                    if(!exist) {
                        count.put(gt.getTree(), 1);
                    }
                }

                for(Tree key : count.keySet()) {
                    System.out.println(key + " " + count.get(key));
                }
            }

            if(logHastings != Utils.INVALID_MOVE && Utils.SAMPLE_EMBEDDINGS) {
                // experimental!
                double embeddingLogHR = _state.getUltrametricNetworkObject().rebuildEmbeddings();

                if(embeddingLogHR == Utils.INVALID_MOVE) {
                    logHastings = Utils.INVALID_MOVE;
                } else {
                    logHastings += embeddingLogHR;
                }
            } else if(logHastings != Utils.INVALID_MOVE) {
            }

            if(logHastings != Utils.INVALID_MOVE) {

                if(Utils.DEBUG_MODE && !_state.isValidState()) {
                    System.out.println("INVALID state after operation and validation!!!"
                            + op + "\n" + _state.toString() + "\n" + _state.getNetwork());
                    throw new RuntimeException("INVALID state after operation and validation!!!"
                            + op + "\n" + _state.toString() + "\n" + _state.getNetwork());
                }

                double logPriorNext = _state.calculatePrior();
                double logLikelihoodNext = _state.calculateLikelihood();
                double logNext = logLikelihoodNext + logPriorNext;
                _logAlpha = (logNext - _logPost) / _temperature + logHastings;

                if(Utils.DEBUG_MODE) {
                    System.out.println("logPosteriorNext: " + logNext + " logPosteriorPrev: " + _logPost + " logAlpha: " + _logAlpha);
                }

                _state.getOperation().optimize(_logAlpha);

                if(op.equals("Slide-SubNet")) {
                    //System.out.println("Likelihood " + logNext + " " + _logPost);
                }

                if( _logAlpha >= Math.log(Randomizer.getRandomDouble()) ) {
                    if(op.equals("Slide-SubNet")) {
                       // System.out.println("!!!!! ");
                        //System.out.println(Networks.getDendroscopeCompatibleString(_state.getNetworkObject()));
                        //System.out.println(Networks.getFullString(_state.getNetworkObject()));
                    }
                    _logLikelihood = logLikelihoodNext;
                    _logPrior = logPriorNext;
                    _logPost = logNext;
                    accept = true;
                    _state.accept(_logAlpha);
                    if(Utils.DEBUG_MODE) {
                        System.out.println("Accepted");
                    }
                } else {
                    _state.undo(_logAlpha);
                    if (Utils.DEBUG_MODE) {
                        System.out.println("Rejected");
                    }
                }
            } else {
                _state.undo(Utils.INVALID_MOVE);
            }
            if(Utils.DEBUG_MODE && !_state.isValidState()) {
                System.out.println("INVALID state!!!"
                        + op + "\n" + _state.toString() + "\n" + _state.getNetwork());
                throw new RuntimeException("INVALID state!!!"
                        + op + "\n" + _state.toString() + "\n" + _state.getNetwork());
            }
            if(Utils.DEBUG_MODE) {
                System.out.println("iteration=" + iteration + " i=" + i);
                System.out.println();
            }

            if(_main) {
                _core.addInfo(accept, op);
            }
        }
        if(!Utils._PRE_BURN_IN && doSample) {
            if(_main) {
                if(Utils.DEBUG_MODE) {
                    System.out.println("Doing sampling");
                }
                _essPost = _core.addPosteriorESS(_logPost);
                double essPost = _core.addPosteriorESS(_logPost);
                double essPrior = _core.addPriorESS(_logPrior);

                long curTime = System.currentTimeMillis();
                System.out.printf("%d;    %2.5f;    %2.5f;    %2.5f;   %2.5f;    %2.5f;    %d;     %2.5f sec/sample;\n",
                        iteration, _logPost, essPost,
//                System.out.printf("%d;    %2.5f;    %2.5f;    %2.5f;   %2.5f;    %2.5f;    %d;\n",
//                        iteration, _logPost, _essPost,
                        _logLikelihood, _logPrior, essPrior,
                        _state.numOfReticulation(),
                        (curTime - _lastSampleTime) / 1000.0);
                _lastSampleTime = curTime;
                System.out.println(Networks.getFullString(_state.getNetworkObject()));
                _core.addSample(_state.toList());
                List<Tuple<String,Double>> netSample = new ArrayList<>();
                netSample.add(new Tuple<>(_state.getNetwork(), _logPost));
                _core.addNetSample(netSample);
                if(Utils.DEBUG_MODE) {
                    System.out.println("Sample added: " + _state.getNetwork().toString());
                }
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

    protected void setPreBurnInParams() {
        _state.setPreBurnInParams();
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

    public String getCurrentNetwork() { return _state.getNetworkObject().toString(); }

    public State getState() { return _state; }

    public double getLikelihood() { return _logLikelihood; }

    public double getPrior() { return _logPrior; }

    public double getESS() {
        return _essPost;
    }
}
