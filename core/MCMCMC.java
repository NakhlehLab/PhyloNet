package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * The class of core MCMCMC chain
 *
 * Created by wendingqiao on 10/17/14.
 */
public class MCMCMC {

    private MC3Organizer organizer;

    private long _sampleFrequency;

    private double k;

    private Random _random;

    private double temperature;

    private boolean main = false;

    public void setTemperature(double temp) {
        this.temperature = temp;
    }

    public void setMain(boolean m) {
        this.main = m;
    }

    public double getTemperature() {
        return this.temperature;
    }

    public boolean getMain() {
        return this.main;
    }


    /**
     * GTT Constructor use defined settings.
     */
    public MCMCMC(MC3Organizer organizer,
                  State start,
                  long sampleFrequency,
                  double k,
                  long seed
    ) {
        this.organizer = organizer;
        this.state = start;
        this._sampleFrequency = sampleFrequency;
        this.k = k;
        this._random = new Random(seed);

        this.logLikelihood = state.calculateLikelihood();
        this.logPrior = state.calculatePrior(k);
        this.logPost = logLikelihood + logPrior;
    }

    /**
     * Metropolis Hastings ratio.
     */
    private double logAlpha;

    /**
     *  State posterior value
     */
    private double logPost;

    /**
     * State Likelihood value
     */
    private double logLikelihood;

    /**
     * State Prior value;
     */
    private double logPrior;


    public double getPosterior() {
        return logPost;
    }

    /**
     * Current state.
     */
    private State state;

    /**
     * Main GTT loop.
     */
    public void run(int iteration, boolean doSample) {

        if(iteration == 1 && doSample) {
            if(main) {
                System.out.printf("%d;    %2.5f;    %2.5f;    %2.5f;   %2.5f;    %2.5f;    %d;\n",
                        0, logPost, 0.0,
                        logLikelihood, logPrior, 0.0,
                        state.numOfReticulation());
                System.out.println(state.toString());
            }
        }

        for (int i = 0; i < _sampleFrequency; i++) {
            String prev = state.getNetwork().toString();
            boolean ac = false;
            double logQ = state.propose();
            double logThisPrior;

            if(logQ != Double.MIN_VALUE && ((logThisPrior = state.calculatePrior(k)) != Double.MIN_VALUE)) {

                double logLikelihoodNext = state.calculateLikelihood();
                double logPriorNext = logThisPrior;
                double logNext = logLikelihoodNext + logPriorNext;
                logAlpha = (logNext - logPost) / temperature + logQ;

                if(logAlpha >= Math.log(_random.nextDouble())) {
                    logLikelihood = logLikelihoodNext;
                    logPrior = logPriorNext;
                    logPost = logNext;
                    ac = true;

                } else {
                    state.setNetwork(Networks.readNetwork(prev));
                }
            } else {
                // nullify the operation, or prior = 0
                // means the proposal or the network is not valid (e.g. has cycle)
                // undo may not be valid
                state.setNetwork(Networks.readNetwork(prev));
            }

            if(main) {
                if(ac) organizer.incrementAC();
                organizer.addInfo(ac, state.getOperation());
                organizer.addLog(i, ac, logPost);
                List<Double> list = new ArrayList<>();
                list.add(logPost);
                list.add(logLikelihood);
                list.add(logPrior);
            }
        }
        if(doSample) {
            if(main) {
                organizer.addSample(state.storeState(this.logPost));
                double essPost = organizer.addPosteriorESS(logPost);
                double essPrior = organizer.addPriorESS(logPrior);
                System.out.printf("%d;    %2.5f;    %2.5f;    %2.5f;   %2.5f;    %2.5f;    %d;\n",
                        iteration, logPost, essPost,
                        logLikelihood, logPrior, essPrior,
                        state.numOfReticulation());
                System.out.println(state.toString());
            }
        }
    }

    public String report() {
        return temperature + " : " + state.toString();
    }

}
