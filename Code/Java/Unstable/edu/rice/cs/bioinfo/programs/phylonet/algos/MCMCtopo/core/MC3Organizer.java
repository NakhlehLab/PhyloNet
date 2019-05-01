package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.core;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary.ESS;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary.OperatorInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary.Summary;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.lang.reflect.Constructor;
import java.util.*;

/**
 * This is the metropolis-coupled MCMC chains organizer:
 * 1. swap MCMC chains samples
 * 2. collect samples from cold chain
 * 3. compute ESS (effective sample size) for posterior, likelihood and prior values
 *
 * Created by wendingqiao on 12/8/14.
 */
public class MC3Organizer {

    private long sampleSize = 100;

    private long _chainLength;

    private long _burnInLength;

    private int burnInIteration;

    private int accept = 0;

    private long _sampleFrequency;

    private MCMCMC[] mc3s;

    private int size;

    private Random _random;

    private List<List<Tuple<String,Double>>> _samples = new ArrayList<List<Tuple<String,Double>>>();

    private ESS posteriorESS = new ESS();

    private ESS priorESS = new ESS();

    private Map<String, OperatorInfo> opMap = new HashMap<String, OperatorInfo>();

    private List<String> logs = new ArrayList<String>();

    private Class _stateClass;

    private Class _mcmcClass;

    private List _trees;

    private double _k = 0.0;

    private long _seed;

    private int maxReti;


    public MC3Organizer(Class stateClass,
                        Class mcmcClass,
                        List<Network> start,
                        List trees,
                        Map<String, List<String>> taxonMap,
                        long chainLength,
                        long burnInLength,
                        long sampleFrequency,
                        long seed,
                        List<Double> temps,
                        double k,
                        int reti,
                        int parallel,
                        double[] weights
    ) {
        this._stateClass = stateClass;
        this._mcmcClass = mcmcClass;
        this._trees = trees;
        this._random = new Random(seed);
        this._seed = seed;
        this._k = k;
        this.maxReti = reti;

        _chainLength = chainLength;
        _burnInLength = burnInLength;
        _sampleFrequency = sampleFrequency;
        size = temps.size();
        burnInIteration = (int)(burnInLength / _sampleFrequency);
        try {
            mc3s = new MCMCMC[size];
            for(int i = 0; i < size; i++) {
                Network snet = start.size() == 1 ? start.get(0) :
                        (start.size() > size ? start.get(i) : null);
                initMC3(i, snet, taxonMap, temps, weights, parallel);
            }
            mc3s[0].setMain(true);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void initMC3(int idx, Network start, Map<String, List<String>> taxonMap,
                         List<Double> temps, double[] weights, int parallel) {
        try{
            Constructor stateClassConstructor = _stateClass.getConstructor(new Class[] {
                    Network.class, List.class, long.class, int.class, int.class, double[].class, Map.class
            });
            State state = (State) stateClassConstructor.newInstance(start, _trees, _seed, maxReti, parallel, weights, taxonMap);

            mc3s[idx] = (MCMCMC) _mcmcClass.getConstructor(new Class[] {
                    MC3Organizer.class, State.class, long.class, double.class, long.class
            }).newInstance(this, state, sampleSize, _k, _seed);

            mc3s[idx].setTemperature(temps.get(idx));
        } catch (Exception ex) {
            ex.getCause().printStackTrace();
            ex.printStackTrace();
        }
    }


    public void run() {
        System.out.println("Iteration;    Posterior;  ESS;    Likelihood;    Prior;  ESS;    #Reticulation");
        for(int i = 1; i <= _chainLength / _sampleFrequency; i++) {//burnInIteration; i++) {
            // swap state
            for(int k = 0; k < _sampleFrequency / sampleSize; k++) {
                runMC3s(i, k==0);
                if(mc3s.length > 1) doSwap();
            }
        }
        summarize();
    }

    private void doSwap() {
        int r1 = _random.nextInt(size);
        int r2 = _random.nextInt(size);
        while(r1 == r2) {
            r2 = _random.nextInt(size);
        }
        double p1 = mc3s[r1].getPosterior(), p2 = mc3s[r2].getPosterior();
        double t1 = mc3s[r1].getTemperature(), t2 = mc3s[r2].getTemperature();
        boolean m1 = mc3s[r1].getMain(), m2 = mc3s[r2].getMain();
        double logRatio = (p2 - p1) * (1.0 / t1 - 1.0 / t2);
        boolean swap = (Math.log(_random.nextDouble()) < logRatio);
        if(swap) {
            mc3s[r1].setMain(m2);
            mc3s[r2].setMain(m1);
            mc3s[r1].setTemperature(t2);
            mc3s[r2].setTemperature(t1);
        }
    }

    private void runMC3s(int iteration, boolean sample) {
        for(MCMCMC mc3 : mc3s) {
            mc3.run(iteration, sample);
        }
    }

    private int topK = 10;

    /**
     * Summarize the results (samples).
     */
    private void summarize() {
        System.out.println("----------------------- Summarization: -----------------------");
        System.out.printf("Burn-in = %d, " +
                        "Chain length = %d, " +
                        "Sample size = %d " +
                        "Acceptance rate = %2.5f \n",
                _burnInLength,
                _chainLength,
                (_chainLength-_burnInLength) / _sampleFrequency,
                (double) (accept) / (double) (_chainLength));

        System.out.println(getOperationDetails());

        Summary summary = new Summary(_samples, false);

        System.out.println("         -------------- 95% credible set of topologies: --------------");
        System.out.println(summary.getTopTopologies(95));
    }

    private int count = 0;
    private boolean add = false;

    public void addSample(List<Tuple<String,Double>> sample) {
        if(add) {
            _samples.add(sample);
            return;
        } else {
            count++;
            if(count == burnInIteration) {
                add = true;
            }
        }
    }

    public double addPosteriorESS(double ess) {
        return add ? posteriorESS.add(ess) : 0;
    }

    public double addPriorESS(double ess) {
        return add ? priorESS.add(ess) : 0;
    }

    public void incrementAC() {
        accept++;
    }

    public void addInfo(boolean ac, String op) {
        if(!opMap.containsKey(op)) opMap.put(op, new OperatorInfo(op));
        OperatorInfo info = opMap.get(op);
        info.count++;
        if(ac) info.accept++;
        opMap.put(op, info);
    }

    public void addLog(long i, boolean ac, double logPost) {
        logs.add(new String("Round-"+i + "; AC-"+ac + "; logPost-"+logPost));
    }

    private String getOperationDetails() {
        StringBuilder sb = new StringBuilder("         --------------- Operations ---------------\n");
        for(String key : opMap.keySet()) {
            sb.append(opMap.get(key).toString());
            sb.append("\n");
        }
        return sb.toString();
    }


}
