package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.summary.ESS;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.summary.Summary;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.util.concurrent.*;

/**
 * Created by wendingqiao on 12/8/14.
 */
public class MC3Core {

    // mcmc chains
    private List<MC3> _mc3s;

    // samples
    private List<List<Tuple<String,Double>>> _netSamples = new ArrayList<>();
    //private List<SampleSummary> _samples = new ArrayList<>();
    private boolean _samplingPhase = false;
    private int _burnInCounter = 0;
    private int _startNumber = 1;

    private List<Double> _likelihoodList = new ArrayList<>();
    private List<String> _networkList = new ArrayList<>();

    // logging info
    private ESS _posteriorESS = new ESS();
    private ESS _priorESS = new ESS();
    private boolean _sampling = false;
    private int _accept = 0;
    private Map<String, OperatorLogger> _opMap = new HashMap<>();
    ExecutorService executor = Executors.newFixedThreadPool(Utils._NUM_THREADS);

    public MC3Core(List<MarkerSeq> markerSeqs, BiAllelicGTR BAGTRModel, String logFileName) {
        if(logFileName != null) {
            System.out.println("Last log file: " + logFileName);
            int burnin = (int) (Utils._BURNIN_LEN / Utils._SAMPLE_FREQUENCY);
            System.out.println("----------------------- Previous logger: -----------------------");
            System.out.println("Iteration;    Posterior;  ESS;    Likelihood;    Prior;  ESS;    #Reticulation");
            try {
                BufferedReader br = new BufferedReader(new FileReader(logFileName));
                String line;
                String inputfile = "";
                String currentChain = "";
                boolean ready = false;

                while ((line = br.readLine()) != null) {
                    if (line.contains("/") && inputfile.equals("")) {
                        inputfile = line.substring(line.lastIndexOf('/') + 1);
                        System.out.println("Last Run Input: " + inputfile);
                    } else if (line.startsWith("Temp") || line.startsWith("(main)")) {
                        currentChain = line;
                        ready = true;
                    } else if (ready) {
                        Scanner lineScanner = new Scanner(line);
                        lineScanner.useDelimiter(";\\s*");
                        int numSamples = lineScanner.nextInt();
                        double logPosterior = lineScanner.nextDouble();
                        double posteriorESS = lineScanner.nextDouble();
                        double likelihood = lineScanner.nextDouble();
                        double logPrior = lineScanner.nextDouble();
                        double priorESS = lineScanner.nextDouble();
                        int numReticulation = lineScanner.nextInt();
                        if (currentChain.contains("main")) {
                            if (numSamples > burnin) {
                                _sampling = true;
                                _samplingPhase = true;
                            }
                            if(_sampling) {
                                posteriorESS = addPosteriorESS(logPosterior);
                                priorESS = addPriorESS(logPrior);
                            } else {
                                posteriorESS = 0.0;
                                priorESS = 0.0;
                            }
                            System.out.println("(main)");
                            System.out.printf("%d;    %2.5f;    %2.5f;    %2.5f;   %2.5f;    %2.5f;    %d;\n",
                                    numSamples, logPosterior, posteriorESS,
                                    likelihood, logPrior, priorESS,
                                    numReticulation);
                            String lastNetwork = br.readLine();
                            System.out.println(lastNetwork);
                            Utils._START_NET = lastNetwork;
                            if(_sampling) {
                                List<Tuple<String, Double>> netSample = new ArrayList<>();
                                netSample.add(new Tuple<>(lastNetwork, logPosterior));
                                addNetSample(netSample);
                            }
                            _startNumber = numSamples + 1;
                        }
                        ready = false;
                    } else if(line.startsWith("Gamma Mean: "))  {
                        if(currentChain.contains("main")) {
                            Utils._POP_SIZE_MEAN = 2.0 * Double.parseDouble(line.substring(12));
                            System.out.println(line);
                        }
                    } else if (line.startsWith("Rank = 0")) {
                        System.out.println("Last Run Finished ");
                    }
                }


            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        initializeChains(markerSeqs, BAGTRModel);
    }

    public MC3Core(List<MarkerSeq> markerSeqs, BiAllelicGTR BAGTRModel) {
        initializeChains(markerSeqs, BAGTRModel);
    }

    private void initializeChains(List<MarkerSeq> markerSeqs, BiAllelicGTR BAGTRModel) {
        try {
            int nChains = Utils._MC3_CHAINS == null ? 1 : 1 + Utils._MC3_CHAINS.size();
            this._mc3s = new ArrayList<>(nChains);
            for(int i = 0; i < nChains; i++) {
                _mc3s.add (new MC3(this,
                        new State(
                                Utils._START_NET,
                                Utils._START_GT_LIST,
                                markerSeqs,
                                Utils._POISSON_PARAM,
                                Utils._TAXON_MAP,
                                BAGTRModel
                        ),
                        i == 0 ? 1.0 : Utils._MC3_CHAINS.get(i - 1),
                        i == 0,
                        Utils.SWAP_FREQUENCY,
                        _networkList,
                        _likelihoodList
                ) );
            }
            _mc3s.get(0).setMain(true);

        } catch (Exception e) {
            e.printStackTrace();
        }

        //_samples.add(new SampleSummary("network", Utils.SampleType.Network));
        /*for(int i = 0; i < markerSeqs.size(); i++) {
            String name = markerSeqs.get(i).getName();
            if (name == null) name = Integer.toString(i);
            _samples.add(new SampleSummary("tree_" +  name, Utils.SampleType.Tree));
        }*/
        //if(Utils._ESTIMATE_POP_SIZE) {
        //    _samples.add(new SampleSummary("popSizePrior", Utils.SampleType.DoubleParam));
        //}
    }


    public String getBestLikelihoodNetwork() {
        int index = 0;
        double bestLikelihood = -1e99;
        for(int i = 0 ; i < _likelihoodList.size() ; i++) {
            if(bestLikelihood < _likelihoodList.get(i)) {
                bestLikelihood = _likelihoodList.get(i);
                index = i;
            }
        }
        return _networkList.get(index);
    }


    public void run() {
        System.out.println("----------------------- Logger: -----------------------");
        System.out.println("Iteration;    Posterior;  ESS;    Likelihood;    Prior;  ESS;    #Reticulation");
        int total = (int) (Utils._CHAIN_LEN / Utils._SAMPLE_FREQUENCY);
        int burnin = (int) (Utils._BURNIN_LEN / Utils._SAMPLE_FREQUENCY);
        int swap = (int) (Utils._SAMPLE_FREQUENCY / Utils.SWAP_FREQUENCY);
        for(int i = _startNumber; i <= total; i++) {
            if(i > burnin) {
                _sampling = true;
            }
            // swap state
            for(int k = 0; k < swap; k++) {
                runMC3s(i, k==0);
                if(_mc3s.size() > 1) doSwap();
            }
        }
        summarize();
    }

    private void doSwap() {
        int r1 = Randomizer.getRandomInt(_mc3s.size());
        int r2 = Randomizer.getRandomInt(_mc3s.size());
        while(r1 == r2) {
            r2 = Randomizer.getRandomInt(_mc3s.size());
        }
        MC3 mc1 = _mc3s.get(r1);
        MC3 mc2 = _mc3s.get(r2);
        double p1 = mc1.getPosterior(), p2 = mc2.getPosterior();
        double t1 = mc1.getTemperature(), t2 = mc2.getTemperature();
        boolean m1 = mc1.getMain(), m2 = mc2.getMain();
        double logRatio = (p2 - p1) * (1.0 / t1 - 1.0 / t2);
        boolean swap = (Math.log(Randomizer.getRandomDouble()) < logRatio);
        if(swap) {
            mc1.setMain(m2);
            mc2.setMain(m1);
            mc1.setTemperature(t2);
            mc2.setTemperature(t1);
        }
    }

    private void runMC3s(int iteration, boolean sample) {

        /*for(MC3 mc3 : _mc3s) {
            executor.execute(new Runnable() {
                public void run() {
                    mc3.run(iteration, sample);
                }
            });
        }
        //executor.shutdown();
        try {
            executor.awaitTermination(1000, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }*/

        for(MC3 mc3 : _mc3s) {
            mc3.run(iteration, sample);
        }
    }

    /**
     * Summarize the results (samples).
     */
    private void summarize() {
        System.out.println("----------------------- Summarization: -----------------------");
        System.out.printf("Burn-in = %d, " +
                        "Chain length = %d, " +
                        "Sample size = %d, " +
                        "Acceptance rate = %2.5f \n",
                Utils._BURNIN_LEN,
                Utils._CHAIN_LEN,
                (Utils._CHAIN_LEN - Utils._BURNIN_LEN) / Utils._SAMPLE_FREQUENCY,
                (double) (_accept) / (double) (Utils._CHAIN_LEN));

        System.out.println(getOperationDetails());

//        for(SampleSummary ss : _samples) {
//            ss.summary();
//            ss.close();
//        }

        Summary summary = new Summary(_netSamples, true);
        System.out.println("         -------------- Top Topologies: --------------");
        System.out.println(summary.getTopK(Utils._TOPK_NETS));
    }

    public void addSample(List<String> sample) {
        if(_samplingPhase) {
            //_samples.get(0).addSample(sample.get(0));
            //for(int i = 0; i < _samples.size(); i++) {
            //    _samples.get(i).addSample(sample.get(i));
            //}
            return;
        } else {
            if(_burnInCounter == Utils._BURNIN_LEN / Utils._SAMPLE_FREQUENCY) {
                _samplingPhase = true;
            }
            _burnInCounter++;
        }
    }

    public void addNetSample(List<Tuple<String,Double>> sample) {
        if(_samplingPhase || _sampling) {
            _netSamples.add(sample);
            return;
        }
    }

    public double addPosteriorESS(double ess) {
        return _sampling ? _posteriorESS.add(ess) : 0;
    }

    public double addPriorESS(double ess) {
        return _sampling ? _priorESS.add(ess) : 0;
    }

    public void addInfo(boolean ac, String op) {
        //System.out.println(op + " " + ac);
        if(ac) _accept++;
        if(!_opMap.containsKey(op)) _opMap.put(op, new OperatorLogger(op));
        OperatorLogger info = _opMap.get(op);
        info.count++;
        if(ac) info.accept++;
    }

    private String getOperationDetails() {
        StringBuilder sb = new StringBuilder("         --------------- Operations ---------------\n");
        for(String key : _opMap.keySet()) {
            sb.append(_opMap.get(key).toString());
            sb.append("\n");
        }
        return sb.toString();
    }

}
