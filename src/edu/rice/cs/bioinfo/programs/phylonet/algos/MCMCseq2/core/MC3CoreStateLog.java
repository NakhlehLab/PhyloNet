//package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.core;
///*
// * @ClassName:   MC3CoreStateLog
// * @Description:
// * @Author:      Zhen Cao
// * @Date:        3/13/23 2:39 PM
// */
//
//import edu.rice.cs.bioinfo.library.programming.Tuple;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.felsenstein.alignment.Alignment;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.structs.NetNodeInfo;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.structs.UltrametricNetwork;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.structs.UltrametricTree;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.summary.SampleSummary;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.util.Randomizer;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.util.Utils;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary.ESS;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary.Summary;
//import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
//import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
//import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
//
//import java.io.BufferedWriter;
//import java.io.FileWriter;
//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.HashMap;
//import java.util.List;
//import java.util.Map;
//
///**
// * Created by wendingqiao on 12/8/14.
// */
//public class MC3CoreStateLog extends MC3Core{
//
//    // mcmc chains
//    private List<MC3> _mc3s;
//
//    // samples
//    private List<List<Tuple<String,Double>>> _netSamples = new ArrayList<>();
//    private List<SampleSummary> _samples = new ArrayList<>();
//    private boolean _samplingPhase = false;
//    private int _burnInCounter = 0;
//
//    // logging info
//    private ESS _posteriorESS = new ESS();
//    private ESS _priorESS = new ESS();
//    private ESS _likelihoodESS = new ESS();
//    private boolean _sampling = false;
//    private int _accept = 0;
//    private Map<String, OperatorLogger> _opMap = new HashMap<>();
//    private int MaxReticulation = Utils._NET_MAX_RETI;
//
//
//    public MC3CoreStateLog(List<Alignment> alignments, String lastStateDir) {
//        super(alignments, lastStateDir);
//    }
//
//    public MC3CoreStateLog(List<Alignment> alignments) {
//        super(alignments);
//    }
//
//
//    private void doSwap() {
//        int r1 = Randomizer.getRandomInt(_mc3s.size());
//        int r2 = Randomizer.getRandomInt(_mc3s.size());
//        while(r1 == r2) {
//            r2 = Randomizer.getRandomInt(_mc3s.size());
//        }
//        MC3 mc1 = _mc3s.get(r1);
//        MC3 mc2 = _mc3s.get(r2);
//        double p1 = mc1.getPosterior(), p2 = mc2.getPosterior();
//        double t1 = mc1.getTemperature(), t2 = mc2.getTemperature();
//        boolean m1 = mc1.getMain(), m2 = mc2.getMain();
//        double logRatio = (p2 - p1) * (1.0 / t1 - 1.0 / t2);
//        boolean swap = (Math.log(Randomizer.getRandomDouble()) < logRatio);
//        if(swap) {
//            mc1.setMain(m2);
//            mc2.setMain(m1);
//            mc1.setTemperature(t2);
//            mc2.setTemperature(t1);
//        }
//    }
//
//    private void runMC3s(int iteration, boolean sample) {
//        for(MC3 mc3 : _mc3s) {
//            mc3.run(iteration, sample);
//        }
//    }
//
//    /**
//     * Summarize the results (samples).
//     */
//    private void summarize() {
//        System.out.println("----------------------- Summarization: -----------------------");
//        System.out.printf("Burn-in = %d, " +
//                        "Chain length = %d, " +
//                        "Sample size = %d, " +
//                        "Acceptance rate = %2.5f \n",
//                Utils._BURNIN_LEN,
//                Utils._CHAIN_LEN,
//                (Utils._CHAIN_LEN - Utils._BURNIN_LEN) / Utils._SAMPLE_FREQUENCY,
//                (double) (_accept) / (double) (Utils._CHAIN_LEN));
//
//        System.out.println(getOperationDetails());
//
//        for(SampleSummary ss : _samples) {
//            ss.summary();
//            ss.close();
//        }
//
//        Summary summary = new Summary(_netSamples, true);
//        System.out.println("         -------------- 95% credible set of topologies --------------");
//        System.out.println(summary.getTopTopologies(95));
//    }
//
//    public void addSample(List<String> sample) {
//        if(_samplingPhase) {
//            for(int i = 0; i < _samples.size(); i++) {
//                _samples.get(i).addSample(sample.get(i));
//            }
//            return;
//        } else {
//            if(_burnInCounter == Utils._BURNIN_LEN / Utils._SAMPLE_FREQUENCY) {
//                _samplingPhase = true;
//            }
//            _burnInCounter++;
//        }
//    }
//
//    public void addNetSample(List<Tuple<String,Double>> sample) {
//        if(_samplingPhase) {
//            _netSamples.add(sample);
//            return;
//        }
//    }
//
//    public double addPosteriorESS(double ess) {
//        return _sampling ? _posteriorESS.add(ess) : 0;
//    }
//
//    public double addLikelihoodESS(double ess){
//        return _sampling ? _likelihoodESS.add(ess) : 0;
//    }
//
//    public double addPriorESS(double ess) {
//        return _sampling ? _priorESS.add(ess) : 0;
//    }
//
//    public void addInfo(boolean ac, String op, double loghastings) {
//        if(ac) _accept++;
//        if(!_opMap.containsKey(op)) _opMap.put(op, new OperatorLogger(op));
//        OperatorLogger info = _opMap.get(op);
//        info.count++;
//        if(ac) info.accept++;
//        if (loghastings == Utils.INVALID_MOVE) {
//            info.rejectOP++;
//        }
//    }
//
//    private String getOperationDetails() {
//        StringBuilder sb = new StringBuilder("         --------------- Operations ---------------\n");
//        for(String key : _opMap.keySet()) {
//            sb.append(_opMap.get(key).toString());
//            sb.append("\n");
//        }
//        return sb.toString();
//    }
//
//
//    /**
//     * Running the MCMC chain from scratch.
//     */
//    public void firstRun() {
//        int total = (int) (Utils._CHAIN_LEN / Utils._SAMPLE_FREQUENCY);
//        int burnin = (int) (Utils._BURNIN_LEN / Utils._SAMPLE_FREQUENCY);
//        int swap = (int) (Utils._SAMPLE_FREQUENCY / Utils.SWAP_FREQUENCY);
//        if (Utils._PRE_BURN_IN) {
//            Utils._NET_MAX_RETI = 0;
//        }
//        for(int i = 1 - Utils._PRE_BURN_IN_ITER; i <= total; i++) {
//            if (i == 1) {
//                Utils._PRE_BURN_IN = false;
//                Utils._NET_MAX_RETI = MaxReticulation;
//                MC3 main = getMain();
//                if(main == null) {
//                    throw new RuntimeException("No main MC3 chain found!!!");
//                }
//                main.setPreBurnInParams();
//                System.out.println("----------------------- Logger: -----------------------");
//                System.out.println("Iteration;    Posterior;  ESS;    Likelihood;    Prior;  ESS;    #Reticulation");
//            }
//            if(i > burnin) {
//                _sampling = true;
//            }
//            // swap state
//            for(int k = 0; k < swap; k++) {
//                runMC3s(i, k==0);
//                if(_mc3s.size() > 1) doSwap();
//            }
//        }
//        outputLastState(Utils._OUT_DIRECTORY);
//        summarize();
//    }
//
//    /**
//     * Resume a chain from where it was stopped.
//     */
//    public void multipleRun() {
//        int total = (int) (Utils._CHAIN_LEN / Utils._SAMPLE_FREQUENCY);
//        int swap = (int) (Utils._SAMPLE_FREQUENCY / Utils.SWAP_FREQUENCY);
//        if (Utils._PreRun) {
//            Utils._PRE_BURN_IN = false;
//            Utils._PRE_BURN_IN_ITER = 0;
//            _sampling = true;
//        }
//        for(int i = 1; i <= total; i++) {
//            if (i == 1) {
//                Utils._PRE_BURN_IN = false;
//                Utils._NET_MAX_RETI = MaxReticulation;
//                MC3 main = getMain();
//                if(main == null) {
//                    throw new RuntimeException("No main MC3 chain found!!!");
//                }
//                main.setPreBurnInParams();
//                System.out.println("----------------------- Logger: -----------------------");
//                System.out.println("Iteration;    Posterior;  ESS;    Likelihood;    Prior;  ESS;    #Reticulation");
//            }
//
//            // swap state
//            for(int k = 0; k < swap; k++) {
//                runMC3s(i, k==0);
//                if(_mc3s.size() > 1) doSwap();
//            }
//        }
//        summarize();
//        outputLastState(Utils._OUT_DIRECTORY);
//    }
//
//    public String getLastNetwork() {
//        for(MC3 mc3 : _mc3s) {
//            if(mc3.getMain()) return mc3.getCurrentNetwork();
//        }
//        return null;
//    }
//
//    private MC3 getMain() {
//        for(MC3 mc3 : _mc3s) {
//            if(mc3.getMain()) return mc3;
//        }
//        return null;
//    }
//
//    public void outputLastState(String outputFile) {
//        for (int index = 0; index < _mc3s.size(); index ++) {
//            MC3 mc3 = _mc3s.get(index);
//            State ss = mc3.getState();
//            String out0 = outputFile;
//            if (_mc3s.size() > 1)
//                out0 += String.valueOf(index) + "_";
//            try {
//                BufferedWriter outAll = new BufferedWriter(new FileWriter(out0 + "lastState.txt"));
//                BufferedWriter outCheck = new BufferedWriter(new FileWriter(out0+"check.txt"));
//
//                // ------------- output gene trees ----------------
//                List<UltrametricTree> gt = ss.getGeneTrees();
//                for(int i =0; i<gt.size();i++){
//                    outAll.write(gt.get(i).toString()+"\n");
//                    //TODO compute nodeheight when read
//                }
//                outAll.flush();
//                for(int i =0; i< gt.size();i++){
//                    outAll.write(gt.get(i).get_mutationRate()+"\n");
//
//                }
//                outAll.flush();
//                // ------------- output population size -----------
//                UltrametricNetwork ultraNet = ss.getUltrametricNetworkObject();
//                BniNetwork<NetNodeInfo> speciesNet = (BniNetwork<NetNodeInfo>)ultraNet.getNetwork();
//                Networks.autoLabelNodes(speciesNet);
//
//                // ------------- output species network with root population size -----------
//                outAll.write("[" +speciesNet.getRoot().getRootPopSize() + "]"
//                        + speciesNet.toString()+ "\n");
//                //------------- output speices network NetNodeInfo -----------
//                for (NetNode<NetNodeInfo> node : speciesNet.dfs()) {
//                    StringBuilder s = new StringBuilder(node.getName());
//                    NetNodeInfo tempInfo = node.getData();
//                    s.append("\t").append(tempInfo.getHeight()).append("\t").append(tempInfo.getIndex()).append("\t").append(tempInfo.getPrevHeight());
//                    if (Utils._ESTIMATE_POP_SIZE && !Utils._CONST_POP_SIZE){
//                        for (NetNode<NetNodeInfo> par : node.getParents()) {
//                            double ps = node.getParentSupport(par);
//                            if(Double.isNaN(ps) || ps <= 0) {
//                                System.err.println("Wrong pop size: " + node.getName() + "<-" + par.getName()+":"+ps);
//                            }
//                            s.append("\t" + "[").append(par.getName()).append(",").append(ps).append("]");
//                        }
//
//                    }
//                    s.append("\n");
//                    outAll.write(s.toString());
//                }
//
//                //------------- output popSizeSample, pop size mean ---------------
//                if(Utils._ESTIMATE_POP_SIZE) {
//                    outAll.write(ss.getPopSize()+"\n");
//                }
//
//
//                outAll.write(Utils._SEED + "\n");
//                outAll.flush();
//                outAll.close();
//                outCheck.write(mc3.getLikelihood() + "\n");
//                outCheck.write(mc3.getPosterior() + "\n");
//                outCheck.write(mc3.getPrior() + "\n");
//                outCheck.write(mc3.getESS() + "\n");
//                outCheck.flush();
//                outCheck.close();
//
//
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }
//    }
//
//}
