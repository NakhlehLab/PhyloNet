package edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.util.Comparator;
import java.util.LinkedList;

/**
 * Created by yunyu on 11/4/14.
 * It implements the simulated annealing in Salter and Pearl 2001, which is the same as the one in STEM
 */

public class SimulatedAnnealingSalterPearL extends SimulatedAnnealingBase{
    private double _U;
    private double _beta;
    private double _burnin = 1000;
    private double _preScore = -1;
    private boolean _doneBurnin = false;
    private LinkedList<Tuple<Network,MutableTuple<Double,Integer>>> _savedNetworks;
    private int _numSavedNetworks = 10;
    private int _maxUninvestigatedProposals;
    private int _consecutiveUninvestigatedProposals=0;
    private int _investigatedProposalThreshold;
    private boolean _initialized = false;

    //Temp: for intermediateFile only
    private long _previousTime = -1;
    private File _intermediateFile;
    private String _currentNetwork;



    public SimulatedAnnealingSalterPearL(Comparator<Double> scoreComparator, NetworkNeighbourhoodGenerator generator, Long seed) {
        super(scoreComparator, generator, seed);

    }

    public void setIntermediateFile(String fileName){
        _intermediateFile = new File(fileName);
        readIntermediateFile();
    }

    public void writeIntermediateFile(){
        if(_intermediateFile == null || !_intermediateFile.exists()) return;
        long currentTime = System.currentTimeMillis();
        if(_previousTime == -1){
            _previousTime = currentTime;
        }
        else if((currentTime - _previousTime) / 60000.0 > 10 && _doneBurnin){
            _previousTime = currentTime;
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(_intermediateFile));
                bw.append("_currentNetwork:" + _currentNetwork+"\n");
                bw.append("_U:" + _U + "\n");
                bw.append("_beta:" + _beta + "\n");
                //bw.append("_preScore:" + _preScore + "\n");
                bw.append("_savedNetworks:\n");
                for(Tuple<Network,MutableTuple<Double,Integer>> tuple: _savedNetworks){
                    bw.append(tuple.Item1.toString() + " " + tuple.Item2.Item1 + " " + tuple.Item2.Item2 + "\n");
                }
                bw.append("_maxUninvestigatedProposals:" + _maxUninvestigatedProposals + "\n");
                bw.append("_consecutiveUninvestigatedProposals:" + _consecutiveUninvestigatedProposals + "\n");
                bw.append("_investigatedProposalThreshold:" + _investigatedProposalThreshold + "\n");
                bw.append("_temperature:" + _temperature + "\n");
                bw.append("_examinationsCount:" + getExaminations() + "\n");
                bw.append("_optimalNetworks:\n");
                for(Tuple<Network,Double> tuple: _optimalNetworks){
                    bw.append(tuple.Item1.toString() + " " + tuple.Item2+ "\n");
                }
                bw.close();
            }catch (Exception e){
                System.err.println(e.getMessage());
                e.getStackTrace();
            }
            //System.out.println("\n\n\nWriting files ........................................................ \n\n");
        }
        else{
            //System.out.println((currentTime - _previousTime) / 60000.0);
        }
    }

    public void readIntermediateFile() {
        try {
            if(!_intermediateFile.exists()) return;
            String line;
            BufferedReader br = new BufferedReader(new FileReader(_intermediateFile));
            if(br.readLine()==null) return;
            _initialized = true;
            _doneBurnin = true;
            line = br.readLine();

            _U = Double.parseDouble(line.substring(line.indexOf(":") + 1));
            line = br.readLine();
            _beta = Double.parseDouble(line.substring(line.indexOf(":") + 1));
            //line = br.readLine();
            //_preScore = Double.parseDouble(line.substring(line.indexOf(":") + 1));
            br.readLine();
            _savedNetworks = new LinkedList<>();
            for(int i=0; i<_numSavedNetworks; i++){
                String[] subs = br.readLine().split(" ");
                _savedNetworks.add(new Tuple<Network, MutableTuple<Double, Integer>>(Networks.readNetwork(subs[0]), new MutableTuple<Double, Integer>(Double.parseDouble(subs[1]), Integer.parseInt(subs[2]))));
            }
            line = br.readLine();
            _maxUninvestigatedProposals = Integer.parseInt(line.substring(line.indexOf(":") + 1));
            line = br.readLine();
            _consecutiveUninvestigatedProposals = Integer.parseInt(line.substring(line.indexOf(":") + 1));
            line = br.readLine();
            _investigatedProposalThreshold = Integer.parseInt(line.substring(line.indexOf(":") + 1));
            line = br.readLine();
            _temperature = Double.parseDouble(line.substring(line.indexOf(":") + 1));
            line = br.readLine();
            setExaminations(Integer.parseInt(line.substring(line.indexOf(":") + 1)));
            br.readLine();
            _optimalNetworks = new LinkedList<>();
            while((line = br.readLine())!=null){
                String[] subs = line.split(" ");
                _optimalNetworks.add(new Tuple<Network, Double>(Networks.readNetwork(subs[0]), Double.parseDouble(subs[1])));
            }
            br.close();
        } catch (Exception e) {
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
    }


    private void clearConsecutiveUninvestigatedProposals(){
        _consecutiveUninvestigatedProposals = 1;
    }

    private void increaseConsecutiveUninvestigatedProposals(){
        _consecutiveUninvestigatedProposals++;
    }

    private boolean reachMaxUninvestigatedProposals(){
        return _consecutiveUninvestigatedProposals>=_maxUninvestigatedProposals;
    }

    protected double computeRandomNeighborScore(Network currentNetwork, Func1<Network, Double> getScore) {
        //System.out.println("Start computing");
        double newScore = super.computeRandomNeighborScore(currentNetwork, getScore);
        if(!_doneBurnin){
            if(!Double.isNaN(newScore)){
                _U = Math.max(_U, Math.abs(newScore-_preScore));
                _preScore = newScore;
            }
            if(getExaminations()==_burnin){
                _doneBurnin = true;
                clearExaminations();
                if(printDetails()){
                    System.out.println("Done Burn-in: " + _U);
                    System.out.println();
                }
            }
        }
        updateTemperature();

        writeIntermediateFile();
        return newScore;
    }

    protected void updateTemperature(){
        _temperature = _U/(1+getExaminations()*_beta);
    }


    protected void search(Network currentNetwork, Func1<Network,Double> getScore, double initialScore) {
        if (!_initialized) {
            initializeParameters(currentNetwork, initialScore);
            _initialized = true;
        }
        super.search(currentNetwork, getScore, initialScore);
    }


    protected void initializeTemperature(){}

    protected void initializeParameters(Network currentNetwork, double initialScore){
        int ntaxa = currentNetwork.getLeafCount();
        _maxUninvestigatedProposals = 100 * ntaxa;
        _consecutiveUninvestigatedProposals=0;
        _investigatedProposalThreshold = (int) (Math.log(0.05) / (ntaxa * Math.log(1.0 - (1.0 / ntaxa)))) + 1;
        _investigatedProposalThreshold *= ntaxa;
        if (scoreEachTopologyOnce()) {
            _investigatedProposalThreshold *= 4;
        } else {
            _investigatedProposalThreshold *= 12;
        }
        _U = Math.abs(initialScore) / 4;
        setBeta(initialScore, ntaxa, 2000);
        _preScore = initialScore;
        if(_numSavedNetworks!=0){
            _savedNetworks = new LinkedList<>();
        }
        _doneBurnin = false;
        if(printDetails()){
            System.out.println("Initial setting:");
            System.out.println("MaxUninvestigatedProposals: " + _maxUninvestigatedProposals);
            System.out.println("InvestigatedProposalThreshold: " + _investigatedProposalThreshold);
            System.out.println("U: " + _U);
            System.out.println("Beta: " + _beta);
            System.out.println();
        }
    }

    protected boolean concludeSearch(){
        if(printDetails()){
            System.out.println("Consecutive unvestigated network: " + _consecutiveUninvestigatedProposals);
        }

        boolean conclude = reachMaxExaminations() || reachMaxUninvestigatedProposals();

        if(conclude){
            /*
            System.out.println("Consecutive unvestigated network: " + _consecutiveUninvestigatedProposals);
            System.out.println("Iterations: " + getExaminations());
            for(Tuple<Network,MutableTuple<Double,Integer>> network: _savedNetworks){
                System.out.println(network.toString());
            }
            */
            _initialized = false;
        }

        return conclude;
    }


    protected void setBeta(double lnl, int n, int m){
        /*
        double c = 0.5;
        double alpha = 0.5;
        _beta = c/((1-alpha)*n + alpha*(-lnl)/m);
        */
        _beta = 0.01;
    }

    private void updateSavedNetworks(Network newProposal, double score){
        double worstScore = 0;
        int worstIndex = -1;
        int index = 0;
        for(Tuple<Network,MutableTuple<Double,Integer>> network: _savedNetworks){
            if(Networks.hasTheSameTopology(newProposal,network.Item1)){
                if(compareTwoScores(score,network.Item2.Item1)>0){
                    network.Item2.Item1 = score;
                    network.Item2.Item2 = 1;
                    if(printDetails()){
                        System.out.println("Exist in saved networks but higher prob");
                    }
                    clearConsecutiveUninvestigatedProposals();
                    return;
                }
                else{
                    network.Item2.Item2++;
                    if(network.Item2.Item2>_investigatedProposalThreshold){
                        if(printDetails()){
                            System.out.println("Exist in saved networks & exceed the threshold ");
                        }
                        increaseConsecutiveUninvestigatedProposals();
                        return;
                    }
                    else{
                        clearConsecutiveUninvestigatedProposals();
                        if(printDetails()){
                            System.out.println("Exist in saved networks: " + network.Item2.Item1);
                        }
                        return;
                    }
                }
            }
            if(index==0 || compareTwoScores(worstScore, network.Item2.Item1)>0){
                worstIndex = index;
                worstScore = network.Item2.Item1;
            }
            index++;
        }

        if(_savedNetworks.size()==_numSavedNetworks){
            _savedNetworks.remove(worstIndex);
        }
        _savedNetworks.add(new Tuple<Network, MutableTuple<Double, Integer>>(newProposal.clone(), new MutableTuple<Double, Integer>(score,1)));
        if(printDetails()){
            System.out.println("Add to saved networks");
        }
        clearConsecutiveUninvestigatedProposals();
    }


    protected void handleAcceptCase(Network currentNetwork, Ref<Double> currentScore, double newScore){
        _currentNetwork = currentNetwork.toString();
        currentScore.set(newScore);
        if(_numSavedNetworks!=0) {
            updateSavedNetworks(currentNetwork, newScore);
        }
    }

    protected Network handleRejectCase(Network currentNetwork, Ref<Double> currentScore, double newScore){
        _networkGenerator.undo();
        increaseConsecutiveUninvestigatedProposals();
        return currentNetwork;
    }



}
