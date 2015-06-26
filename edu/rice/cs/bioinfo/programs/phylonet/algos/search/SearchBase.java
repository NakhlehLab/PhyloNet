package edu.rice.cs.bioinfo.programs.phylonet.algos.search;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by yunyu on 10/27/14.
 */
public abstract class SearchBase {
    protected LinkedList<Tuple<Network,Double>> _optimalNetworks;
    protected int _numOptimums;
    private File _logFile = null;
    private File _resultFile = null;
    private List<Tuple<Network,Double>> _networksTried = new ArrayList<>();
    protected boolean _scoreEachTopologyOnce;
    private boolean _printDetails = false;
    private Comparator<Double> _scoreComparator;

    public SearchBase(Comparator<Double> scoreComparator){
        _scoreComparator = scoreComparator;
    }

    protected int compareTwoScores(double score1, double score2){
        return _scoreComparator.compare(score1, score2);
    }

    protected boolean printDetails(){
        return _printDetails;
    }

    protected boolean scoreEachTopologyOnce(){
        return _scoreEachTopologyOnce;
    }

    public void setLogFile(File file) {
        _logFile = file;
        if (_logFile != null)
            readPreviousTriedNetwork();
    }


    protected int getCashedNetworkListSize(){
        return _networksTried.size();
    }

    protected List<Tuple<Network,Double>> getAllCashedNetworks(){
        return _networksTried;
    }

    private void readPreviousTriedNetwork(){
        try{
            if(_logFile.exists()) {
                BufferedReader br = new BufferedReader(new FileReader(_logFile));
                String line;
                while ((line=br.readLine())!=null){
                    String[] fields = line.split(" ");
                    _networksTried.add(new Tuple<Network, Double>(Networks.readNetwork(fields[1].trim()),Double.parseDouble(fields[0])));
                }
                br.close();
            }
        }catch(Exception e) {
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
    }

    protected void writeLogFile(Network network, double score){
        if(_logFile!=null) {
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(_logFile, true));
                bw.append(score + " " + network.toString() + "\n");
                bw.close();
            } catch (Exception e) {
                System.err.println(e.getMessage());
                e.getStackTrace();
            }
        }
    }

    protected Double getPreviousResult(Network network){
        for(Tuple<Network,Double> networkTried: _networksTried){
            if(Networks.hasTheSameTopology(networkTried.Item1, network)){
                return networkTried.Item2;
            }
        }
        return null;
    }

    public void setResultFile(File file){
        _resultFile = file;
    }


    protected int readResultsFromPreviousRuns(){
        int numRuns = 0;
        if(_resultFile!=null){
            try {
                BufferedReader br = new BufferedReader(new FileReader(_resultFile));
                String line;
                List<String> results = new ArrayList<>();
                while((line=br.readLine()) != null){
                    if(line.startsWith("Results after run")){
                        results.clear();
                        while(!(line=br.readLine()).startsWith("===============")){
                            results.add(line);
                        }
                        numRuns++;
                    }
                }
                for(String result: results) {
                    int index = result.indexOf(":");
                    updateOptimalNetworks(Networks.readNetwork(result.substring(index + 1).trim()),Double.parseDouble(result.substring(0, index).trim()));
                }
                br.close();
            }catch (Exception e) {
                System.err.println(e.getMessage());
                e.getStackTrace();
            }
        }
        return numRuns;
    }


    protected void writeResultFile(int numRun){
        try{
            if(_resultFile != null) {
                BufferedWriter bw = new BufferedWriter(new FileWriter(_resultFile, true));
                bw.append("Results after run #" + numRun + "\n");
                for (Tuple<Network, Double> result: _optimalNetworks) {
                    bw.append(result.Item2 + ": " + result.Item1.toString() + "\n");
                }
                //System.out.println("Running Time: " + (System.currentTimeMillis()-startingTime)/60000.0);
                bw.append("===============================\n\n");
                bw.close();
            }
        }catch(Exception e) {
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
    }




    protected boolean updateOptimalNetworks(Network network, double score){
        if(_optimalNetworks.isEmpty()){
            _optimalNetworks.add(new Tuple<>(network.clone(), score));
            return true;
        }

        /*
        Network trueNetwork = Networks.readNetwork("((((Zrou:1.0,(Kpol:1.0,(Cgla:1.0,(Scas:1.0,(Sbay:1.0,(Skud:1.0,(Smik:1.0,(Scer:1.0,Spar:1.0)I21:1.0)I13:1.0)I23:1.0)I11:1.0)I8:1.0)I15:1.0)I24:1.0)I4:1.0)I2#H1:1.0::0.74,((Agos:1.0,Klac:1.0)I17:1.0,((Kthe:1.0,Kwal:1.0)I25:1.0,(Sklu:1.0,I2#H1:1.0::0.26)I10:1.0)I7:1.0)I18:1.0)I0:1.0,((((Ctro:1.0,(Calb:1.0,Cdub:1.0)I22:1.0)I19:1.0,(Lelo:1.0,Cpar:1.0)I12:1.0)I16:1.0)I5#H2:1.0::0.44,(Clus:1.0,(Cgui:1.0,(Dhan:1.0,(Psti:1.0,I5#H2:1.0::0.56)I6:1.0)I14:1.0)I9:1.0)I20:1.0)I3:1.0)I1;");
        boolean discovered = false;
        if(Networks.hasTheSameTopology(trueNetwork, network)){
            System.out.println("===================================");
            System.out.println("True network discovered: " + score);
            discovered = true;
        }
        */

        if(compareTwoScores(score, _optimalNetworks.getLast().Item2)>0) {
            for (int i = 0; i < _optimalNetworks.size(); i++) {
                Tuple<Network, Double> entry = _optimalNetworks.get(i);
                if (Networks.hasTheSameTopology(network, entry.Item1)) {
                    if (compareTwoScores(score, entry.Item2)>0){
                        if ((i == 0 || (i != 0 && compareTwoScores(_optimalNetworks.get(i-1).Item2, score)>0)) && (i == _optimalNetworks.size() - 1 || (i != _optimalNetworks.size() - 1 && compareTwoScores(_optimalNetworks.get(i + 1).Item2, score)<0))) {
                            _optimalNetworks.set(i, new Tuple<>(network.clone(), score));
                            return true;
                        }
                        _optimalNetworks.remove(i);
                        break;
                    } else {
                        return false;
                    }
                }
            }
            for (int i = 0; i < _optimalNetworks.size(); i++) {
                if(compareTwoScores(score, _optimalNetworks.get(i).Item2)>0){
                    _optimalNetworks.add(i, new Tuple<>(network.clone(), score));
                    break;
                }

            }

            if(_optimalNetworks.size()>_numOptimums){
                _optimalNetworks.removeLast();
            }

            /*
            if(discovered){
                System.out.println("Added to optimum ");
                System.out.println("===================================");
            }
            */

            return true;
        }
        else{
            if(_optimalNetworks.size()<_numOptimums){
                for (int i = 0; i < _optimalNetworks.size(); i++) {
                    Tuple<Network, Double> entry = _optimalNetworks.get(i);
                    if (Networks.hasTheSameTopology(network, entry.Item1)) {
                        return false;
                    }
                }
                _optimalNetworks.add(new Tuple<>(network.clone(), score));

                /*
                if(discovered){
                    System.out.println("Added to optimum ");
                    System.out.println("===================================");
                }
                */
                return true;
            }
            else{
                /*
                if(discovered){
                    System.out.println("Worse than tops ");
                    System.out.println("===================================");
                }
                */
                return false;
            }
        }


    }

    protected void updateCashedResults(Network network, double score) {
        if (_scoreEachTopologyOnce) {
            _networksTried.add(new Tuple<Network, Double>(network.clone(), score));
        }
    }


    protected abstract boolean concludeSearch();
}
