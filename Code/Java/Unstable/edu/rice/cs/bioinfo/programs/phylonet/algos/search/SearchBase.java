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
 * Created by yy9
 * Date: 5/10/12
 * Time: 2:31 PM
 *
 * This class is the base of search
 */
public abstract class SearchBase {
    protected LinkedList<Tuple<Network,Double>> _optimalNetworks;
    protected int _numOptimums;
    private File _logFile = null;
    private File _intermediateResultFile = null;
    private List<Tuple<Network,Double>> _networksTried = new ArrayList<>();    //only used when _scoreEachTopologyOnce is true
    protected boolean _scoreEachTopologyOnce;
    private boolean _printDetails = false;
    private Comparator<Double> _scoreComparator;


    /**
     * Constructor / Initializer
     */
    public SearchBase(Comparator<Double> scoreComparator){
        _scoreComparator = scoreComparator;
    }


    /**
     * This function is to compare two scores of two candidate species networks
     *
     * @return 1 if score1 is better; 0 if they are equal; -1 if score2 is better
     */
    protected int compareTwoScores(double score1, double score2){
        return _scoreComparator.compare(score1, score2);
    }


    /**
     * This function is to set the option of printing details
     */
    protected boolean printDetails(){
        return _printDetails;
    }


    /**
     * This function is to set the option of scoring each topology once.
     * It can be used in MP, or ML when branch lengths and inheritance probabilities are optimized for each topology
     */
    protected boolean scoreEachTopologyOnce(){
        return _scoreEachTopologyOnce;
    }


    /**
     * This function is to set the log file which saves all the networks that have been tried during the search
     */
    public void setLogFile(File file) {
        _logFile = file;
        if (_logFile != null)
            readPreviousTriedNetwork();
    }


    /**
     * This function is to return the number of networks that have been tried
     */
    protected int getCashedNetworkListSize(){
        return _networksTried.size();
    }


    /**
     * This function is to return all the networks that have been tried
     */
    protected List<Tuple<Network,Double>> getAllCashedNetworks(){
        return _networksTried;
    }



    /**
     * This function is to read all the networks that have been tried from the log file
     */
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



    /**
     * This function is to write the network that has been tried to the log file
     *
     * @param network   the species network that has been tried
     * @param score     the score of the network
     */
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



    /**
     * This function is to get the score of the species network if this network has already been visited before
     */
    protected Double getPreviousResult(Network network){
        for(Tuple<Network,Double> networkTried: _networksTried){
            if(Networks.hasTheSameTopology(networkTried.Item1, network)){
                return networkTried.Item2;
            }
        }
        return null;
    }



    /**
     * This function is to set the intermediate result file which saves the results from every run
     */
    public void setIntermediateResultFile(File file){
        _intermediateResultFile = file;
    }


    /**
     * This function is to read the results of previous runs from intermediate result file
     */
    protected int readResultsFromPreviousRuns(){
        int numRuns = 0;
        if(_intermediateResultFile!=null){
            try {
                BufferedReader br = new BufferedReader(new FileReader(_intermediateResultFile));
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


    /**
     * This function is to write the results of the current run to the intermediate result file
     */
    protected void writeIntermediateResultFile(int numRun){
        try{
            if(_intermediateResultFile != null) {
                BufferedWriter bw = new BufferedWriter(new FileWriter(_intermediateResultFile, true));
                bw.append("Results after run #" + numRun + "\n");
                for (Tuple<Network, Double> result: _optimalNetworks) {
                    bw.append(result.Item2 + ": " + result.Item1.toString() + "\n");
                }
                bw.append("===============================\n\n");
                bw.close();
            }
        }catch(Exception e) {
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
    }



    /**
     * This function is to update the list of optimal networks if needed
     *
     * @param network   the current network that is being visited
     * @param score     the score of the current network
     */
    protected boolean updateOptimalNetworks(Network network, double score){
        if(_optimalNetworks.isEmpty()){
            _optimalNetworks.add(new Tuple<>(network.clone(), score));
            return true;
        }

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
                return true;
            }
            else{
                return false;
            }
        }
    }


    /**
     * This function is to update the list of networks that have been tried along with the likelihood
     *
     * @param network   the current network that is being visited
     * @param score     the score of the current network
     */
    protected void updateCashedResults(Network network, double score) {
        if (_scoreEachTopologyOnce) {
            _networksTried.add(new Tuple<Network, Double>(network.clone(), score));
        }
    }


    /**
     * This function is to check if the search should be terminated
     */
    protected abstract boolean concludeSearch();
}
