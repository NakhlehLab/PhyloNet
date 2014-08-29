package edu.rice.cs.bioinfo.programs.phylonet.algos.integration;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.FelsensteinAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

/**
 * Created by yunyu on 6/3/14.
 */
public class FullLikelihoodFromSequence {
    List<Tree> _gts;
    double _minHeight;
    double _maxHeight;
    double _theta;
    SubstitutionModel _model;
    //Network _network;
    //String[] _networkTaxa;
    String[] _gtTaxa;
    List<Tuple<char[],Integer>> _sequences;
    Long _seed = null;
    int _counter = 0;
    int _currentID = 0;
    int _numThreads = 1;
    boolean _printDetails = false;

    int _numSamples;
    int _numSamplesPerBin;
    int _numBins;
    int _dim;

    short[][] _samples;
    //double[][] _samples;
    float[][][] _gtLikelihoods;

    public void setRandomSeed(Long seed){
        _seed = seed;
    }

    public void setPrintDetails(boolean ifPrint){
        _printDetails = ifPrint;
    }

    public void setParallel(int numThreads){
        _numThreads = numThreads;
    }

    public FullLikelihoodFromSequence(List<Tree> gts, String[] gtTaxa, List<Tuple<char[], Integer>> sequences, SubstitutionModel model, double theta, double minHeight, double maxHeight){
        //_network = network;
        _gts = gts;
        //_species2alleles = species2alleles;
        _minHeight = minHeight;
        _maxHeight = maxHeight;
        _model = model;
        _theta = theta;
        _gtTaxa = gtTaxa;
        _sequences = sequences;
    }



    private boolean setTreeBranchLength(Tree tr, short[] nodeHeight, double scalor){
        int index = 0;
        Map<TNode, Double> node2height = new HashMap<TNode, Double>();
        for (TNode node : tr.postTraverse()) {
            double height = 0;
            if (!node.isLeaf()) {
                height = nodeHeight[index++]/1000.0;
            }
            for (TNode child : node.getChildren()) {
                double bl = height - node2height.get(child);
                if(bl < 0){
                    return false;
                }
                child.setParentDistance(bl * scalor);
            }
            node2height.put(node, height);
        }
        return true;
    }




    private float[] computeGTLikelihood(Tree gt){
        float[] results = new float[_sequences.size()];
        int index = 0;
        for(Tuple<char[], Integer> seq: _sequences) {
            Map<String, Character> sequenceMap = new HashMap<String, Character>();
            for(int i=0; i<_gtTaxa.length; i++){
                sequenceMap.put(_gtTaxa[i], seq.Item1[i]);
            }

            OneNucleotideObservation converter = new OneNucleotideObservation(sequenceMap);

            FelsensteinAlgorithm fcalc = new FelsensteinAlgorithm(gt,_model);
            results[index++] = (float)fcalc.getProbability(converter);
            //Felsenstein fcalc = new Felsenstein(_model);
            //results[index++] = (float)fcalc.getLikelihoodtree(gt, new ObservationMap(sequenceMap));
        }
        return results;
    }

    private void clearParallelIndex(){
        _currentID = 0;
    }

    public void computeGTLikelihoodWithIntegral(int numSamplesPerBin, int numBins){
        _numBins = numBins;
        _numSamplesPerBin = numSamplesPerBin;
        _dim = _gts.get(0).getNodeCount() - _gts.get(0).getLeafCount();
        _numSamples = (int)Math.pow(_numBins, _dim) * numSamplesPerBin;
        _samples = new short[_numSamples][_dim];
        createSamples();
        _gtLikelihoods = new float[_gts.size()][_numSamples][_sequences.size()];


        clearParallelIndex();
        Thread[] myThreads = new Thread[_numThreads];
        for(int i=0; i<_numThreads; i++){
            myThreads[i] = new GTLikelihoodThread();
            myThreads[i].start();
        }

        for(int i=0; i<_numThreads; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }
    }


    public double computeNetworkLikelihoodWithIntegral(Network network, Map<String, List<String>> species2alleles) {
        clearParallelIndex();
        Thread[] myThreads = new Thread[_numThreads];
        String networkExp = network.toString();
        double[][] networkLikelihood = new double[_gts.size()][_numSamples];
        for (int i = 0; i < _numThreads; i++) {
            myThreads[i] = new NetworkLikelihoodThread(networkExp, species2alleles, networkLikelihood);
            myThreads[i].start();
        }

        for (int i = 0; i < _numThreads; i++) {
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {
            }
        }

        double finalResult = 0;
        for(int i=0; i<_sequences.size(); i++){
            double siteProb = 0;
            for(int j=0; j<_gts.size(); j++){
                double[] integralValues = new double[_numSamples];
                for(int k=0; k<_numSamples; k++){
                    integralValues[k] = _gtLikelihoods[j][k][i] * networkLikelihood[j][k];
                    if(_printDetails && _gtLikelihoods[j][k][i]!=0){
                        Tree gt = _gts.get(j);
                        setTreeBranchLength(gt, _samples[k], 1);
                        //System.out.println(gt + " : " + _gtLikelihoods[j][k][i] + " * " + networkLikelihood[j][k]);
                    }

                }
                siteProb += MultivariateMonteCarloIntegral.integrate(_dim, _minHeight, _maxHeight, integralValues);
            }
            //System.out.println(Math.log(siteProb) + " * " + _sequences.get(i).Item2);
            //System.out.println(siteProb);
            //finalResult += Math.log(siteProb) * _sequences.get(i).Item2;
            finalResult += siteProb * _sequences.get(i).Item2;
            System.out.println(siteProb);
        }

        //System.out.println(finalResult);
        return finalResult;
    }




    private void createSamples(){
        MultivariateMonteCarloIntegral integration = new MultivariateMonteCarloIntegral(_numSamplesPerBin, _numBins);
        integration.setRandomSeed(_seed);
        int j = 0;
        for (double[] sample : integration.getSamplePoints(_dim, _minHeight, _maxHeight)) {
            for(int k=0; k<sample.length; k++) {
                _samples[j][k] = (short)(sample[k]*1000);
            }
            j++;
        }
    }

    private class GTLikelihoodThread extends Thread{

        public void run() {
            int myID;
            while((myID = getNextID(_gts.size()))!=-1){
                float[][] localLikelihoods = new float[_numSamples][_sequences.size()];
                Tree gt = _gts.get(myID);
                //GTBranchLengthsValidator validator = new GTBranchLengthsValidator(gt, _dim);
                int j = 0;

                for (short[] sample : _samples) {
                    boolean valid = setTreeBranchLength(gt, sample, _theta);
                    //_samples[myID][j] = gtCopy;

                    if (!valid) {
                        Arrays.fill(localLikelihoods[j], 0);
                        //Arrays.fill(_gtLikelihoods[myID][j], 0);
                    } else {
                        localLikelihoods[j] = computeGTLikelihood(gt);
                        //_gtLikelihoods[myID][j] = computeGTLikelihood(gt);

                        //System.out.println(gt);
                    }
                    //System.out.println(myID + ": " + j);
                    j++;
                }
                _gtLikelihoods[myID] = localLikelihoods;
            }
        }

        private synchronized int getNextID(int total){

            if(_currentID >= total){
                return -1;
            }
            else{
                _currentID++;
                return _currentID-1;
            }
        }
    }



    private class NetworkLikelihoodThread extends Thread{
        Network _network;
        Map<String, List<String>> _species2alleles;
        double[][] _results;

        public NetworkLikelihoodThread(String network, Map<String, List<String>> species2alleles, double[][] results){
            _network = Networks.readNetwork(network);
            _species2alleles = species2alleles;
            _results = results;
        }


        public void run() {
            int myID;
            while((myID = getNextID(_gts.size()))!=-1) {
                List<Tree> gtList = new ArrayList<Tree>();
                Tree gt = _gts.get(myID);
                gtList.add(gt);
                GeneTreeWithBranchLengthProbabilityYF gtp = new GeneTreeWithBranchLengthProbabilityYF(_network, gtList, _species2alleles);

                int index = 0;
                for (short[] sample : _samples) {
                    if(_gtLikelihoods[myID][index][0] != 0) {
                        setTreeBranchLength(gt, sample, 1);
                        double[] probs = new double[1];
                        gtp.calculateGTDistribution(probs);
                        _results[myID][index] = probs[0];
                    }
                    index++;
                }
            }
        }



        private synchronized int getNextID(int total){
            if(_currentID >= total){
                return -1;
            }
            else{
                _currentID++;
                return _currentID-1;
            }
        }
    }


    public static void main(String[] args){
        double sum =0;
        char[] nucleotides = {'A','C','T','G'};

        List<Tuple<char[],Integer>> sequences = new ArrayList<Tuple<char[], Integer>>();
        for (char a : nucleotides)
            for (char b: nucleotides)
                for (char c: nucleotides)
                    for (char d: nucleotides)
                {
                    //Network speciesNetworkTopology = Networks.string2network("(A:.5,(B:.2,C:.2)n1:.3)n0;");
                    char[] seq = new char[4];
                    seq[0] = a;
                    seq[1] = b;
                    seq[2] = c;
                    seq[3] = d;
                    sequences.add(new Tuple<char[], Integer>(seq,1));
                }


        String[] taxa = {"A","B","C","D"};
        GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});
        //Network speciesNetworkTopology = Networks.string2network("((A:.2,X#H1:0::0.2)n2:0.3,((B:.1)X#H1:0::0.8,C:.2)n1:.3)root;");
        Network speciesNetworkTopology = Networks.readNetwork("((A:1,X#H1:1::0.3)n2:1,(((B:1,C:1)n3:1)X#H1:1::0.7,D:1)n1:1)root;");
        //Network speciesNetworkTopology = Networks.string2network("((A:2,X#H1:1::0.7)n3:1,((B:1)X#H1:1::0.3,C:2)n1:1)root;");
        FullLikelihoodFromSequence likelihood = new FullLikelihoodFromSequence(Trees.generateAllBinaryTrees(taxa),taxa,sequences,gtrModel,0.02,0,10);
        likelihood.computeGTLikelihoodWithIntegral(1000,10);
        //likelihood.computeNetworkLikelihoodWithIntegral(speciesNetworkTopology,null);
        System.out.println(likelihood.computeNetworkLikelihoodWithIntegral(speciesNetworkTopology,null));
        //System.out.println("Sum: " + sum);
    }

}
