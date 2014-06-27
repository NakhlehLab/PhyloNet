package edu.rice.cs.bioinfo.programs.phylonet.algos.integration;

import edu.rice.cs.bioinfo.library.programming.Func1;
  import edu.rice.cs.bioinfo.library.programming.MutableTuple;
  import edu.rice.cs.bioinfo.library.programming.Tuple;
  import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
  import edu.rice.cs.bioinfo.programs.phmm.src.phylogeny.Felsenstein;
  import edu.rice.cs.bioinfo.programs.phmm.src.substitutionModel.GTRSubstitutionModel;
  import edu.rice.cs.bioinfo.programs.phmm.src.substitutionModel.SubstitutionModel;
  import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbability;
  import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
  import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
  import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
  import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
  import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
  import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
  import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
  import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
  import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
  import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
  import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

  import java.io.StringReader;
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

      Tree[][] _gtSamples;
      double[][][] _gtLikelihoods;

      public void setRandomSeed(Long seed){
          _seed = seed;
      }

      public void setPrintDetails(boolean ifPrint){
          _printDetails = ifPrint;
      }

      public void setParallel(int numThreads){
          _numThreads = numThreads;
      }

      public FullLikelihoodFromSequence(List<Tree> gts, String[] gtTaxa, List<Tuple<char[],Integer>> sequences, SubstitutionModel model, double theta, double minHeight, double maxHeight){
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



      private Tree cloneTreeWithBranchLength(Tree tr, double[] nodeHeight, double scalor){
          Tree trCopy = new STITree(tr);
          int index = 0;
          Map<TNode, Double> node2height = new HashMap<TNode, Double>();
          for (TNode node : trCopy.postTraverse()) {
              double height = 0;
              if (!node.isLeaf()) {
                  height = nodeHeight[index++];
              }
              for (TNode child : node.getChildren()) {
                  double bl = height - node2height.get(child);
                  if(bl < 0){
                      return null;
                  }
                  child.setParentDistance(bl * scalor);
              }
              node2height.put(node, height);
          }
          return trCopy;
      }

      private double[] computeGTLikelhood(Tree gt){
          double[] results = new double[_sequences.size()];
          int index = 0;
          for(Tuple<char[], Integer> seq: _sequences) {
              Map<String, Character> sequenceMap = new HashMap<String, Character>();
              for(int i=0; i<_gtTaxa.length; i++){
                  sequenceMap.put(_gtTaxa[i], seq.Item1[i]);
              }
              Felsenstein fcalc = new Felsenstein(_model);
              results[index++] = fcalc.getLikelihoodtree(gt, new ObservationMap(sequenceMap));
          }
          return results;
      }

      private void clearParallelIndex(){
          _currentID = 0;
      }


      public void computeGTLikelhoodWithIntegral(int numSamplesPerBin, int numBins){
          _numBins = numBins;
          _numSamplesPerBin = numSamplesPerBin;
          _dim = _gts.get(0).getNodeCount() - _gts.get(0).getLeafCount();
          _numSamples = (int)Math.pow(_numBins, _dim) * numSamplesPerBin;
          _gtSamples = new Tree[_gts.size()][_numSamples];
          _gtLikelihoods = new double[_gts.size()][_numSamples][_sequences.size()];


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


      public double computeNetworkLikelhoodWithIntegral(Network network, Map<String, List<String>> species2alleles) {
          clearParallelIndex();
          Thread[] myThreads = new Thread[_numThreads];
          String networkExp = Networks.network2string(network);
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

                      if(_printDetails && _gtSamples[j][k]!=null){
                          System.out.println(_gtSamples[j][k] + " : " + _gtLikelihoods[j][k][i] + " * " + networkLikelihood[j][k]);
                      }

                  }
                  siteProb += MultivariateMonteCarloIntegral.integrate(_dim, _minHeight, _maxHeight, integralValues);
              }
              //System.out.println(Math.log(siteProb) + " * " + _sequences.get(i).Item2);
              //System.out.println(siteProb);
              finalResult += Math.log(siteProb) * _sequences.get(i).Item2;
          }

          //System.out.println(finalResult);
          return finalResult;
      }


     private void printArray(double array[]){
         for(double value: array){
             System.out.print(value + " ");
         }
         System.out.println();
     }

      private class GTLikelihoodThread extends Thread{

          public void run() {
              int myID;
              while((myID = getNextID(_gts.size()))!=-1){
                  Tree gt = _gts.get(myID);
                  //GTBranchLengthsValidator validator = new GTBranchLengthsValidator(gt, _dim);
                  MultivariateMonteCarloIntegral integration = new MultivariateMonteCarloIntegral(_numSamplesPerBin, _numBins);
                  integration.setRandomSeed(_seed);
                  int j = 0;
                  for (double[] sample : integration.getSamplePoints(_dim, _minHeight, _maxHeight)) {
                      //printArray(sample);
                      Tree gtCopy = cloneTreeWithBranchLength(gt, sample, _theta);
                      _gtSamples[myID][j] = gtCopy;
                      if(gtCopy == null){
                           Arrays.fill(_gtLikelihoods[myID][j], 0);
                      }
                      else {
                          _gtLikelihoods[myID][j] = computeGTLikelhood(gtCopy);
                          Trees.scaleBranchLengths((MutableTree) gtCopy, 1 / _theta);
                      }

                      j++;
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


      private class NetworkLikelihoodThread extends Thread{
                Network _network;
                Map<String, List<String>> _species2alleles;
                double[][] _results;

                public NetworkLikelihoodThread(String network, Map<String, List<String>> species2alleles, double[][] results){
                    _network = Networks.string2network(network);
                    _species2alleles = species2alleles;
                    _results = results;
                }


          public void run() {
              int myID;
              while((myID = getNextID(_gtSamples.length))!=-1) {
                  List<MutableTuple<Tree, Double>> gtList = new ArrayList<MutableTuple<Tree, Double>>();
                  List<Integer> originalIDs = new ArrayList<Integer>();
                  int index = 0;
                  for (Tree gt : _gtSamples[myID]) {
                      if(gt != null) {
                          gtList.add(new MutableTuple<Tree, Double>(gt, 1.0));
                          originalIDs.add(index);
                      }
                      index++;

                  }
                  double[] probs = new double[gtList.size()];
                  if(probs.length!=0) {
                      GeneTreeWithBranchLengthProbabilityYF gtp = new GeneTreeWithBranchLengthProbabilityYF(_network, gtList, _species2alleles);
                      gtp.calculateGTDistribution(probs);
                      Iterator<Integer> idIt = originalIDs.iterator();
                      for(double prob: probs){
                          _results[myID][idIt.next()] = prob;
                      }
                  }
                  //System.out.println(Networks.network2string(_network));
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



 /*
           class GTBranchLengthsValidator implements Func1<double[],Boolean>{
               Tree _gt;
               boolean[][] _R;

               public GTBranchLengthsValidator(Tree gt, int dim){
                   _gt = gt;
                   findGTNodesAncestralRelationship(dim);
               }


               public void findGTNodesAncestralRelationship(int numInternalNodes){
                   _R = new boolean[numInternalNodes][numInternalNodes];
                   Map<TNode, Integer> node2index = new HashMap<TNode, Integer>();
                   int index = 0;
                   for(TNode node: _gt.postTraverse()){
                       for(TNode child: node.getChildren()){
                           if(!child.isLeaf()) {
                               _R[index][node2index.get(child)] = true;
                           }
                       }
                       if(!node.isLeaf()) {
                           node2index.put(node, index++);
                       }
                   }
               }



               public Boolean execute(double[] argument) {
                   for (int i = 0; i < argument.length; i++) {
                       for (int j = 0; j < argument.length; j++) {
                           if (_R[i][j]) {
                               if (argument[i] <= argument[j]) {
                                   return false;
                               }
                           }
                       }
                   }


                   return true;
               }
           }
      */



  }
