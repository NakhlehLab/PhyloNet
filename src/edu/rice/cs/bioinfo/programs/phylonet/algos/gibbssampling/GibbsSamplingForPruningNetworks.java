package edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling;

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling.util.FindMaximumForEnvelopeUsingBrent;
import edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling.util.FindMaximumForEnvelopeUsingMinuit;
import edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling.util.PruneNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import org.apache.commons.math3.distribution.*;


import java.awt.geom.Point2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.StringWriter;
import java.util.*;

/**
 * Created by yunyu on 10/1/15.
 */
//public abstract class GibbsSamplingForPruningNetworks extends GibbsSamplingBase<double[], Tuple<List<List<MutableTuple<Tree,Double>>>,Map<String, List<String>>>> {

public abstract class GibbsSamplingForPruningNetworks{
    private double _epsilon = 0.001;
    private boolean _printDetails = false;
    protected double[] _branchLengthBound = {0.001, 6};
    protected double[] _inheritanceProbBound = {1.0/Math.pow(10, 6), 1- 1.0/Math.pow(10, 6)};
    protected int _numIterations = 11000;
    protected int _burnIn = 1000;
    protected int _sampleInterval = 100;
    protected int _maxFailureForRejectSampling = 10;
    protected List<Tuple<NetNode,NetNode>> _parametersToSample;
    protected int _numProcessors = 16;
    protected double _pruneThreshold = 0.01;
    protected File _outputResultFile;
    protected File _outputProbFile;
    protected StringWriter _outputString;

    private int _startIteration = 1;

    public void setStartIteration(int start){
        _startIteration = start;
    }

    public void setNumProcessors(int numProcessors){
        _numProcessors = numProcessors;
    }

    public void setOutputResultFile(String fileName){
        _outputResultFile = new File(fileName);
        if(_outputResultFile.exists() && _startIteration==1){
            _outputResultFile.delete();
        }
    }

    public void setOutputString(StringWriter writer){
        _outputString = writer;
    }

    public void setParameters(double maxBranchLength, int numIterations, int burnIn, int sampleInterval, double pruneThreshold, int numProcessor){
        _branchLengthBound[1] = maxBranchLength;
        _numIterations = numIterations;
        _burnIn = burnIn;
        _sampleInterval = sampleInterval;
        _pruneThreshold = pruneThreshold;
        _numProcessors = numProcessor;
    }

    public void setOutputProbFile(String fileName){
        _outputProbFile = new File(fileName);
    }

    public List<MutableTuple<Network,Integer>> sample(Network network, List<List<MutableTuple<Tree,Double>>> geneTrees, Map<String, List<String>> species2alleles){
        Networks.autoLabelNodes(network);
        Map<String,String> allele2species = null;
        if(species2alleles!=null){
            allele2species = new HashMap<String, String>();
            for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                for(String allele: entry.getValue()){
                    allele2species.put(allele, entry.getKey());
                }
            }
        }

        List summarizedData = new ArrayList<>();
        List dataCorrespondences = new ArrayList<>();
        summarizeGTs(geneTrees, allele2species, summarizedData, dataCorrespondences);

        _parametersToSample = getAllParametersToSample(network, summarizedData, allele2species);
        if(_printDetails){
            System.out.println("Number of parameters: " + _parametersToSample.size() + "\n");
        }
        double totalLnPrior;
        if(_startIteration == 1){
            totalLnPrior = generateInitialSample();
        }
        else{
            totalLnPrior = getInitialSample();
        }
        if(_printDetails){
            System.out.println("Initial network: " + network.toString());
        }

        double lnLikelihood = computeInitialLikelihood(network, summarizedData, species2alleles, dataCorrespondences);

        if(_printDetails){
            System.out.println("TotalLnPrior: " + totalLnPrior + "  LnLikelihood: " + lnLikelihood + "\n");
        }

        MutableTuple<Double,Double> lnLikelihoodPrior = new MutableTuple<>(lnLikelihood, totalLnPrior);
        List<Double> lnPosteriorProbPlot = new ArrayList<>();
        List<MutableTuple<Network,Integer>> results = new ArrayList<>();
        for(int i=_startIteration; i<=_numIterations; i++){
            Collections.shuffle(_parametersToSample);

            if(_printDetails){
                System.out.println("\n\nIteration #" + i + ":");
            }
            double[] newSample = getNextSample(network, summarizedData, species2alleles, dataCorrespondences, lnLikelihoodPrior);
            if(_printDetails) {
                System.out.println("Network: " + network.toString());
                System.out.println("TotalLnPrior: " + lnLikelihoodPrior.Item2 + "  LnLikelihood: " + lnLikelihoodPrior.Item1 + "\n\n");
            }
            lnPosteriorProbPlot.add(lnLikelihoodPrior.Item1 + lnLikelihoodPrior.Item2);

            if(i > _burnIn && i % _sampleInterval == 0){
                Network copy = Networks.readNetwork(network.toString());
                PruneNetwork.prune(copy, _pruneThreshold);

                if(_outputString != null){
                    _outputString.append("Iteration #" + i + ":\n");
                    _outputString.append("Network: " + network.toString() + "\n");
                    _outputString.append("Posterior: " + (lnLikelihoodPrior.Item1 + lnLikelihoodPrior.Item2) + "\n\n");
                }

                if(_outputResultFile != null) {
                    try {
                        BufferedWriter bw = new BufferedWriter(new FileWriter(_outputResultFile, true));
                        bw.append("#" + i + ":\n");
                        bw.append(network.toString()+"\n");
                        bw.append(copy.toString()+"\n");
                        bw.append((lnLikelihoodPrior.Item1 + lnLikelihoodPrior.Item2) + "");
                        bw.append("\n\n");
                        bw.close();
                    } catch (Exception e) {
                    }
                }

                if(_outputProbFile != null) {
                    try {
                        BufferedWriter bw = new BufferedWriter(new FileWriter(_outputProbFile, true));
                        for (double prob : lnPosteriorProbPlot) {
                            bw.append(prob + ", ");
                        }
                        bw.close();
                        lnPosteriorProbPlot.clear();
                    } catch (Exception e) {
                    }
                }


                boolean exist = false;
                for(MutableTuple<Network,Integer> tuple: results){
                    if(Networks.hasTheSameTopology(tuple.Item1, copy)){
                        tuple.Item2++;
                        exist = true;
                        break;
                    }
                }
                if(!exist){
                    results.add(new MutableTuple<>(copy, 1));
                }
                if(_printDetails){
                    for(MutableTuple<Network,Integer> tuple: results){
                        System.out.println(tuple.Item2 + ": " + tuple.Item1.toString());
                    }
                }
            }

        }


        for(int i=1;i<results.size();i++){
            MutableTuple<Network, Integer> tuple1 = results.get(i);
            for(int j=0;j<i;j++){
                MutableTuple<Network, Integer> tuple2 = results.get(j);
                if(tuple2.Item2 < tuple1.Item2){
                    results.remove(tuple1);
                    results.add(j,tuple1);
                    break;
                }
            }
        }

        if(_outputResultFile != null) {
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(_outputResultFile, true));
                bw.append("\n\nFinal result:\n");
                for (MutableTuple<Network, Integer> tuple : results) {
                    bw.append(tuple.Item2 + ": " + tuple.Item1.toString() + "\n");
                }
                bw.close();
            } catch (Exception e) {
            }
        }

        if(_outputString != null){
            _outputString.append("\nFinal result after pruning:\n");
            double totalSamples = (_numIterations - _burnIn) / _sampleInterval;
            for (MutableTuple<Network, Integer> tuple : results) {
                _outputString.append(tuple.Item2 / totalSamples + ": " + tuple.Item1.toString() + "\n");
            }
        }

        return results;
    }

    protected double generateInitialSample(){
        double totalLnProb = 0;
        for(Tuple<NetNode,NetNode> edge: _parametersToSample){
            if(edge.Item1 != null){
                ExponentialDistribution distribution = new ExponentialDistribution(1);
                double branchLength;
                do{
                    branchLength = distribution.sample();
                }while (branchLength < _branchLengthBound[0] || branchLength > _branchLengthBound[1]);
                double prob = distribution.density(branchLength);
                totalLnProb += Math.log(prob);
                edge.Item2.setParentDistance(edge.Item1, branchLength);
            }
            else{
                BetaDistribution distribution = new BetaDistribution(0.1, 0.1);
                double inheritanceProb;
                do{
                    //inheritanceProb = distribution.sample();
                    inheritanceProb = 0.5;
                }while (inheritanceProb < _inheritanceProbBound[0] || inheritanceProb > _inheritanceProbBound[1]);
                double prob = distribution.density(inheritanceProb);
                totalLnProb += Math.log(prob);
                NetNode reticulationNode = edge.Item2;
                for(Object parentNode: reticulationNode.getParents()){
                    reticulationNode.setParentProbability((NetNode)parentNode, inheritanceProb);
                    inheritanceProb = 1 - inheritanceProb;

                }
            }
        }
        return totalLnProb;
    }

    protected double getInitialSample(){
        double totalLnProb = 0;
        for(Tuple<NetNode,NetNode> edge: _parametersToSample){
            if(edge.Item1 != null){
                double branchLength = edge.Item2.getParentDistance(edge.Item1);
                ExponentialDistribution distribution = new ExponentialDistribution(1);
                double prob = distribution.density(branchLength);
                totalLnProb += Math.log(prob);
            }
            else{
                double inheritanceProb = edge.Item2.getParentProbability((NetNode)edge.Item2.getParents().iterator().next());
                BetaDistribution distribution = new BetaDistribution(0.1, 0.1);
                double prob = distribution.density(inheritanceProb);
                totalLnProb += Math.log(prob);
            }
        }
        return totalLnProb;
    };


    protected  double[] getNextSample(Network network, List summarizedData, Map<String, List<String>> species2alleles, List dataCorrespondences, MutableTuple<Double,Double> lnLikelihoodPrior){
        double[] newSample = new double[_parametersToSample.size()];
        int index = 0;
        for(Tuple<NetNode,NetNode> edge: _parametersToSample){
            if(_printDetails){
                if (edge.Item1 != null) {
                    System.out.println("Sample branch length of (" + edge.Item1.getName() + "," + edge.Item2.getName() + ")");
                }else{
                    System.out.println("Sample inheritance probability of node " + edge.Item2.getName());
                }
            }

            double originalValue;
            if (edge.Item1 != null) {
                originalValue = edge.Item2.getParentDistance(edge.Item1);
            } else {
                originalValue = edge.Item2.getParentProbability((NetNode)edge.Item2.getParents().iterator().next());
            }
            double totalLnPriorMinus = lnLikelihoodPrior.Item2 - computeLnPrior(originalValue, edge.Item1 != null);
            if(Double.isNaN(totalLnPriorMinus)){
                System.out.println("Nan");
            }

            Func<Double> lnLikelihoodFunc = likelihoodUpdater(network, summarizedData, species2alleles, dataCorrespondences, edge);

            Point2D maxPoint = findMaximumValue(edge, lnLikelihoodFunc, totalLnPriorMinus);
            double[][] envelopes = new double[1][3];
            envelopes[0][2] = maxPoint.getY();
            double[] boundary = edge.Item1 != null? _branchLengthBound: _inheritanceProbBound;
            envelopes[0][0] = boundary[0];
            envelopes[0][1] = boundary[1];
            Point2D maxPoint1 = null, maxPoint2 = null;
            if(maxPoint.getX() - boundary[0] < _epsilon){
                maxPoint2 = new Point2D.Double(maxPoint.getX(), maxPoint.getY());
            }
            else if(boundary[1] - maxPoint.getX() < _epsilon){
                maxPoint1 = new Point2D.Double(maxPoint.getX(), maxPoint.getY());
            }
            else{
                maxPoint1 = new Point2D.Double(maxPoint.getX(), maxPoint.getY());
                maxPoint2 = new Point2D.Double(maxPoint.getX(), maxPoint.getY());
            }


            double newSamplePoint;
            do{
                newSamplePoint = doRejectionSampling(edge, envelopes, lnLikelihoodFunc, totalLnPriorMinus, lnLikelihoodPrior, _maxFailureForRejectSampling);

                if(newSamplePoint == -1){
                    envelopes = buildStairEnvelope(edge, envelopes, lnLikelihoodFunc, maxPoint1, maxPoint2, totalLnPriorMinus);
                }

            }while(newSamplePoint == -1);

            newSample[index++] = newSamplePoint;
        }

        return newSample;
    }


    private boolean checkEnvelopes(double[][] envelopes, double[] boundary){
        for(int i=0; i<envelopes.length; i++){
            if(envelopes[i][0] < boundary[0] || envelopes[i][1] > boundary[1] || envelopes[i][0] == envelopes[i][1]){
                return false;
            }
        }
        return true;
    }


    private double[][] buildStairEnvelope(Tuple<NetNode,NetNode> edge, double[][] originalEnvelopes, Func<Double> lnLikelihoodFunc, Point2D maxPoint1, Point2D maxPoint2, double totalLnPriorMinus){
        if(maxPoint1 == null){
            return addStairEnvelopeAtRightSide(edge, originalEnvelopes, lnLikelihoodFunc, maxPoint2, totalLnPriorMinus);
        }
        else if(maxPoint2 == null){
            return addStairEnvelopeAtLeftSide(edge, originalEnvelopes, lnLikelihoodFunc, maxPoint1, totalLnPriorMinus);
        }
        else{
            return addStairEnvelopeAtBothSide(edge, originalEnvelopes, lnLikelihoodFunc, maxPoint1, maxPoint2, totalLnPriorMinus);
        }
    }


    private double[][] addStairEnvelopeAtRightSide(Tuple<NetNode,NetNode> edge, double[][] originalEnvelopes, Func<Double> lnLikelihoodFunc, Point2D maxPoint, double totalLnPriorMinus){
        double delta = 0.1;
        boolean findDelta = false;
        double subMaximum = -1;
        double addedX;
        do{
            addedX = maxPoint.getX() + delta;
            if (edge.Item1 != null) {
                if(addedX >= _branchLengthBound[1]){
                    delta = delta/2;
                    continue;
                }
                edge.Item2.setParentDistance(edge.Item1, addedX);
            } else {
                if(addedX >= _inheritanceProbBound[1]){
                    delta = delta/2;
                    continue;
                }
                for (Object parentNode : edge.Item2.getParents()) {
                    edge.Item2.setParentProbability((NetNode) parentNode, addedX);
                    addedX = 1 - addedX;
                }
            }
            double lnLikelihood = lnLikelihoodFunc.execute();
            subMaximum = lnLikelihood + totalLnPriorMinus + computeLnPrior(addedX, edge.Item1 != null);
            if(Math.exp(subMaximum - maxPoint.getY()) < 1.0/100){
                delta = delta/2;
            }else{
                findDelta = true;
                if(_printDetails){
                    System.out.println("Find delta " + delta);
                }
            }
        }while(!findDelta);

        double[][] newEnvelopes = new double[originalEnvelopes.length+1][3];
        for(int i=0; i<originalEnvelopes.length; i++){
            newEnvelopes[i] = Arrays.copyOf(originalEnvelopes[i], 3);
        }
        newEnvelopes[originalEnvelopes.length][1] = newEnvelopes[originalEnvelopes.length-1][1];
        newEnvelopes[originalEnvelopes.length-1][1] = addedX;
        newEnvelopes[originalEnvelopes.length][0] = addedX;
        newEnvelopes[originalEnvelopes.length][2] = subMaximum;


        if(_printDetails) {
            System.out.println("Find subMaximum " + addedX + " resulting in " + subMaximum);
        }
        if(!checkEnvelopes(newEnvelopes, edge.Item1 != null?_branchLengthBound:_inheritanceProbBound)){
            throw new RuntimeException("Wrong");
        }
        maxPoint.setLocation(addedX, subMaximum);
        return newEnvelopes;
    }


    private double[][] addStairEnvelopeAtLeftSide(Tuple<NetNode,NetNode> edge, double[][] originalEnvelopes, Func<Double> lnLikelihoodFunc, Point2D maxPoint, double totalLnPriorMinus){
        double delta = 0.1;
        boolean findEpsilon = false;
        double subMaximum = -1;
        double addedX;
        do{
            addedX = maxPoint.getX() - delta;
            if (edge.Item1 != null) {
                if(addedX <= _branchLengthBound[0]){
                    delta = delta/2;
                    continue;
                }
                edge.Item2.setParentDistance(edge.Item1, addedX);
            } else {
                if(addedX <= _inheritanceProbBound[0]){
                    delta = delta/2;
                    continue;
                }
                for (Object parentNode : edge.Item2.getParents()) {
                    edge.Item2.setParentProbability((NetNode) parentNode, addedX);
                    addedX = 1 - addedX;
                }
            }
            double lnLikelihood = lnLikelihoodFunc.execute();
            subMaximum = lnLikelihood + totalLnPriorMinus + computeLnPrior(addedX, edge.Item1 != null);
            if(Math.exp(subMaximum - maxPoint.getY()) < 1.0/100){
                delta = delta/2;
            }else{
                findEpsilon = true;
                if(_printDetails){
                    System.out.println("Find delta " + delta);
                }
            }
        }while(!findEpsilon);

        double[][] newEnvelopes = new double[originalEnvelopes.length+1][3];
        for(int i=0; i<originalEnvelopes.length; i++){
            newEnvelopes[i+1] = Arrays.copyOf(originalEnvelopes[i], 3);
        }
        newEnvelopes[0][0] = originalEnvelopes[0][0];
        newEnvelopes[0][1] = addedX;
        newEnvelopes[0][2] = subMaximum;
        newEnvelopes[1][0] = addedX;


        if(_printDetails) {
            System.out.println("Find subMaximum " + addedX + " resulting in " + subMaximum);
        }
        if(!checkEnvelopes(newEnvelopes, edge.Item1 != null?_branchLengthBound:_inheritanceProbBound)){
            throw new RuntimeException("Wrong");
        }
        maxPoint.setLocation(addedX, subMaximum);
        return newEnvelopes;
    }


    private double[][] addStairEnvelopeAtBothSide(Tuple<NetNode,NetNode> edge, double[][] originalEnvelopes, Func<Double> lnLikelihoodFunc, Point2D maxPoint1, Point2D maxPoint2, double totalLnPriorMinus){
        double[][] newEnvelopes = addStairEnvelopeAtLeftSide(edge, originalEnvelopes, lnLikelihoodFunc, maxPoint1, totalLnPriorMinus);
        newEnvelopes = addStairEnvelopeAtRightSide(edge, newEnvelopes, lnLikelihoodFunc, maxPoint2, totalLnPriorMinus);
        return newEnvelopes;
    }


    private Point2D findMaximumValue(Tuple<NetNode,NetNode> edge, Func<Double> lnLikelihoodFunc, double totalLnPriorMinus){
        Func2<Double, Boolean, Double> lnPriorFunc = new Func2<Double, Boolean, Double>() {
            @Override
            public Double execute(Double value, Boolean isBranchLength) {
                return computeLnPrior(value, isBranchLength);
            }
        };
        //find the maximum value of ln(p(theta,N)p(G|theta,N))
        //FindMaximumForEnvelopeUsingBrent findMaximum = new FindMaximumForEnvelopeUsingBrent();
        FindMaximumForEnvelopeUsingMinuit findMaximum = new FindMaximumForEnvelopeUsingMinuit();
        findMaximum.setBounds(_branchLengthBound, _inheritanceProbBound);
        Point2D maxPoint = findMaximum.getMaximum(edge, lnLikelihoodFunc, lnPriorFunc, totalLnPriorMinus);
        return maxPoint;
    }


    private double doRejectionSampling(Tuple<NetNode,NetNode> edge, double[][] envelopes, Func<Double> lnLikelihoodFunc, double totalLnPriorMinus, MutableTuple<Double,Double> lnLikelihoodPrior, int maxRejection){
        double[] weights = new double[envelopes.length];

        for(int i=0; i<envelopes.length; i++){
            double denominator = 0;
            double lnArea = envelopes[i][2] + Math.log(envelopes[i][1] - envelopes[i][0]);
            for(int j=0; j<envelopes.length; j++){
                if(i==j){
                    denominator += 1;
                }else{
                    denominator += Math.exp(envelopes[j][2] + Math.log(envelopes[j][1] - envelopes[j][0]) - lnArea);
                }
            }
            weights[i] = 1/denominator;
            if(i!=0){
                weights[i] += weights[i-1];
            }
        }

        if(Math.abs(weights[weights.length-1]-1)>0.0000001){
            throw new RuntimeException("Wrong weights");
        }

        double resultSample = -1;
        int rejectCounter = 0;
        do{
            double randomForArea = Math.random();
            int areaIndex = 0;
            for(; areaIndex<weights.length; areaIndex++){
                if(randomForArea <= weights[areaIndex]){
                    break;
                }
            }

            if(areaIndex>=envelopes.length){
                System.out.println("Wrong");
            }
            double randomSample = Math.random() * (envelopes[areaIndex][1]-envelopes[areaIndex][0]) + envelopes[areaIndex][0];
            if (edge.Item1 != null) {
                edge.Item2.setParentDistance(edge.Item1, randomSample);
            } else {
                for(Object parentNode: edge.Item2.getParents()){
                    edge.Item2.setParentProbability((NetNode)parentNode, randomSample);
                    randomSample = 1 - randomSample;
                }
            }
            double lnLikelihood = lnLikelihoodFunc.execute();
            double lnPrior = totalLnPriorMinus + computeLnPrior(randomSample, edge.Item1 != null);
            double sampleLnProb = lnLikelihood + lnPrior;
            double randomUniformLog = Math.log(Math.random());
            if(randomUniformLog < (sampleLnProb - envelopes[areaIndex][2])){
                lnLikelihoodPrior.Item1 = lnLikelihood;
                lnLikelihoodPrior.Item2 = lnPrior;
                resultSample = randomSample;
            }
            else{
                rejectCounter++;
            }
            if(_printDetails) {
                if (resultSample!=-1) {
                    System.out.println("Accept " + randomSample + ": " + Math.exp(randomUniformLog) + " < " + sampleLnProb + "/" + envelopes[areaIndex][2]);
                } else {
                    System.out.println("Reject " + randomSample + ": " + Math.exp(randomUniformLog) + " > " + sampleLnProb + "/" + envelopes[areaIndex][2]);
                }
            }
        }while(resultSample==-1 && rejectCounter<maxRejection);


        return resultSample;
    }


    protected double computeLnPrior(double value, boolean isBranchLength){
        AbstractRealDistribution distribution;
        if (isBranchLength) {
            distribution = new ExponentialDistribution(1);
        } else {
            distribution = new BetaDistribution(0.1, 0.1);
        }
        return Math.log(distribution.density(value));
    }

    abstract protected double computeInitialLikelihood(Network network, List summarizedData, Map<String, List<String>> species2alleles, List dataCorrespondences);

    abstract protected void summarizeGTs(List<List<MutableTuple<Tree,Double>>> originalGTs, Map<String, String> allele2species, List summarizedData, List dataCorrespondences);

    abstract protected List<Tuple<NetNode,NetNode>> getAllParametersToSample(Network network, List<Tree> gts, Map<String, String> allele2species);

    abstract protected Func<Double> likelihoodUpdater(Network network, List summarizedData, Map<String, List<String>> species2alleles,List dataCorrespondences, Tuple<NetNode, NetNode> editedEdge);

}
