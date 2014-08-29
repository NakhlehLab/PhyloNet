package edu.rice.cs.bioinfo.programs.phylonet.algos.integration;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.FelsensteinAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;


import java.io.StringReader;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Map;

/**
 * Created by yunyu on 6/3/14.
 */
public class GTBranchLengthsIntegrationForSequence {
    Tree _gt;
    double _minBL;
    double _maxBL;
    Map<String,Character> _leaf2char;
    SubstitutionModel _model;
    int _counter = 0;


    public GTBranchLengthsIntegrationForSequence(Tree gt, SubstitutionModel model, Map<String,Character> omap, double minbl, double maxbl){
        _gt = gt;
        _minBL = minbl;
        _maxBL = maxbl;
        _leaf2char = omap;
        _model = model;

    }

    public double computeLikelihoodWithIntegral(int sampleSize){
        return computeLikelihoodWithIntegral(sampleSize, 1);
    }

    public double computeLikelihoodWithIntegral(int sampleSize, int bins){
        MultivariateMonteCarloIntegral integration = new MultivariateMonteCarloIntegral(sampleSize, bins);
        FelsensteinProbabilityFunction function = new FelsensteinProbabilityFunction(_gt, _model, _leaf2char, _minBL, _maxBL);
        int dimension = function.getNumArguments();
        double[] mins = new double[dimension];
        Arrays.fill(mins, _minBL);
        double[] maxs = new double[dimension];
        Arrays.fill(maxs, _maxBL);
        double result = integration.integrate(function, dimension, mins, maxs);
        return result;
    }


    class FelsensteinProbabilityFunction implements Func1<double[], Double> {
        Tree _gt;
        int _numArguments;
        double _lowerBound;
        double _upperBound;
        Map<String,Character> _leaf2char;
        SubstitutionModel _model;

        public FelsensteinProbabilityFunction(Tree gt, SubstitutionModel model, Map<String,Character> omap, double lowerBound, double upperBound){
            _gt = gt;
            _numArguments = _gt.getNodeCount() - 1;
            _leaf2char = omap;
            _model = model;
            _lowerBound = lowerBound;
            _upperBound = upperBound;
        }

        public Double execute(double[] argument){
            if(_counter++ % 1000 == 0){
                System.out.println(_counter);
            }

            int index = 0;
            for(TNode node: _gt.postTraverse()){
                if(!node.isRoot())
                    node.setParentDistance(argument[index++]);
            }

            OneNucleotideObservation converter = new OneNucleotideObservation(_leaf2char);

            FelsensteinAlgorithm fcalc = new FelsensteinAlgorithm(_gt,_model);
            double result = fcalc.getProbability(converter);
            //System.out.println(_gt.toString()+": " + result);
            return result;
        }

        public int getNumArguments(){
            return _numArguments;
        }

    }

    public static void main(String[] args){
        double[] rates = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        double[] freqs = {0.25, 0.25, 0.25, 0.25};
        GTRModel gtrsm = new GTRModel(freqs,rates);
        NewickReader nr = new NewickReader(new StringReader("((SA:1.0,SB:1.0)SAB:0.5,(SC:1.0,SD:1.0)SCD:0.5)root"));
        STITree<Double> tree = new STITree<Double>(true);
        try {
            nr.readTree(tree);
        }
        catch(Exception e) {
            System.err.println(e);
            e.printStackTrace();
            return;
        }

        Map<String,Character> omap = new Hashtable<String,Character>();
        omap.put("SA", 'A');
        omap.put("SB", 'A');
        omap.put("SC", 'T');
        omap.put("SD", 'A');

        GTBranchLengthsIntegrationForSequence integrator = new GTBranchLengthsIntegrationForSequence(tree, gtrsm, omap, 0, 6);
        System.out.println(integrator.computeLikelihoodWithIntegral(10000));
    }
}
