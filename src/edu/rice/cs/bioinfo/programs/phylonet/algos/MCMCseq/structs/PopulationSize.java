package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.prior.ChangePopSizeParam;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.prior.ScalePopSizeParam;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import org.apache.commons.math3.distribution.GammaDistribution;

/**
 * Prior of the population sizes where hyper parameter is a state node that might be changed.
 * Created by dw20 on 8/19/16.
 */
public class PopulationSize extends StateNode {

    private GammaDistribution _popSize;
    private double _gammaMean = Utils._POP_SIZE_MEAN / 2;
    private Operator[] _operators;
    private double[] _opWeights;

    public PopulationSize() {
        _popSize = new GammaDistribution(Utils.GAMMA_SHAPE, _gammaMean);
        _operators = new Operator[]{new ChangePopSizeParam(this), new ScalePopSizeParam(this)};
        _opWeights = Utils.PopSize_Op_Weights;
    }

    public double getGammaMean() {
        return _gammaMean;
    }

    public void setGammaMean(double mean) {
        if(!Utils._ESTIMATE_POP_SIZE) {
            throw new RuntimeException("Don't allow to change gamma mean parameter for population size estimation");
        }
        _gammaMean = mean;
        _popSize = new GammaDistribution(Utils.GAMMA_SHAPE, _gammaMean);
    }

    @Override
    public double propose() {
        this._operator = getOp(_operators, _opWeights);
        return this._operator.propose();
    }

    @Override
    public void undo() {
        if(this._operator == null) throw new IllegalArgumentException("null operator");
        this._operator.undo();
    }

    @Override
    public void accept() {

    }

    @Override
    public void reject() {

    }

    @Override
    public double logDensity() {
        return -Math.log(_gammaMean);
    }

    @Override
    public boolean mayViolate() {
        return false;
    }

    @Override
    public boolean isValid() {
        return true;
    }

    public String toString() {
        return Double.toString(_gammaMean);
    }

    public double density(double popSize) {
        return _popSize.density(popSize);
    }


}
