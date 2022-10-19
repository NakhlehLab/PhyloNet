package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.all;
/*
 * @ClassName:   DeltaExchange
 * @Description: Migrated by zhen from BEAST
 * @Date:        8/17/21 8:09 PM
 */

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.Operator;

import java.util.*;


public class DeltaExchange extends Operator {
    private double _logHR;
    protected boolean _violate;

    private double _delta = Utils._DELTA;
    private boolean _isIntegerOperator = Utils._IS_INTEGER_OPERATOR;
    private List<UltrametricTree> _alltrees = null;
    public List<Double> _parameter_input = Utils._PARAMETER_INPUT;
    public List<Integer> _weight_vector = Utils._WEIGHT_VECTOR;
    private double _lowerbound = Utils._LOWERBOUND;
    private double _upperbound = Double.MAX_VALUE;
    public double _scalar1;
    public double _scalar2;
    private double _prev_scalar1;
    private double _prev_scalar2;
    public int _dim1;
    public int _dim2;


    /* Constructor */
    public DeltaExchange(List<UltrametricTree> trees) {
        _alltrees = trees;
        if (!_parameter_input.isEmpty()){
            if (_parameter_input.size() != _alltrees.size()){
                throw new IllegalArgumentException(
                        "number of parameters "+_parameter_input.size()+" should have the same length as number of trees "+_alltrees.size());
            }
            for(int i = 0; i < _parameter_input.size(); i++){
                UltrametricTree ut = _alltrees.get(i);
                ut.updateMutationRate(_parameter_input.get(i));
            }
        }
        else{
            for(int i = 0; i < _alltrees.size(); i++){
                _parameter_input.add(_alltrees.get(i).get_mutationRate());
            }
        }
    }

    public DeltaExchange() {


    }

    private int[] weights() {
        int[] weights;
        weights = new int[_parameter_input.size()];
        if (!_weight_vector.isEmpty()){
            if (_weight_vector.size() != weights.length){
                throw new IllegalArgumentException(
                        "Weights vector should have the same length as parameter dimension");
            }
            for(int i = 0; i < weights.length; i++){
                weights[i] = _weight_vector.get(i);
            }
        }
        else{
            for(int i = 0; i < weights.length; i++){
                weights[i] = 1;
            }
        }
        return weights;
    }

    @Override
    public double propose() {
        int[] parameterWeights = weights();
        final int dim = parameterWeights.length;

        // Find the number of weights that are nonzero
        int nonZeroWeights = 0;
        for (int i: parameterWeights) {
            if (i != 0) {
                ++nonZeroWeights;
            }
        }

        if (nonZeroWeights <= 1) {
            // it is impossible to select two distinct entries in this case, so there is nothing to propose
            return 0.0;
        }

        // Generate indices for the values to be modified
        _dim1 = Randomizer.getRandomInt(nonZeroWeights);
        _dim2 = Randomizer.getRandomInt(nonZeroWeights-1);
        if (_dim2 >= _dim1) {
            ++_dim2;
        }
        if (nonZeroWeights<dim) {
            // There are zero weights, so we need to increase dim1 and dim2 accordingly.
            int nonZerosBeforeDim1 = _dim1;
            int nonZerosBeforeDim2 = _dim2;
            _dim1 = 0;
            _dim2 = 0;
            while (nonZerosBeforeDim1 > 0 | parameterWeights[_dim1] == 0 ) {
                if (parameterWeights[_dim1] != 0) {
                    --nonZerosBeforeDim1;
                }
                ++_dim1;
            }
            while (nonZerosBeforeDim2 > 0 | parameterWeights[_dim2] == 0 ) {
                if (parameterWeights[_dim2] != 0) {
                    --nonZerosBeforeDim2;
                }
                ++_dim2;
            }
        }

        _logHR = 0.0;


        _scalar1 =  _parameter_input.get(_dim1);
        _scalar2 =  _parameter_input.get(_dim2);


        if (_isIntegerOperator) {
            final int d = Randomizer.getRandomInt((int) Math.round(_delta)) + 1;

            if (parameterWeights[_dim1] != parameterWeights[_dim2]) throw new RuntimeException();
            _scalar1 = Math.round(_scalar1 - d);
            _scalar2 = Math.round(_scalar2 + d);
        }
        else {


            // exchange a random delta
            final double d = Randomizer.getRandomDouble() * _delta;
            _scalar1 -= d;
            if (parameterWeights[_dim1] != parameterWeights[_dim2]) {
                _scalar2 += d * parameterWeights[_dim1] / parameterWeights[_dim2];
            }
            else {
                _scalar2 += d;
            }
        }
        if(Utils.DEBUG_MODE){
            System.out.println(_dim1 +","+ _scalar1);
            System.out.println(_dim2 +","+ _scalar2);
        }

        if(_scalar1 < _lowerbound || _scalar1 > _upperbound
                || _scalar2 < _lowerbound || _scalar2 > _upperbound){
            _logHR = Double.NEGATIVE_INFINITY;
        }
        else{
            _prev_scalar1 = _parameter_input.get(_dim1);
            _prev_scalar2 = _parameter_input.get(_dim2);
            _parameter_input.set(_dim1, _scalar1);
            _parameter_input.set(_dim2, _scalar2);
            if (_alltrees != null){

                UltrametricTree ut1 = _alltrees.get(_dim1);
//                ut1.scale(_scalar1/_prev_scalar1, false);
                _alltrees.get(_dim1).updateMutationRate(_scalar1);
                ut1.setDirty(true);

                UltrametricTree ut2 = _alltrees.get(_dim2);
//                ut2.scale(_scalar2/_prev_scalar2, false);
                _alltrees.get(_dim2).updateMutationRate(_scalar2);
                ut2.setDirty(true);

                if (Utils.DEBUG_MODE) {
                    System.out.println("Delta Exchange: "+_dim1 + ", "+ _dim2);
                    System.out.println(_parameter_input);
                }

            }
        }

        // symmetrical move so return a zero hasting ratio
        return _logHR;
    }

    @Override
    public void undo() {
        if (_logHR != Double.NEGATIVE_INFINITY){
            _parameter_input.set(_dim1, _prev_scalar1);
            _parameter_input.set(_dim2, _prev_scalar2);

            UltrametricTree ut1 = _alltrees.get(_dim1);
//            ut1.scale(_prev_scalar1/_scalar1, true);
            _alltrees.get(_dim1).updateMutationRate(_prev_scalar1);
            ut1.setDirty(true);

            UltrametricTree ut2 = _alltrees.get(_dim2);
//            ut2.scale(_prev_scalar2/_scalar2, true);
            _alltrees.get(_dim2).updateMutationRate(_prev_scalar2);
            ut2.setDirty(true);
        }
    }

    @Override
    public String getName() {
        return "DeltaExchange";
    }

    @Override
    public void optimize(double logAlpha) {
        double delta = Utils.calcDelta(this, logAlpha);
        delta += Math.log(_delta);
        _delta = Math.exp(delta);
        if (_isIntegerOperator) {
            // when delta < 0.5
            // Randomizer.nextInt((int) Math.round(delta)) becomes
            // Randomizer.nextInt(0) which results in an exception
            _delta = Math.max(0.5000000001, _delta);
        }
//        System.out.println("delta: "+_delta);
    }

    public double getCoercableParameterValue() {
        return _delta;
    }

    public void setCoercableParameterValue(final double value) {
        _delta = value;
    }

    public Utils.MOVE_TYPE getCategory() {
        return Utils.MOVE_TYPE.ALL;
    }

    public boolean mayViolate() {
        return _violate;
    }
}
