package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.prior;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.PopulationSize;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;

/**
 * Created by dw20 on 8/19/16.
 * Scale the hyper parameter of the Gamma distribution for population size
 */
public class ScalePopSizeParam extends Operator{

    private double _scaleFactor = 0.92;
    private double _upperLimit = 1.0 - 1e-6;
    private double _lowerLimit = 0.8;

    private PopulationSize _popSize;
    private Double scale;
    private double _prevParam = -1;

    public ScalePopSizeParam(PopulationSize ps) {
        this._popSize = ps;
    }

    @Override
    public double propose() {
        scale = getScaler();
        _prevParam = _popSize.getGammaMean();
        _popSize.setGammaMean(scale * _prevParam);
        return -Math.log(scale);
    }

    @Override
    public void undo() {
        _popSize.setGammaMean(_prevParam);
        _prevParam = -1;
    }

    @Override
    public void optimize(double logAlpha) {
        double delta = Utils.calcDelta(this, logAlpha);
        delta += Math.log(1.0 / _scaleFactor - 1.0);
        setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
    }

    @Override
    public Utils.MOVE_TYPE getCategory() {
        return Utils.MOVE_TYPE.PRIOR;
    }

    @Override
    public boolean mayViolate() {
        return false;
    }

    @Override
    public String getName() {
        return "Scale-PopSize-Prior-Param";
    }

    private void setCoercableParameterValue(final double value) {
        _scaleFactor = Math.max(Math.min(value, _upperLimit), _lowerLimit);
    }

    private double getScaler() {
        return _scaleFactor + Randomizer.getRandomDouble() * (1.0 / _scaleFactor - _scaleFactor);
    }

    public static void main(String[] args) {
        {
            Utils._CONST_POP_SIZE = true;
            Utils._ESTIMATE_POP_SIZE = true;
            UltrametricNetwork net = new UltrametricNetwork("(O:0.5,((((A:0.5,(C:0.5,(G:0.424733668924115)I10#H3:0.27526633107588505::0.95)I6:1.2148953368828783)I5:2.0902392667517367,(((R:0.23468712)I8#H2:0.57719209683010764::0.76)I2#H1:0.4723366810412248::0.29,(Q:0.5,I10#H3:0.24368001872495645::0.05)I4:1.4561500003545642)I9:2.081453711192161)I3:0.8324858688602822,(L:0.30183286,I8#H2:0.8699237693946538::0.24)I7:0.7369103845180647)I1:1.427856628628852,I2#H1:0.8195085232672825::0.71)I0:3.0);");
            int runs = 100000;
            int counter = 0;
            PopulationSize ps = new PopulationSize();
            System.out.println(ps.getGammaMean());
            for(int i = 0; i < runs; i++) {
                Operator op = new ScalePopSizeParam(ps);
                double logHR = op.propose();
                if(logHR != Utils.INVALID_MOVE) counter++;
                if(logHR < Math.log(Randomizer.getRandomDouble())) op.undo();
            }
            System.out.println(counter == runs);
            System.out.printf("%d out of %d\n", counter, runs);
            System.out.println(ps.getGammaMean());
        }
        {
            UltrametricNetwork net = new UltrametricNetwork("(O:0.5,((((A:0.5,(C:0.5,(G:0.424733668924115)I10#H3:0.27526633107588505::0.95)I6:1.2148953368828783)I5:2.0902392667517367,(((R:0.23468712)I8#H2:0.57719209683010764::0.76)I2#H1:0.4723366810412248::0.29,(Q:0.5,I10#H3:0.24368001872495645::0.05)I4:1.4561500003545642)I9:2.081453711192161)I3:0.8324858688602822,(L:0.30183286,I8#H2:0.8699237693946538::0.24)I7:0.7369103845180647)I1:1.427856628628852,I2#H1:0.8195085232672825::0.71)I0:3.0);");
            int runs = 10000;
            int test = 0;
            PopulationSize ps = new PopulationSize();
            System.out.println(ps.getGammaMean());
            double original = ps.getGammaMean();
            for(int i = 0; i < runs; i++) {
                Operator op = new ScalePopSizeParam(ps);
                double logHR = op.propose();
                op.undo();
                if(ps.getGammaMean() != original) {
                    test--;
                    System.err.println(ps.getGammaMean());
                }
                test++;
            }
            System.out.println(test == runs);
            System.out.printf("%d out of %d\n", test, runs);
            System.out.println(ps.getGammaMean());
        }
    }
}