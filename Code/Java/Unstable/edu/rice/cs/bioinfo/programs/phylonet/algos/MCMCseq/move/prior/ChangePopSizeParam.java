package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.prior;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.PopulationSize;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;

/**
 * Created by dw20 on 8/19/16.
 * Change the hyper parameter of the Gamma distribution for population size
 */
public class ChangePopSizeParam extends Operator {

    private double _windowSize = Utils._POP_SIZE_MEAN;

    private PopulationSize _popSize;
    private Double _oldParam;

    public ChangePopSizeParam(PopulationSize ps) {
        this._popSize = ps;
    }

    @Override
    public double propose() {
        _oldParam = _popSize.getGammaMean();
        double newParam = Math.abs(_oldParam.doubleValue() + (Randomizer.getRandomDouble() - 0.5) * _windowSize);
        _popSize.setGammaMean(newParam);
        return 0;
    }

    @Override
    public void undo() {
        _popSize.setGammaMean(_oldParam);
    }

    @Override
    public void optimize(double logAlpha) {
        _windowSize *= Math.exp(Utils.calcDelta(this, logAlpha));
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
        return "Change-PopSize-Prior-Param";
    }

    public static void main(String[] args) {
        {
            Utils._CONST_POP_SIZE = true;
            Utils._ESTIMATE_POP_SIZE = true;
            UltrametricNetwork net = new UltrametricNetwork("(O:0.5,((((A:0.5,(C:0.5,(G:0.424733668924115)I10#H3:0.27526633107588505::0.95)I6:1.2148953368828783)I5:2.0902392667517367,(((R:0.23468712)I8#H2:0.57719209683010764::0.76)I2#H1:0.4723366810412248::0.29,(Q:0.5,I10#H3:0.24368001872495645::0.05)I4:1.4561500003545642)I9:2.081453711192161)I3:0.8324858688602822,(L:0.30183286,I8#H2:0.8699237693946538::0.24)I7:0.7369103845180647)I1:1.427856628628852,I2#H1:0.8195085232672825::0.71)I0:3.0);");
            int runs = 10000;
            int counter = 0;
            for(int i = 0; i < runs; i++) {
                Operator op = new ChangePopSizeParam(new PopulationSize());
                double logHR = op.propose();
                if(logHR != Utils.INVALID_MOVE) counter++;
            }
            System.out.println(counter == runs);
            System.out.printf("%d out of %d\n", counter, runs);
            System.out.println(net.getNetwork().getRoot().getRootPopSize());
        }
        {
            UltrametricNetwork net = new UltrametricNetwork("(O:0.5,((((A:0.5,(C:0.5,(G:0.424733668924115)I10#H3:0.27526633107588505::0.95)I6:1.2148953368828783)I5:2.0902392667517367,(((R:0.23468712)I8#H2:0.57719209683010764::0.76)I2#H1:0.4723366810412248::0.29,(Q:0.5,I10#H3:0.24368001872495645::0.05)I4:1.4561500003545642)I9:2.081453711192161)I3:0.8324858688602822,(L:0.30183286,I8#H2:0.8699237693946538::0.24)I7:0.7369103845180647)I1:1.427856628628852,I2#H1:0.8195085232672825::0.71)I0:3.0);");
            int runs = 10000;
            int test = 0;
            for(int i = 0; i < runs; i++) {
                Operator op = new ChangePopSizeParam(new PopulationSize());
                double logHR = op.propose();
                op.undo();
                if(Math.abs(net.getNetwork().getRoot().getRootPopSize() - Utils._POP_SIZE_MEAN) > 0.000001) {
                    test--;
                    System.err.println(net.getNetwork().getRoot().getRootPopSize());
                }
                test++;
            }
            System.out.println(test == runs);
            System.out.printf("%d out of %d\n", test, runs);
        }
    }
}
