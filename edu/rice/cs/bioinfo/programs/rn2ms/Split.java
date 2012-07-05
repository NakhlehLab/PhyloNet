package edu.rice.cs.bioinfo.programs.rn2ms;
import java.math.BigDecimal;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/20/12
 * Time: 2:54 PM
 * To change this template use File | Settings | File Templates.
 */
class Split extends PopCommand
{
    public final int OldPopulationNumber;

    public final int NewPopulationNumber;

    public final BigDecimal OldPopulationProbability;

    Split(BigDecimal backTime, int oldPopulationNumber, int newPopulationNumber, BigDecimal oldPopulationProb) {
        super(backTime);
        OldPopulationNumber = oldPopulationNumber;
        NewPopulationNumber = newPopulationNumber;
        OldPopulationProbability = oldPopulationProb;
    }


    public int compareTo(PopCommand o) {
        return o.execute(new PopCommandAlgo<Integer,Object,RuntimeException>()
        {
            public Integer forSplit(Split split, Object input) throws RuntimeException {
                int backComp = BackTime.compareTo(split.BackTime);
                if(backComp != 0)
                    return backComp;
                else
                    return new Integer(NewPopulationNumber).compareTo(split.NewPopulationNumber);
            }

            public Integer forJoin(Join join, Object input) throws RuntimeException {
               return BackTime.compareTo(join.BackTime);
            }
        }, null);
    }

    @Override
    public <R, T, E extends Exception> R execute(PopCommandAlgo<R, T, E> algo, T input) throws E {
        return algo.forSplit(this, input);
    }
}
