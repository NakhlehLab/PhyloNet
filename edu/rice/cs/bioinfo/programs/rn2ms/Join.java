package edu.rice.cs.bioinfo.programs.rn2ms;
import java.math.BigDecimal;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/20/12
 * Time: 2:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class Join extends PopCommand
{
    public final int PrimaryPopulationNumber;

    public final int SecondaryPopulationNumber;

    Join(BigDecimal backTime, int primaryPopulationNumber, int secondaryPopulationNumber) {
        super(backTime);
        PrimaryPopulationNumber = primaryPopulationNumber;
        SecondaryPopulationNumber = secondaryPopulationNumber;
    }

    @Override
    public <R, T, E extends Exception> R execute(PopCommandAlgo<R, T, E> algo, T input) throws E {
        return algo.forJoin(this, input);
    }

    public int compareTo(PopCommand o) {
        return this.BackTime.compareTo(o.BackTime);
    }
}
