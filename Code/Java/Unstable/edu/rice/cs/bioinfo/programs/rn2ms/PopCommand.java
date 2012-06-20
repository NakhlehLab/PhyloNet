package edu.rice.cs.bioinfo.programs.rn2ms;
import java.math.BigDecimal;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/20/12
 * Time: 2:53 PM
 * To change this template use File | Settings | File Templates.
 */
abstract class PopCommand implements Comparable<PopCommand>
{
    public final BigDecimal BackTime;

    PopCommand(BigDecimal backTime)
    {
        BackTime = backTime;
    }

    public abstract <R,T,E extends Exception> R execute(PopCommandAlgo<R,T,E> algo, T input) throws E;
}
