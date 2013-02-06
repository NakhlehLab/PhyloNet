package edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/5/13
 * Time: 7:02 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NeighborJoinerTemplateIntegerDistance <N,E,G> extends NeighborJoinerTemplate<N,E,G,Integer>
{
    @Override
    protected Integer makeD(int value) {
        return value;
    }

    @Override
    protected Integer add(Integer d1, Integer d2) {
        return d1 + d2;
    }

    @Override
    protected Integer subtract(Integer minuend, Integer subtrahend) {
        return minuend - subtrahend;
    }

    @Override
    protected Integer multiply(Integer d1, Integer d2) {
        return d1 * d2;
    }

    @Override
    protected Integer divide(Integer dividend, Integer divisor) {
        return dividend / divisor;
    }
}
