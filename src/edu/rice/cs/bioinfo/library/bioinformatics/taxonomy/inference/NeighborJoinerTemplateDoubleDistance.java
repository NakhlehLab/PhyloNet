package edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/29/13
 * Time: 7:05 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NeighborJoinerTemplateDoubleDistance<N,E,G> extends NeighborJoinerTemplate<N,E,G,Double>
{
    @Override
    protected Double makeD(int value) {
        return new Double(value);
    }

    @Override
    protected Double add(Double d1, Double d2) {
        return d1 + d2;
    }

    @Override
    protected Double subtract(Double minuend, Double subtrahend) {
        return minuend - subtrahend;
    }

    @Override
    protected Double multiply(Double d1, Double d2) {
        return d1 * d2;
    }

    @Override
    protected Double divide(Double dividend, Double divisor) {
        return dividend / divisor;
    }
}
