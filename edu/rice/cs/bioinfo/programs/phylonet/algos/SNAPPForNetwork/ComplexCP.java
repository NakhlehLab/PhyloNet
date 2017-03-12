package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;


/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 3/11/17
 * Time: 3:10 PM
 * To change this template use File | Settings | File Templates.
 */

//a compact implementation of complex number - for performance
public class ComplexCP {
    public double _re;
    public double _im;

    public ComplexCP() {
        _re = 0.0;
        _im = 0.0;
    }

    public ComplexCP(double re) {
        _re = re;
        _im = 0.0;
    }

    public ComplexCP(double re, double im) {
        _re = re;
        _im = im;
    }

    public ComplexCP(ComplexCP c) {
        _im = c._im;
        _re = c._re;
    }

    public void setToProductOf(ComplexCP a, ComplexCP b) {
        double re = b._re * a._re - b._im * a._im;
        double im = b._re * a._im + b._im * a._re;
        _re = re;
        _im = im;
    }

    public void addProductOf(ComplexCP a, ComplexCP b) {
        double re = b._re * a._re - b._im * a._im;
        double im = b._re * a._im + b._im * a._re;
        _re += re;
        _im += im;
    }

    public void add(double re) {
        _re += re;
    }

    public void setToSumOfScaled(ComplexCP a, double x, ComplexCP b, double y) {
        double re = a._re * x + b._re * y;
        double im = a._im * x + b._im * y;
        _re = re;
        _im = im;
    }

    public void setToQuotientOf(ComplexCP num, ComplexCP div) {
        double f = div._re * div._re + div._im * div._im;
        double re = (num._re * div._re + num._im * div._im) / f;
        double im = (num._im * div._re - num._re * div._im) / f;
        _re = re;
        _im = im;
    }

    public String toString() {
        return "(" +_re + "," + _im + ")";
    }



}
