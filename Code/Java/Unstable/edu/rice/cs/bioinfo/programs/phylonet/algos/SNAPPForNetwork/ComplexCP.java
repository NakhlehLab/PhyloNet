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

    public void setToProductOf(double a_re, double a_im, double b_re, double b_im) {
        _re = b_re * a_re - b_im * a_im;
        _im = b_re * a_im + b_im * a_re;
    }

    public void addProductOf(double a_re, double a_im, double b_re, double b_im) {
        _re += b_re * a_re - b_im * a_im;
        _im += b_re * a_im + b_im * a_re;
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

    public void setToQuotientOf(double num_re, double num_im, double div_re, double div_im) {
        double f = div_re * div_re + div_im * div_im;
        double re = (num_re * div_re + num_im * div_im) / f;
        double im = (num_im * div_re - num_re * div_im) / f;
        _re = re;
        _im = im;
    }

    public String toString() {
        return "(" +_re + "," + _im + ")";
    }



}
