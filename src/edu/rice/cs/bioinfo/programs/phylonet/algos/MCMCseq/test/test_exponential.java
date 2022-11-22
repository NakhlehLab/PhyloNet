package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.test;
/*
 * @ClassName:   test_exponential
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        11/4/22 9:34 PM
 */

import org.apache.commons.math3.distribution.ExponentialDistribution;

public class test_exponential {
    /* Constructor */
    public test_exponential() {
        ExponentialDistribution exp = new ExponentialDistribution(6.6);
        System.out.println(exp.getMean());
        double y = exp.density(100);
        System.out.println(y);

    }

    public static void main(String[] args) {
        test_exponential x = new test_exponential();

    }
}
