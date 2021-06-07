package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPP.prototypes;


import convolutionlib.JNIConvolution;
import org.apache.commons.math3.util.ArithmeticUtils;

public class JNIConvolutionSpeedTester
{

    public static void main(String[] args)
    {
        JNIConvolution f = new JNIConvolution(4);

        for (int i = 1; i < 100; i++)
        {
            f.prepareEngine(i);


            for (int j = 0;j < 100; j++)
                f.convolute(new double[ArithmeticUtils.pow(i,4)],new double[ArithmeticUtils.pow(i,4)],i);

            long startTime = System.nanoTime();
            for (int j = 0;j < 100; j++)
                f.convolute(new double[ArithmeticUtils.pow(i,4)],new double[ArithmeticUtils.pow(i,4)],i);

            long totalTime = System.nanoTime() - startTime;

            System.out.println(i + "," + totalTime/(1e11));
        }

    }
}
