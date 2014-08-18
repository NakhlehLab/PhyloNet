package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPP.prototypes;

import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPP.R;

public class TestingNewSizeOfR
{

    public static double getSize(int m, int numberOfLoci)
    {
        double total = 0;
        for (int i = 1; i <= m;i++)
        {
            total += Math.pow(R.getMatrixSizeWithoutN(i),numberOfLoci);
        }
        return Math.log(total);
    }


    public static void main(String[] args)
    {
        for (int m = 1; m <= 40; m++)
            for (int numberOfLoci = 1; numberOfLoci <= 1; numberOfLoci++)
            {
                System.out.println(m + "," + numberOfLoci +","+ getSize(m,numberOfLoci));
            }
    }
}
