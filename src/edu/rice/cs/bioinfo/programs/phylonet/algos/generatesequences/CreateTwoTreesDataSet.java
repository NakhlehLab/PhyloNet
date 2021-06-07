package edu.rice.cs.bioinfo.programs.phylonet.algos.generatesequences;

import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.ConfigurationBuilder;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

public class CreateTwoTreesDataSet
{
    public static String getModifiedNetwork(double a, double b, double gamma)
    {
        String networkString = "(((B:BBOT):0.0)ANC#H1:BTOP::GAMMA,(ANC#H1:0.01::NOTGAMMA,A:ABOT):ATOP)N;";

        String actualNetwork = networkString
                .replace("BBOT",""+a)
                .replace("BTOP", "" + b)
                .replace("NOTGAMMA",""+(1-gamma))
                .replace("GAMMA", "" + gamma)
                .replace("ABOT", "" + (a + .01))
                .replace("ATOP",""+(b-.01));

        System.out.println(actualNetwork);

        return  actualNetwork;
    }

    public static void main(String[] args)
    {
        double[] aValues = {.1,1,3,5};
        double[] bValues = {.1,1,3,5};

        double[] sValues = {.01,.05,.1,.5};
        double[] gammas = {.1,.3,.5};

        Map<String,Integer> f = new HashMap<>();
        f.put("A",1);
        f.put("B",1);
        for (double a : aValues)
        {
            for (double b : bValues)
                for (double s : sValues)
                    for (double g : gammas)
                    {
                        Path output = Paths.get("foo",""+a,""+b,""+s,""+g);
                        String networkString = getModifiedNetwork(a,b,g);
                        CreateDatasetInput i = new CreateDatasetInput();
                        i.iterations = 1000;
                        i.threads = 3;
                        i.numberOfLoci = 10000;
                        i.mutationRate = s;
                        i.network = networkString;
                        i.networkAlleleCounts = f;

                        CreateDataset.createDataset(output,i);
                    }
        }

    }
}
