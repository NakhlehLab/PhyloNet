package edu.rice.cs.bioinfo.programs.networksearchgen;

import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.math.BigDecimal;
import java.util.LinkedList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/21/12
 * Time: 2:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program
{
    public static void main(String[] args)
    {
        if(args[0].equals("ms"))
        {
            int num_gt = Integer.parseInt(args[1]);
            String networkNewick = args[2];
            BigDecimal ultrametricThreshold = new BigDecimal(args[3]);
            Tuple<String,Map<String,Integer>> toMSScriptResult = edu.rice.cs.bioinfo.programs.rn2ms.Program.toMSScript(num_gt, networkNewick, ultrametricThreshold);
            System.out.println(toMSScriptResult.Item1);
        }
        else if(args[0].equals("st"))
        {
            String[] lines = args[1].split("//");
            for(String line : lines)
            {
                line = line.trim();
                if(line.startsWith("("))
                {
                    System.out.println(line);
                }
            }
        }
        else
        {
            throw new IllegalArgumentException(args[0]);
        }

    }
}
