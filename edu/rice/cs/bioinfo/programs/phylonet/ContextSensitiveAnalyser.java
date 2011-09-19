package edu.rice.cs.bioinfo.programs.phylonet;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.csa.ASTContextAnalyser;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.csa.CSAError;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.commands.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:07 PM
 * To change this template use File | Settings | File Templates.
 */
class ContextSensitiveAnalyser {

    public static void analyseNetworks(Map<String,Network> sourceIdentToNetwork, Proc3<String,Integer,Integer> errorDetected)
    {


       for(String rNewickStringIdent : sourceIdentToNetwork.keySet())
       {
           Network network = sourceIdentToNetwork.get(rNewickStringIdent);


               for(CSAError error : ASTContextAnalyser.analyse(network))
               {
                    errorDetected.execute(error.Message, -1, -1);
               }

       }
    }
}
