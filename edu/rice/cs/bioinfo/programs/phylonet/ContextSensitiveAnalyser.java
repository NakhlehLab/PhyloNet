package edu.rice.cs.bioinfo.programs.phylonet;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.BlockContents;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.RichNewickAssignment;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.*;
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

    public static void analyseNetworks(Map<String,NetworkNonEmpty> sourceIdentToNetwork, BlockContents blockContents, final Proc3<String,Integer,Integer> errorDetected)
    {
       for(String rNewickStringIdent : sourceIdentToNetwork.keySet())
       {
           final NetworkNonEmpty network = sourceIdentToNetwork.get(rNewickStringIdent);

           for(CSAError error : ASTContextAnalyser.analyse(network))
           {
               errorDetected.execute(error.Message, -1, -1);
           }

           final RichNewickAssignment assigment = blockContents.getRichNewickAssigment(rNewickStringIdent);

           if(assigment.isDefinedByNetworksBlock())
           {
               network.RootageQualifier.execute(new RootageQualifierAlgo<Object, Object, RuntimeException>() {
                   public Object forEmptyQualifier(RootageQualifierEmpty rootageQualifierEmpty, Object o) throws RuntimeException {
                       return null;  //To change body of implemented methods use File | Settings | File Templates.
                   }

                   public Object forNonEmptyQualifier(RootageQualifierNonEmpty rootageQualifierNonEmpty, Object o) throws RuntimeException {

                       if(!rootageQualifierNonEmpty.isRooted())
                       {
                           errorDetected.execute("Networks may not be unrooted.", assigment.getRichNewickStringLine(), assigment.getRichNewickStringColumn());
                       }
                       return null;
                   }
               }, null);
           }


       }
    }
}
