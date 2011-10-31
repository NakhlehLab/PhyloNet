package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.util.ArrayList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class CommandFactory {

    public static Command make(SyntaxCommand directive, Map<String,NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected)
    {
        String lowerCommandName = directive.getName().toLowerCase();

        ArrayList<Parameter> params = new ArrayList<Parameter>();
        for(Parameter p : directive.getParameters())
        {
            params.add(p);
        }

        if(lowerCommandName.equals("symmetricdifference") || lowerCommandName.equals("rf"))
        {
            return new SymmetricDifference(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("lca"))
        {
            return new LCA(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("mast"))
        {
            return new MAST(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("riatahgt"))
        {
            return new RIATAHGT(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("deepcoalcount"))
        {
            return new DeepCoalCount(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("processgt"))
        {
            return new ProcessGT(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("charnet"))
        {
            return new CharNet(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else
        {
             throw new IllegalArgumentException(String.format("Unknown command name '%s'.", directive.getName()));
        }

    }
}
