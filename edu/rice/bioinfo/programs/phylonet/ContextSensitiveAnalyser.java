package edu.rice.bioinfo.programs.phylonet;

import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.Parameter;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.StringsAndCommands;
import edu.rice.bioinfo.library.programming.*;
import edu.rice.bioinfo.programs.phylonet.commands.*;
import sun.text.normalizer.IntTrie;

import java.lang.reflect.Array;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:07 PM
 * To change this template use File | Settings | File Templates.
 */
class ContextSensitiveAnalyser {

    public static void Analyse(StringsAndCommands ir, Proc3<String,Integer,Integer> errorDetected)
    {
       for(edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.Command irCommand : ir.getCommands())
       {
           Command command;
           try
           {
                command = CommandFactory.make(irCommand.getName());
           }
           catch(IllegalArgumentException e)
           {
               errorDetected.execute(e.getMessage(), irCommand.getLine(), irCommand.getColumn());
               continue;
           }

           ArrayList<Parameter> params = new ArrayList<Parameter>();
           for(Parameter p : irCommand.getParameters())
           {
               params.add(p);
           }

           if(command.getNumParameters() != params.size())
           {
                errorDetected.execute(String.format("Expected %s parameters for command but found %s.", command.getNumParameters(), params.size()),
                                                    irCommand.getLine(), irCommand.getColumn());
           }
       }
    }
}
