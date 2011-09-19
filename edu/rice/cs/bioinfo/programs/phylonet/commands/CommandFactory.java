package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SyntaxCommand;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class CommandFactory {

    public static Command make(SyntaxCommand directive)
    {
        String lowerCommandName = directive.getName().toLowerCase();

        ArrayList<Parameter> params = new ArrayList<Parameter>();
        for(Parameter p : directive.getParameters())
        {
            params.add(p);
        }

        if(lowerCommandName.equals("symmetricdifference"))
        {
            return new SymmetricDifference(directive, params);
        }
        else
        {
             throw new IllegalArgumentException(String.format("Unknown command name '%s'.", directive.getName()));
        }

    }
}
