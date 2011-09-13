package edu.rice.bioinfo.programs.phylonet.commands;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class CommandFactory {

    public static Command make(String sourceName)
    {
        String lowerCommandName = sourceName.toLowerCase();

        if(lowerCommandName.equals("symmetricdifference"))
        {
            return new SymmetricDifference();
        }
        else
        {
             throw new IllegalArgumentException(String.format("Unknown command name '%s'.", sourceName));
        }

    }
}
