package edu.rice.cs.bioinfo.programs.phylonet;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.Blocks;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SACFactoryFromAST;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.StringsAndCommands;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickParser;
import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.commands.Command;
import edu.rice.cs.bioinfo.programs.phylonet.commands.CommandFactory;
import org.omg.Dynamic.Parameter;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 1:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program {

    // tracks whether the program should execute commands from script during execution phase.
    // flip to false on syntax or context errors.
    private static boolean _allowCommandExecution = true;

    // this function called if an error is detected. display error to user.
    private static final  Proc3<String, Integer, Integer> _errorDetected = new Proc3<String, Integer, Integer>()
    {
        public void execute(String message, Integer line, Integer col) {

            int oneBasedColumn = col + 1;
            System.err.println(String.format("Error at [%s,%s]: %s", line, oneBasedColumn, message));
            _allowCommandExecution = false;
        }
    };

    private static final Proc<String> _displayResult = new Proc<String>()
    {
        public void execute(String s) {

            System.out.println(s);
        }
    };

    public static void main(String[] args) throws FileNotFoundException, IOException
    {

        if(args.length != 1) // we only expect one parameter, the input nexus file
        {
            showUsage();
            return;
        }

        File nexusFile = new File(args[0]);

        if(!nexusFile.exists()) // assert the file the user gave us actually exists
        {
            showFileDoesNotExist(nexusFile);
            return;
        }

        /*
         * Parse input nexus file to AST. Report and terminate on syntax errors.
         */
        Blocks blocks;
        try
        {
            blocks = parseToBlocks(new FileInputStream(nexusFile)); // Exception thrown on syntax errors.
        }
        catch(CoordinateParseErrorsException e)  // syntax errors detected
        {
            for(CoordinateParseError error : e.Errors)
            {
                _errorDetected.execute(error.getMessage(), error.getLineNumber(), error.getColumnNumber() );
            }
            return;
        }

        // no syntax errors at this point

        StringsAndCommands sac = SACFactoryFromAST.make(blocks); // covert AST to our IR

        /*
         * Perform Context Sensitive Analysis
         */
        Map<String,Network> sourceIdentToNetwork = makeNetworks(sac, _errorDetected); // make a Network representation of each defined Rich Newick string
                                                                                        // map is keyed by source code identifier of the string

        ContextSensitiveAnalyser.analyseNetworks(sourceIdentToNetwork, _errorDetected); // check the Networks for any context errors

        LinkedList<Command> commands = new LinkedList<Command>();

        for(SyntaxCommand sCommand : sac.getCommands())
        {
            commands.add(CommandFactory.make(sCommand));
        }

        for(Command command : commands)
        {
            command.checkParams(_errorDetected);
            command.checkContext(sourceIdentToNetwork, _errorDetected);
        }

        /*
         * Execute commands if no errors detected
         */

        if(_allowCommandExecution)
        {
           for(Command command : commands)
           {
               try
               {
                    showCommand(command.getDefiningSyntaxCommand());
                    command.executeCommand(_displayResult);
               }
               catch(IOException e)
               {
                   SyntaxCommand motivatingSyntaxCommand = command.getDefiningSyntaxCommand();
                   _errorDetected.execute(String.format("Error executing command '%s' (%s).", motivatingSyntaxCommand.getName(), e.getMessage()),
                                          motivatingSyntaxCommand.getLine(), motivatingSyntaxCommand.getColumn());
               }
           }
        }


        return;



    }

    private static void showCommand(SyntaxCommand definingSyntaxCommand) {

        StringBuffer accum = new StringBuffer("\n" + definingSyntaxCommand.getName());

        for(edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.Parameter p : definingSyntaxCommand.getParameters())
        {
              accum.append(" " + p.getValue());
        }

        _displayResult.execute(accum.toString());
    }


    private static Map<String, Network> makeNetworks(StringsAndCommands sac, Proc3<String, Integer, Integer> errorDetected) throws IOException {

        HashMap<String,Network> tbr = new HashMap<String, Network>();
        for(String richNewickSourceIdent : sac.getRickNewickStringIdentifiers())
        {
            String richNewickString = sac.getRickNewickString(richNewickSourceIdent);


            try
            {
                Iterator<Network> oneNetwork = RichNewickParser.parse(
                        new ByteArrayInputStream(richNewickString.getBytes())).Networks.iterator();

                if(!oneNetwork.hasNext())
                {
                    errorDetected.execute(
                            String.format("Rich Newick string %s does not define a network.", richNewickSourceIdent),
                            -1, -1);
                }

                Network network = oneNetwork.next();

                if(oneNetwork.hasNext())
                {
                     errorDetected.execute(
                            String.format("Rich Newick string %s does defines more than one network.", richNewickSourceIdent),
                            -1, -1);
                    continue;
                }

                tbr.put(richNewickSourceIdent, network);
            }
            catch(CoordinateParseErrorsException e)
            {
                for(CoordinateParseError error : e.Errors)
                {
                    errorDetected.execute(error.getMessage(), -1, -1);
                }
            }
        }

        return tbr;
    }

    private static Blocks parseToBlocks(FileInputStream fileInputStream) throws IOException, CoordinateParseErrorsException {

        return edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.Parser.parse(fileInputStream);

    }

    private static void showFileDoesNotExist(File inFile) {
        //To change body of created methods use File | Settings | File Templates.
    }

    private static void showUsage() {

        System.out.println("\nUsage: java -jar [phylonet.jar] [nexus file]");

    }
}
