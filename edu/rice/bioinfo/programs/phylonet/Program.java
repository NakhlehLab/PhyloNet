package edu.rice.bioinfo.programs.phylonet;

import com.sun.org.apache.xerces.internal.util.SymbolTable;
import edu.rice.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.bioinfo.library.language.pyson._1_0.ast.Blocks;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SACFactoryFromAST;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SyntaxCommand;
import edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.*;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.StringsAndCommands;
import edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.PySONLexer;
import edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.PySONParser;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.Networks;
import edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickParser;
import edu.rice.bioinfo.library.programming.Proc3;
import edu.rice.bioinfo.programs.phylonet.commands.Command;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;

import javax.swing.text.html.parser.Parser;
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

    private static int _numErrorsDetected = 0;

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

        Proc3<String, Integer, Integer> errorDetected = new Proc3<String, Integer, Integer>()
        {

                // this function called if an error is detected. display error to user.
                public void execute(String message, Integer line, Integer col) {

                    System.err.println(String.format("Error at [%s,%s]: %s", line, col, message));
                    _numErrorsDetected++;

                }
        };

        Blocks blocks;
        try
        {
            // parse input nexus file to AST. Exception thrown on syntax error.
            blocks = parseToBlocks(new FileInputStream(nexusFile));
        }
        catch(CoordinateParseErrorsException e)  // syntax error detected
        {
            for(CoordinateParseError error : e.Errors)
            {
                errorDetected.execute(e.getMessage(), error.getLineNumber(), error.getColumnNumber() );
            }
            return;
        }

        // no syntax errors at this point

        StringsAndCommands sac = SACFactoryFromAST.make(blocks); // covert AST to our IR

        // perform CSA to make sure the input nexus file doesn't contain non-grammar errors
        List<AbstractMap.SimpleImmutableEntry<SyntaxCommand,Command>> commandToCommand = ContextSensitiveAnalyser.analyseCommands(sac.getCommands(), errorDetected);

        Map<String,Network> sourceIdentToNetwork = makeNetworks(sac, errorDetected);

        ContextSensitiveAnalyser.analyseNetworks(sourceIdentToNetwork, errorDetected);

        if(_numErrorsDetected == 0)
        {
            executeCommands(commandToCommand, sourceIdentToNetwork, errorDetected);
        }


        return;



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

    private static void executeCommands(List<AbstractMap.SimpleImmutableEntry<SyntaxCommand,Command>> commandToCommand,
                                        Map<String,Network> sourceIdentToNetwork, Proc3<String,Integer,Integer> onError)
    {
        for(AbstractMap.SimpleImmutableEntry<SyntaxCommand,Command> kvp : commandToCommand)
        {
            SyntaxCommand sCommand = kvp.getKey();
            Command command = kvp.getValue();
            CommandExecuter.execute(sCommand, command, sourceIdentToNetwork, System.out, onError);

        }
    }



    private static void showParseError(CoordinateParseErrorsException e) {
    }

    private static Blocks parseToBlocks(FileInputStream fileInputStream) throws IOException, CoordinateParseErrorsException {

        return edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.Parser.parse(fileInputStream);

    }

    private static void showFileDoesNotExist(File inFile) {
        //To change body of created methods use File | Settings | File Templates.
    }

    private static void showUsage() {
        //To change body of created methods use File | Settings | File Templates.
    }
}
