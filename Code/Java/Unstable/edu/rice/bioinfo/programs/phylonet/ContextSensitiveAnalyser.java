package edu.rice.bioinfo.programs.phylonet;

import edu.rice.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.Parameter;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.StringsAndCommands;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SyntaxCommand;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.NetworkInfo;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.Networks;
import edu.rice.bioinfo.library.language.richnewick._1_0.csa.ASTContextAnalyser;
import edu.rice.bioinfo.library.language.richnewick._1_0.csa.CSAError;
import edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickParser;
import edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickReaderAST_ANTLR;
import edu.rice.bioinfo.library.programming.*;
import edu.rice.bioinfo.programs.phylonet.commands.*;
import org.antlr.runtime.RecognitionException;
import sun.text.normalizer.IntTrie;

import java.awt.image.ImagingOpException;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:07 PM
 * To change this template use File | Settings | File Templates.
 */
class ContextSensitiveAnalyser {

    static List<AbstractMap.SimpleImmutableEntry<SyntaxCommand,Command>> analyseCommands(Iterable<SyntaxCommand> commands,
                                                      final Proc3<String,Integer,Integer> errorDetected)
    {
       LinkedList<AbstractMap.SimpleImmutableEntry<SyntaxCommand,Command>> tbr = new LinkedList<AbstractMap.SimpleImmutableEntry<SyntaxCommand,Command>>();
       for(final SyntaxCommand syntaxCommand : commands)
       {
           Command command = CommandFactory.make(syntaxCommand.getName());
           tbr.add(new AbstractMap.SimpleImmutableEntry<SyntaxCommand, Command>(syntaxCommand, command));

           final ArrayList<Parameter> params = new ArrayList<Parameter>();
           for(Parameter p : syntaxCommand.getParameters())
           {
               params.add(p);
           }

           command.checkParams(syntaxCommand, params, errorDetected);

       }

        return tbr;
    }



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
