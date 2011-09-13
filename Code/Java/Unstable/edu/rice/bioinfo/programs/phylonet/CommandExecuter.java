package edu.rice.bioinfo.programs.phylonet;

import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.Parameter;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SyntaxCommand;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.bioinfo.library.programming.Proc3;
import edu.rice.bioinfo.library.programming.extensions.java.lang.iterable.Proc;
import edu.rice.bioinfo.programs.phylonet.commands.Command;
import edu.rice.bioinfo.programs.phylonet.commands.CommandAlgo;
import edu.rice.bioinfo.programs.phylonet.commands.SymmetricDifference;

import java.io.PrintStream;
import java.util.Iterator;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/8/11
 * Time: 5:38 PM
 * To change this template use File | Settings | File Templates.
 */
public class CommandExecuter {

    public static void execute(final SyntaxCommand sCommand, Command command, final Map<String, Network> sourceIdentToNetwork,
                               final PrintStream defaultOutput, final Proc3<String,Integer,Integer> onError) {

        command.execute(new CommandAlgo<Object, Object, RuntimeException>() {

            public Object forSymmetricDifference(SymmetricDifference command, Object input) throws RuntimeException {

                executeSymmetricDifference(sCommand, command, sourceIdentToNetwork, defaultOutput, onError);
                return null;
            }
        }, null);

    }

    private static void executeSymmetricDifference(final SyntaxCommand sCommand, SymmetricDifference command, Map<String, Network> sourceIdentToNetwork,
                                                   PrintStream defaultOutput, final Proc3<String,Integer,Integer> onError) {

        Iterator<Parameter> params = sCommand.getParameters().iterator();

        String modelTreeIdent = params.next().getValue();
        String experimentalTreeIdent = params.next().getValue();

        Network modelTree = sourceIdentToNetwork.get(modelTreeIdent);
        Network experimentalTree = sourceIdentToNetwork.get(experimentalTreeIdent);

        boolean noError = command.checkContext(modelTree, modelTreeIdent, experimentalTree, experimentalTreeIdent, new Proc<String>()
        {
            public void execute(String s) {
                onError.execute(s, sCommand.getLine(), sCommand.getColumn());
            }
        });

        if(noError)
        {
            SymmetricDifference.SymmetricDifferenceResult result = command.execute(modelTree, experimentalTree);
            defaultOutput.println(result.getFalseNegativesCount() + " " + result.getFalsePositivesCount() + " " +
                                  result.getModelInternalEdgesCount() + " " + result.getExperamentalInternalEdgesCount());
        }







    }
}
