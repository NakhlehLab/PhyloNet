package edu.rice.bioinfo.programs.phylonet.commands;

import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.Parameter;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SyntaxCommand;
import edu.rice.bioinfo.library.programming.Proc3;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:26 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Command {

    public void checkParams(SyntaxCommand command, ArrayList<Parameter> params, Proc3<String,Integer,Integer> errorDetected);

    public <R,T,E extends Exception> R execute(CommandAlgo<R,T,E> algo, T input) throws E;
}
