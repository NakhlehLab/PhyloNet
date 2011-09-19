package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:26 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Command {

    public boolean checkParams(Proc3<String,Integer,Integer> errorDetected);

    public boolean checkContext(Map<String,Network> sourceIdentToNetwork, Proc3<String,Integer,Integer> errorDetected);

    public void executeCommand(Proc<String> onResult);

    public <R,T,E extends Exception> R execute(CommandAlgo<R,T,E> algo, T input) throws E;
}
