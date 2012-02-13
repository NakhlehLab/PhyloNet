package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.io.IOException;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:26 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Command {

    public SyntaxCommand getDefiningSyntaxCommand();

    public boolean checkParams();

    public void executeCommand(Proc<String> displayResult) throws IOException;

    public void addSTTreeGeneratedListener(Proc<String> listener);

}
