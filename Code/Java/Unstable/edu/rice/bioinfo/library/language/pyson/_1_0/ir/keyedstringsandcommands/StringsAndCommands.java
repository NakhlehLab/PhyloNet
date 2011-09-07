package edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands;

import javax.print.DocFlavor;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 1:35 PM
 * To change this template use File | Settings | File Templates.
 */
public interface StringsAndCommands {

    String getRickNewickString(String identifier);

    List<Command> getCommands();
}
