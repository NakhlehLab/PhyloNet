package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 1:35 PM
 * To change this template use File | Settings | File Templates.
 */
public interface StringsAndCommands {

    String getRickNewickString(String identifier);

    Iterable<String> getRickNewickStringIdentifiers();

    Iterable<SyntaxCommand> getCommands();


}
