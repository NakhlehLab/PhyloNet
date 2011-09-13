package edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 1:36 PM
 * To change this template use File | Settings | File Templates.
 */
public interface SyntaxCommand {

    public int getLine();

    public int getColumn();

    public String getName();

    public Iterable<Parameter> getParameters();


}
