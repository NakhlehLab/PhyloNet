package edu.rice.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 5:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickString implements PySONNode{

    public final String String;

    public RichNewickString(String string)
    {
        String = string;
    }
}
