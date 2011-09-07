package edu.rice.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 5:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class Blocks implements PySONNode{

    public final Iterable<Block> Contents;

    public Blocks(Iterable<Block> blocks)
    {
        Contents = blocks;
    }
}
