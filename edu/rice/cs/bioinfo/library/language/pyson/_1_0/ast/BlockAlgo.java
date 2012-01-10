package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 2:55 PM
 * To change this template use File | Settings | File Templates.
 */
public interface BlockAlgo<R, T, E extends Exception> {

    public R forTreesBlock(TreesBlockBody treesBlock, T input) throws E;

    public R forNetworksBlock(NetworksBlockBody treesBlock, T input) throws E;

    public R forPhylonetBlockBody(PhyloNetBlockBody phyloBlock, T input) throws E;

    public R forDataBlock(DataBlockBody dataBlock, T input) throws E;
}
