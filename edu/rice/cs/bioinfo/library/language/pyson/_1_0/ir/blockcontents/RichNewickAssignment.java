package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/20/11
 * Time: 4:27 PM
 * To change this template use File | Settings | File Templates.
 */
public interface RichNewickAssignment {

    String getIdentifier();

    boolean isDefinedByNetworksBlock();

    boolean isDefinedByTreesBlock();

    String getRichNewickString();

    int getRichNewickStringLine();

    int getRichNewickStringColumn();
}
