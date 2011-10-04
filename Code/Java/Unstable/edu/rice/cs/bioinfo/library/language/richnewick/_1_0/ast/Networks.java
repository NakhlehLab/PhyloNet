package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

import java.nio.ReadOnlyBufferException;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/16/11
 * Time: 1:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class Networks implements AbstractSyntaxNode {

    public final Iterable<NetworkNonEmpty> Networks;

    public Networks(Iterable<NetworkNonEmpty> networks)
    {
        Networks = networks;
    }
}
