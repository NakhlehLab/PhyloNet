package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

import java.nio.ReadOnlyBufferException;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 4:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class Network implements AbstractSyntaxNode {

    public final DescendantList PrincipleDescendants;

    public final NetworkInfo PrincipleInfo;

    public Network(DescendantList principleDescendants, NetworkInfo principleInfo)
    {
        PrincipleDescendants = principleDescendants;
        PrincipleInfo = principleInfo;
    }
}
