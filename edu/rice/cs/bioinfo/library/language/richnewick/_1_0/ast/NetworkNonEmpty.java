package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/27/11
 * Time: 1:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkNonEmpty implements Network{

    public final DescendantList PrincipleDescendants;

    public final NetworkInfo PrincipleInfo;

    public final RootageQualifier RootageQualifier;

    public NetworkNonEmpty(RootageQualifier rootageQualifier, DescendantList principleDescendants, NetworkInfo principleInfo)
    {
        RootageQualifier = rootageQualifier;
        PrincipleDescendants = principleDescendants;
        PrincipleInfo = principleInfo;
    }

    public <R, T, E extends Exception> R execute(NetworkAlgo<R, T, E> algo, T input) throws E {
        return algo.forNetworkNonEmpty(this, input);
    }
}
