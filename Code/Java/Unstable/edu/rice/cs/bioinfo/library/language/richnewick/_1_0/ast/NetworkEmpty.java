package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/27/11
 * Time: 1:28 PM
 * To change this template use File | Settings | File Templates.
 */
public final class NetworkEmpty implements Network{

    public static final NetworkEmpty Singleton = new NetworkEmpty();

    private NetworkEmpty()
    {

    }

    public <R, T, E extends Exception> R execute(NetworkAlgo<R, T, E> algo, T input) throws E {
        return algo.forNetworkEmpty(this, input);
    }
}
