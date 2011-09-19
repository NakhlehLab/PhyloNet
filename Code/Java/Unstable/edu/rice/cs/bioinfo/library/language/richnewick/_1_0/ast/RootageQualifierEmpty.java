package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;


/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/3/11
 * Time: 4:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class RootageQualifierEmpty implements RootageQualifier {

    public static final  RootageQualifierEmpty Singleton = new RootageQualifierEmpty();

    private RootageQualifierEmpty()
    {

    }


    public <R, T, E extends Exception> R execute(RootageQualifierAlgo<R, T, E> algo, T input) throws E {
        return algo.forEmptyQualifier(this, input);
    }
}
