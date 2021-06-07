package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/3/12
 * Time: 10:54 AM
 * To change this template use File | Settings | File Templates.
 */
public class IsRooted implements RootageQualifierAlgo<Boolean, Object, RuntimeException> {
    public Boolean forEmptyQualifier(RootageQualifierEmpty rootage, Object input) throws RuntimeException {
       return true;
    }

    public Boolean forNonEmptyQualifier(RootageQualifierNonEmpty rootage, Object input) throws RuntimeException {
        return rootage.isRooted();
    }
}
