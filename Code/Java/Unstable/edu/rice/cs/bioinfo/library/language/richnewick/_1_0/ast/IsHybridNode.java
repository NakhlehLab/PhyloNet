package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/24/11
 * Time: 4:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class IsHybridNode implements HybridNodeQualifierAlgo<Boolean,Object, RuntimeException>
{

    public Boolean forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Object input) throws RuntimeException {
       return false;
    }

    public Boolean forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Object input) throws RuntimeException {
       return true;
    }

    public Boolean forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Object input) throws RuntimeException {
       return true;
    }
}
