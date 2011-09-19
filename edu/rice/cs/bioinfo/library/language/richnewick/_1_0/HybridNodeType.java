package edu.rice.cs.bioinfo.library.language.richnewick._1_0;


/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/26/11
 * Time: 12:18 PM
 * To change this template use File | Settings | File Templates.
 */
public enum HybridNodeType
{
    Recombination, Hybridization, LateralGeneTransfer, Unspecified;

    public static HybridNodeType fromString(String string)
    {
        if(string.toLowerCase().equals("r"))
        {
            return HybridNodeType.Recombination;
        }
        else if(string.toLowerCase().equals("lgt"))
        {
            return HybridNodeType.LateralGeneTransfer;
        }
        else if(string.toLowerCase().equals("h"))
        {
            return HybridNodeType.Hybridization;
        }
        else
        {
            throw new IllegalArgumentException("Unknown hybrid type: " + string);
        }
    }
}
