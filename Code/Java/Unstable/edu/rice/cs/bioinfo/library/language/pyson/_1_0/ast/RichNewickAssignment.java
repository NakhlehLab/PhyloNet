package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 5:48 PM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickAssignment implements PySONNode{

    public final boolean IsDefault;

    public final Identifier LHSIdentifier;

    public final RichNewickString RHSRichNewickString;


    public RichNewickAssignment(boolean isDefault, Identifier lhsIdentifier, RichNewickString rhsRichNewickString)
    {
        IsDefault = isDefault;
        LHSIdentifier = lhsIdentifier;
        RHSRichNewickString = rhsRichNewickString;
    }
}
