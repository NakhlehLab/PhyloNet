package edu.rice.cs.bioinfo.programs.phylonet.commands;


import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.util.ArrayList;
import java.util.Map;

public abstract class CommandBaseFileOutWithDNA extends CommandBaseFileOut
{
    private Map<String,String> sourceIdentToDNA;

    public CommandBaseFileOutWithDNA(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Map<String,String> sourceIdentToDNA, Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader)
    {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
        this.sourceIdentToDNA = sourceIdentToDNA;
    }

    protected Map<String,String> getSourceIdentToDNA()
    {
        return sourceIdentToDNA;
    }
}
