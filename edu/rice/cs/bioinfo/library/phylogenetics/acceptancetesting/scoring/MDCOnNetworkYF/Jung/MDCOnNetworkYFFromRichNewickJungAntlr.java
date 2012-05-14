package edu.rice.cs.bioinfo.library.phylogenetics.acceptancetesting.scoring.MDCOnNetworkYF.Jung;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickReaderAST_ANTLR;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.MDCOnNetworkYFFromRichNewick;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/14/12
 * Time: 2:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class MDCOnNetworkYFFromRichNewickJungAntlr extends MDCOnNetworkYFFromRichNewickJung
{

    public MDCOnNetworkYFFromRichNewickJungAntlr()
    {
        super(new RichNewickReaderAST_ANTLR());
    }
}
