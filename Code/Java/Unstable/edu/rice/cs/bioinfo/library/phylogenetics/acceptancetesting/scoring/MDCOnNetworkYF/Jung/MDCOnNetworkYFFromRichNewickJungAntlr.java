package edu.rice.cs.bioinfo.library.phylogenetics.acceptancetesting.scoring.MDCOnNetworkYF.Jung;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;

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
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }
}
