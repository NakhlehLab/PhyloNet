package edu.rice.cs.bioinfo.library.phylogenetics.acceptancetesting.search.hillclimbing.network.rea.jung;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea.acceptancetesting.jung.ReaSearchFromRichNewickJung;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 5:27 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReaSearchFromRichNewickJungAntlr extends ReaSearchFromRichNewickJung
{
    public ReaSearchFromRichNewickJungAntlr() {
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }
}
