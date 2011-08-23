package edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;

import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReadException;
import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReader;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.RichNewickReaderAST;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.RichNewickReaderAST_Test;
import org.junit.Assert;
import org.junit.Test;

import java.io.ByteArrayInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/1/11
 * Time: 5:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickReaderAST_ANTLR_Test extends RichNewickReaderAST_Test
{

    @Override
    protected RichNewickReaderAST makeReader() {
        return new RichNewickReaderAST_ANTLR();
    }
}
