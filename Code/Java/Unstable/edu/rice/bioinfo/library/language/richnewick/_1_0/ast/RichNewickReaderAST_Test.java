package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReadException;
import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReader;
import edu.rice.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilderNoAction;
import org.junit.Assert;
import org.junit.Test;

import java.io.ByteArrayInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/1/11
 * Time: 5:34 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class RichNewickReaderAST_Test {

    protected abstract RichNewickReaderAST makeReader();

    @Test
    public void testRead() throws Exception
    {
        RichNewickReaderAST reader = makeReader();

        try
        {
            reader.read(new ByteArrayInputStream("((H#1:::.6)A,(H#1:::.6)B)C;".getBytes()), GraphBuilderNoAction.Singleton);
            Assert.fail("Expected exception.");
        }
        catch(RichNewickReadException e)
        {

        }

           reader.read(
                    new ByteArrayInputStream("((1, ((2, (3, (4)Y#H1)g)e, (((Y#H1, 5)h, 6)f)X#H2)c)a, ((X#H2, 7)d, 8)b)r;".getBytes()),
                    GraphBuilderNoAction.Singleton);


    }
}
