package edu.rice.cs.bioinfo.library.language.vaal._1_1.out.reading.parsers;

import edu.rice.cs.bioinfo.library.language.vaal._1_1.out.reading.ir.LineDifferenceRecord;
import edu.rice.cs.bioinfo.library.language.vaal.out.reading.OutReadResult;
import junit.framework.Assert;
import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.io.InputStream;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 3/21/13
 * Time: 2:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class VAALOutStringSplitParserTest
{
    @Test
    public void testParse() throws Exception
    {
        String input = "header1\nheader2\n0 20895 left=TGGATCAGACCTTTAGCAGC sample=C ref=G right=TGACGGTCCACGATCGGTTG";
        InputStream inStream = new ByteArrayInputStream(input.getBytes());

        OutReadResult<List<LineDifferenceRecord>> result = new VAALOutStringSplitParser().parse(inStream);
        List<LineDifferenceRecord> differences = result.getDifferences();

        Assert.assertTrue(1 == differences.size());
        Assert.assertEquals("C", differences.get(0).Sample);
        Assert.assertEquals("G", differences.get(0).Ref);
    }
}
